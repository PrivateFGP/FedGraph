#ifndef ALGO_KERNELS_EDGE_CENTRIC_SSSP_SSSP_H_
#define ALGO_KERNELS_EDGE_CENTRIC_SSSP_SSSP_H_

#include "graph.h"
#include "ss_edge_centric_algo_kernel.h"
#include "task.h"
#include "TaskUtil.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>

// template <typename Archive> 
//     void serialize(Archive &ar, CipherEntry& ce, const unsigned int version) {
//     ar & ce.data;
//     ar & ce.nac;
//     ar & ce.mac;
// }

#define INF(type) (std::numeric_limits<type>::max()/2)
#define INV_VID -1uL

/*
 * Graph types definitions.
 */
template<typename EdgeWeightType = uint32_t>
struct SSSPData {
//    union {
//        struct {
            EdgeWeightType distance;
            GraphGASLite::VertexIdx predecessor;
//        };
        CipherEntry enc;
//    };
    bool isEncrypted;

    GraphGASLite::IterCount activeIter;

    SSSPData(const GraphGASLite::VertexIdx&)
        : distance(INF(EdgeWeightType)), predecessor(INV_VID), activeIter(-1), isEncrypted(false)
    {
        // Nothing else to do.
    }

    void intoShareVec(ShareVec& sv0, ShareVec& sv1) {
        int64_t distance_share0, distance_share1;
        int64_t predecessor_share0, predecessor_share1;
        CryptoUtil::intoShares((uint32_t)distance, distance_share0, distance_share1);
        CryptoUtil::intoShares((uint64_t)predecessor, predecessor_share0, predecessor_share1);
        // if (distance == 0) {
        //     std::cout<<distance<<" "<<distance_share0<<" "<<distance_share1<<std::endl;
        //     std::cout<<predecessor<<" "<<predecessor_share0<<" "<<predecessor_share1<<std::endl;
        // }
        sv0.shares.push_back(distance_share0);
        sv0.shares.push_back(predecessor_share0);
        sv1.shares.push_back(distance_share1);
        sv1.shares.push_back(predecessor_share1);
    }

    void fromShareVec(const ShareVec& sv0, const ShareVec& sv1) {
        distance = CryptoUtil::mergeShareAsUint32(sv0.shares[0], sv1.shares[0]);
        predecessor = CryptoUtil::mergeShareAsUint64(sv0.shares[1], sv1.shares[1]);
    }
};

template<typename EdgeWeightType = uint32_t>
struct SSSPUpdate {
    EdgeWeightType distance;
    GraphGASLite::VertexIdx predecessor;
    CipherEntry enc;
    bool isEncrypted;

    SSSPUpdate(const EdgeWeightType distance_ = INF(EdgeWeightType), const GraphGASLite::VertexIdx& predecessor_ = INV_VID)
        : distance(distance_), predecessor(predecessor_), isEncrypted(false)
    {
        // Nothing else to do.
    }

    SSSPUpdate(const CipherEntry& c)
        : enc(c), isEncrypted(true)
    {
        // Nothing else to do.
    }

    SSSPUpdate& operator+=(const SSSPUpdate& update) {
        // Min of distance.
        if (distance > update.distance) {
            distance = update.distance;
            predecessor = update.predecessor;
        }
        return *this;
    }

    friend class boost::serialization::access;
    
    template <typename Archive> 
    void serialize(Archive &ar, const unsigned int version) 
    {
        ar & distance;
        ar & predecessor;
        ar & enc;
        ar & isEncrypted;
    } 

    // 'final == true' writes to the vertex, else to the vertex update queue.
    static Task genTask(const SSSPUpdate& update1, const SSSPUpdate& update2, uint64_t vertexIndex, uint64_t srcTid, uint64_t dstTid, bool isFinal=false, bool isDummy=false) {
        MinWithParent<uint32_t, uint64_t>* mwp = new MinWithParent<uint32_t, uint64_t>();

        if (!update1.isEncrypted) {
            mwp->operands[0].plain_val = update1.distance;
            mwp->operands[0].plain_parent = update1.predecessor;
        } else {
            mwp->operands[0].enc = update1.enc;
            if (mwp->operands[0].enc.plainNum != 2) {
                printf("Error mwp->operands[0].enc.plainNum %d %d\n", mwp->operands[0].enc.plainNum);
                exit(0);
            }
        }     

        if (!update2.isEncrypted) {
            mwp->operands[1].plain_val = update2.distance;
            mwp->operands[1].plain_parent = update2.predecessor;
        } else {
            mwp->operands[1].enc = update2.enc;
            if (mwp->operands[1].enc.plainNum != 2) {
                printf("Error mwp->operands[1].enc.plainNum %d %d\n", mwp->operands[1].enc.plainNum);
                exit(0);
            }
        }
        
        mwp->isOperandEncrypted[0] = update1.isEncrypted;
        mwp->isOperandEncrypted[1] = update2.isEncrypted;
        
        struct Task task = {MIN_UINT_WITH_PARENT, vertexIndex, false, (void*)mwp, (uint64_t)-1, dstTid, srcTid, 0, isFinal, isDummy};
        return task;
    }

    // static Task genEmptyTask(const SSSPUpdate& update1, const SSSPUpdate& update2, uint64_t vertexIndex, uint64_t srcTid, uint64_t dstTid, bool isFinal=false) {
    //     MinWithParent<uint32_t, uint64_t>* mwp = new MinWithParent<uint32_t, uint64_t>();
    //     struct Task task = {MIN_UINT_WITH_PARENT, -1, false, (void*)mwp, -1, -1, -1, 0, isFinal};
    //     return task;
    // }

    static SSSPUpdate<EdgeWeightType> getTaskResult(const struct Task& task) {
        if (!task.finished) {
            throw ResultNotReadyException("Task result not ready yet!\n");
        }
        MinWithParent<uint32_t, uint64_t>* mwp = (MinWithParent<uint32_t, uint64_t>*)task.buf;

        if (!mwp->isOperandEncrypted[0]) {
            SSSPUpdate<EdgeWeightType> result = SSSPUpdate<EdgeWeightType>(mwp->operands[0].plain_val, GraphGASLite::VertexIdx(mwp->operands[0].plain_parent));
            return result;
        } else {
            SSSPUpdate<EdgeWeightType> result = SSSPUpdate<EdgeWeightType>(mwp->operands[0].enc);
            if (result.enc.plainNum != 2) {
                printf("Error results update.enc.plainNum %d %d\n", result.enc.plainNum, mwp->operands[0].enc.plainNum);
                exit(0);
            }
            return result;
        }
    }

    static void getTaskResult(const struct Task& task, ShareVec& sv) {
        if (!task.finished) {
            throw ResultNotReadyException("Task result not ready yet!\n");
        }
        MinWithParent<uint32_t, uint64_t>* mwp = (MinWithParent<uint32_t, uint64_t>*)task.buf;

        if (!mwp->isOperandEncrypted[0]) {
            printf("Pre-merging results ought to be shared!\n");
            exit(-1);
        } else {
            if (mwp->operands[0].enc.plainNum != 2) {
                printf("Error results update.enc.plainNum %d\n", mwp->operands[0].enc.plainNum);
                exit(0);
            }
            sv.shares.clear();
            sv.shares.push_back(mwp->operands[0].enc.share[0]);
            sv.shares.push_back(mwp->operands[0].enc.share[1]);
        }
    }

    static SSSPUpdate<EdgeWeightType> splitRandomShareFromUpdate(const uint32_t srcTid, const uint32_t dstTid, SSSPUpdate<EdgeWeightType>& update) {
        // The update variable must be encrypted by the public key of dstTid
        SSSPUpdate<EdgeWeightType> update_share;
        CryptoUtil& cryptoUtil = CryptoUtil::getInstance();
        if (!update.isEncrypted) {
            std::vector<int64_t> plains = {(uint64_t)update.distance, (uint64_t)update.predecessor};
            // std::cout<<"plain transfer "<<plains[0]<<" "<<plains[1]<<std::endl;
            update.enc = cryptoUtil.encryptToCipherEntry(plains, srcTid);
            update.isEncrypted = true;
            update.distance = 0;
            update.predecessor = 0;
            update_share.isEncrypted = false;
            update_share.distance = 0;
            update_share.predecessor = 0;
        } else {
            std::vector<int64_t> new_shares;
            cryptoUtil.splitRandomShareFromCipherEntry(update.enc, 2, new_shares);
            update_share.enc = update.enc; // Legal?
            update_share.isEncrypted = true;
            update.enc = cryptoUtil.encryptToCipherEntry(new_shares, srcTid);
            update.isEncrypted = true;
        }

        return update_share;
    }

    static void mergeUpdateShare(const uint32_t dstTid, SSSPUpdate<EdgeWeightType>& update, SSSPUpdate<EdgeWeightType>& update_share) {
        CryptoUtil& cryptoUtil = CryptoUtil::getInstance();
        if (update_share.isEncrypted) {
            std::vector<int64_t> shares;
            cryptoUtil.decryptCipherEntryToInt64Vec(update_share.enc, shares);
            // std::cout<<"++++++"<<std::endl;
            // std::cout<<dstTid<<std::endl;
            // std::cout<<update.enc.tid<<std::endl;
            cryptoUtil.mergeShareIntoCipherEntry(update.enc, shares, update.enc.tid);
            if (update.enc.plainNum != 2) {
                printf("Error update.enc.plainNum %d\n", update.enc.plainNum);
                exit(0);
            }
        }
        return;
    }
};

void print(const SSSPUpdate<uint32_t>& update) {
    // printf("SSSPUpdate---------\n");
    // printf("count %d\n", update.count.cnt());
    // printf(update.isEncrypted? "True":"False");
    // printf("---------\n");
}


/*
 * Algorithm kernel definition.
 */
template<typename GraphTileType>
class SSSPEdgeCentricAlgoKernel : public GraphGASLite::SSEdgeCentricAlgoKernel<GraphTileType> {
public:
    static Ptr<SSSPEdgeCentricAlgoKernel> instanceNew(const string& name,
            const GraphGASLite::VertexIdx& src) {
        return Ptr<SSSPEdgeCentricAlgoKernel>(new SSSPEdgeCentricAlgoKernel(name, src));
    }

protected:
    typedef typename GraphTileType::UpdateType UpdateType;
    typedef typename GraphTileType::VertexType VertexType;
    typedef typename GraphTileType::EdgeType::WeightType EdgeWeightType;

    std::pair<UpdateType, bool> scatter(const GraphGASLite::IterCount& iter, Ptr<VertexType>& src, EdgeWeightType& weight) const {
        auto& data = src->data();
        if (data.activeIter == iter) {
            return std::make_pair(SSSPUpdate<EdgeWeightType>(data.distance + weight, src->vid()), true);
        } else {
            return std::make_pair(SSSPUpdate<EdgeWeightType>(), false);
        }
    }

    struct Task genScatterTask(const GraphGASLite::IterCount& iter, Ptr<VertexType>& src, EdgeWeightType& weight, const GraphGASLite::VertexIdx& dstId) const {
        auto& data = src->data();
        AddWithReplaceParent<uint32_t, uint64_t>* awrp = new AddWithReplaceParent<uint32_t, uint64_t>();

        if (!data.isEncrypted) {
            awrp->operands[0].plain_val = (uint32_t)data.distance;
            awrp->operands[0].plain_parent = (uint64_t)data.predecessor;
        } else {
            awrp->operands[0].enc = data.enc;
        }
        awrp->operands[1].plain_val = (uint32_t)weight;
        awrp->operands[1].plain_parent = (uint64_t)(src->vid());

        awrp->isOperandEncrypted[0] = data.isEncrypted;
        awrp->isOperandEncrypted[1] = false;

        struct Task task = {ADD_UINT_WITH_REPLACE_PARENT, (uint64_t)dstId, false, (void*)awrp, (uint64_t)(src->vid()), this->getVertexTid(dstId), this->getVertexTid(src->vid())};
        // struct Task task = {ADD_UINT, (uint64_t)dstId, false, (void*)add};
        return task;
    }

    struct Task genScatterTask(Ptr<VertexType> src, const ShareVec& vertexData, EdgeWeightType weight, uint64_t srcId, uint64_t srcTid, uint64_t dstId, uint64_t dstTid, bool isDummy=false) const {
        AddWithReplaceParent<uint32_t, uint64_t>* awrp = new AddWithReplaceParent<uint32_t, uint64_t>();
        if (vertexData.shares.size() != 2) {
            printf("Wrong vertex data shares length for SSSP Scatter");
        }

        awrp->setCipherEntryPlainNum(0, 2);
        awrp->writeShareToOperand(vertexData.shares[0], 0, 0);
        awrp->writeShareToOperand(vertexData.shares[1], 0, 1);
        awrp->operands[1].plain_val = (uint32_t)weight;
        awrp->operands[1].plain_parent = (uint64_t)(srcId);

        awrp->isOperandEncrypted[0] = true;
        awrp->isOperandEncrypted[1] = false;

        struct Task task = {ADD_UINT_WITH_REPLACE_PARENT, dstId, false, (void*)awrp, srcId, dstTid, srcTid};
        task.isDummy = isDummy;
        return task;
    }

    struct Task genDummyScatterTask(const Task& src_task) const {
        AddWithReplaceParent<uint32_t, uint64_t>* awrp = new AddWithReplaceParent<uint32_t, uint64_t>();
        AddWithReplaceParent<uint32_t, uint64_t>* src_awrp = (AddWithReplaceParent<uint32_t, uint64_t>*)src_task.buf;
        *awrp = *src_awrp;
        struct Task task = src_task;
        task.buf = (void*)awrp;
        task.isDummy = true;
        return task;
    }

    // struct Task genEmptyScatterTask() const {
    //     AddWithReplaceParent<uint32_t, uint64_t>* awrp = new AddWithReplaceParent<uint32_t, uint64_t>();
    //     struct Task task = {ADD_UINT_WITH_REPLACE_PARENT, -1, false, (void*)awrp, -1, -1, -1};
    //     return task;
    // }

    UpdateType getScatterTaskResult(const Task& task) const {
        if (!task.finished) {
            throw ResultNotReadyException("Task result not ready yet!\n");
        }
        AddWithReplaceParent<uint32_t, uint64_t>* awrp = (AddWithReplaceParent<uint32_t, uint64_t>*)task.buf;

        if (!awrp->isOperandEncrypted[0]) {
            // printf("Scatter result = %d, %lld\n", add->operands[0].plain, task.srcIndex);
            SSSPUpdate<EdgeWeightType> result = SSSPUpdate<EdgeWeightType>(awrp->operands[0].plain_val, awrp->operands[0].plain_parent);
            return result;
        } else {
            SSSPUpdate<EdgeWeightType> result = SSSPUpdate<EdgeWeightType>(awrp->operands[0].enc);
            return result;
        }        
    }

    bool gather(const GraphGASLite::IterCount& iter, Ptr<VertexType>& dst, const UpdateType& update) const {
        auto& data = dst->data();
        if (data.distance > update.distance) {
            data.distance = update.distance;
            data.predecessor = update.predecessor;
            data.activeIter = iter + 1;
            // Not converged.
            return false;
        }
        return true;
    }

    struct Task genGatherTask(const GraphGASLite::IterCount& iter, Ptr<VertexType>& dst, const UpdateType& update) const {
        auto& data = dst->data();
        MinWithParent<uint32_t, uint64_t>* mwp = new MinWithParent<uint32_t, uint64_t>();

        if (!data.isEncrypted) {
            mwp->operands[0].plain_val = data.distance;
            mwp->operands[0].plain_parent = data.predecessor;
        } else {
            mwp->operands[0].enc = data.enc;
        }

        if (!update.isEncrypted) {
            mwp->operands[1].plain_val = update.distance;
            mwp->operands[1].plain_parent = update.predecessor;
        } else {
            mwp->operands[1].enc = update.enc;
        }

        mwp->isOperandEncrypted[0] = data.isEncrypted;
        mwp->isOperandEncrypted[1] = update.isEncrypted;

        struct Task task = {MIN_UINT_WITH_PARENT, (uint64_t)dst->vid(), false, (void*)mwp, (uint64_t)-1, this->getVertexTid(dst->vid()), this->getVertexTid(dst->vid())};
        return task;        
    }

    struct Task genGatherTask(const ShareVec& vertexData, const ShareVec& update, uint64_t dstId, uint64_t dstTid, bool isDummy=false) const {
        MinWithParent<uint32_t, uint64_t>* mwp = new MinWithParent<uint32_t, uint64_t>();

        mwp->setCipherEntryPlainNum(0, 2);
        mwp->setCipherEntryPlainNum(1, 2);
        mwp->writeShareToOperand(vertexData.shares[0], 0, 0);
        mwp->writeShareToOperand(vertexData.shares[1], 0, 1);
        mwp->writeShareToOperand(update.shares[0], 1, 0);
        mwp->writeShareToOperand(update.shares[1], 1, 1);

        mwp->isOperandEncrypted[0] = true;
        mwp->isOperandEncrypted[1] = true;
        mwp->useMask = true;
        mwp->operandMask = true;

        struct Task task = {MIN_UINT_WITH_PARENT, dstId, false, (void*)mwp, (uint64_t)-1, dstTid, dstTid};
        task.isDummy = isDummy;
        return task;        
    }

    void writeGatherTaskResult(const Task& task, ShareVec& dst) const {
        if (!task.finished) {
            throw ResultNotReadyException("Task result not ready yet!\n");
        }
        MinWithParent<uint32_t, uint64_t>* mwp = (MinWithParent<uint32_t, uint64_t>*)task.buf;

        dst.shares.clear();
        dst.shares.push_back(mwp->operands[0].enc.share[0]);
        dst.shares.push_back(mwp->operands[0].enc.share[1]);
    }

    // struct Task genEmptyGatherTask() const {
    //     MinWithParent<uint32_t, uint64_t>* mwp = new MinWithParent<uint32_t, uint64_t>();
    //     struct Task task = {MIN_UINT_WITH_PARENT, -1, false, (void*)mwp, -1, -1, -1};
    //     return task;        
    // }

    void writeGatherTaskResult(const Task& task, Ptr<VertexType>& dst) const {
        if (!task.finished) {
            throw ResultNotReadyException("Task result not ready yet!\n");
        }
        MinWithParent<uint32_t, uint64_t>* mwp = (MinWithParent<uint32_t, uint64_t>*)task.buf;
        auto& data = dst->data();

        if (!mwp->isOperandEncrypted[0]) {
            data.distance = mwp->operands[0].plain_val;
            data.predecessor = GraphGASLite::VertexIdx(mwp->operands[0].plain_parent);
        } else {
            data.enc = mwp->operands[0].enc;
        }

        data.isEncrypted = mwp->isOperandEncrypted[0];
        
        return;
    }

    void onAlgoKernelStart(Ptr<GraphTileType>& graph) const {
        // Set source vertex.
        auto vsrc = graph->vertex(src_);
        // Source vertex does not exist.
        if (vsrc == nullptr) return;
        vsrc->data().distance = 0;
        vsrc->data().activeIter = 0;
    }

    void onAlgoKernelEnd(Ptr<GraphTileType>& graph) const {
        // Nothing to do
    }

    uint32_t getPlainNumPerOperand() const {
        return 2;
    }

protected:
    SSSPEdgeCentricAlgoKernel(const string& name, const GraphGASLite::VertexIdx& src)
        : GraphGASLite::SSEdgeCentricAlgoKernel<GraphTileType>(name),
          src_(src)
    {
        // Nothing else to do.
    }

private:
    const GraphGASLite::VertexIdx src_;
};

#endif // ALGO_KERNELS_EDGE_CENTRIC_SSSP_SSSP_H_

