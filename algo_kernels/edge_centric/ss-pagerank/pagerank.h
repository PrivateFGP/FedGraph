#ifndef ALGO_KERNELS_EDGE_CENTRIC_PAGERANK_PAGERANK_H_
#define ALGO_KERNELS_EDGE_CENTRIC_PAGERANK_PAGERANK_H_

#include "graph.h"
#include "ss_edge_centric_algo_kernel.h"
#include "task.h"
#include "TaskUtil.h"

#include <queue>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>

/*
 * Graph types definitions.
 */
struct PageRankData {
    double rank;
    CipherEntry enc;
    bool isEncrypted;

    PageRankData(const GraphGASLite::VertexIdx&)
        : rank(0), isEncrypted(false) {
        // Nothing else to do.
    }

    void intoShareVec(ShareVec& sv0, ShareVec& sv1) {
        int64_t rank_share0, rank_share1;
        CryptoUtil::intoShares(rank, rank_share0, rank_share1);
        sv0.shares.push_back(rank_share0);
        sv1.shares.push_back(rank_share1);
    }

    void fromShareVec(const ShareVec& sv0, const ShareVec& sv1) {
        rank = CryptoUtil::mergeShareAsDouble(sv0.shares[0], sv1.shares[0]);
        // std::cout<<"Reconstruct share: "<<sv0.shares[0]<<std::endl;
        // std::cout<<"rank: "<<rank<<std::endl;
    }
};

struct PageRankUpdate {
    double contribute;
    CipherEntry enc;
    bool isEncrypted;

    PageRankUpdate(const double contribute_ = 0)
        : contribute(contribute_), isEncrypted(false)
    {
        // Nothing else to do.
    }

    PageRankUpdate(const CipherEntry& cipher)
        : enc(cipher), isEncrypted(true)
    {
        // Nothing else to do.
    }

    PageRankUpdate& operator+=(const PageRankUpdate& update) {
        contribute += update.contribute;
        return *this;
    }

    friend class boost::serialization::access;

    template <typename Archive> 
    void serialize(Archive &ar, const unsigned int version) 
    { 
        ar & contribute; 
        ar & enc;
        ar & isEncrypted;
    }

    static Task genTask(const PageRankUpdate& update1, const PageRankUpdate& update2, uint64_t vertexIndex, uint64_t srcTid, uint64_t dstTid, bool isFinal=false, bool isDummy=false) {
        Add<double>* add = new Add<double>();

        if (!update1.isEncrypted) {
            add->operands[0].plain = update1.contribute;
        } else {
            add->operands[0].enc = update1.enc;
        }

        if (!update2.isEncrypted) {
            add->operands[1].plain = update2.contribute;
        } else {
            add->operands[1].enc = update2.enc;
        }

        add->isOperandEncrypted[0] = update1.isEncrypted;
        add->isOperandEncrypted[1] = update2.isEncrypted;

        struct Task task = {ADD_DOUBLE, vertexIndex, false, (void*)add, (uint64_t)-1, dstTid, srcTid, 0, isFinal, isDummy};
        return task;
    }   

    static PageRankUpdate getTaskResult(const struct Task& task) {
        if (!task.finished) {
            throw ResultNotReadyException("Task result not ready yet!\n");
        }

        Add<double>* add = (Add<double>*)task.buf;
        if (!add->isOperandEncrypted[0]) {
            PageRankUpdate result = PageRankUpdate(add->operands[0].plain);
            return result;
        } else {
            PageRankUpdate result = PageRankUpdate(add->operands[0].enc);
            return result;
        }
    }

    static void getTaskResult(const struct Task& task, ShareVec& sv) {
        if (!task.finished) {
            throw ResultNotReadyException("Task result not ready yet!\n");
        }
        Add<double>* add = (Add<double>*)task.buf;

        if (!add->isOperandEncrypted[0]) {
            printf("Pre-merging results ought to be shared!\n");
            exit(-1);
        } else {
            if (add->operands[0].enc.plainNum != 1) {
                printf("Error results update.enc.plainNum %d\n", add->operands[0].enc.plainNum);
                exit(0);
            }
            sv.shares.clear();
            sv.shares.push_back(add->operands[0].enc.share[0]);
        }
    }

    static PageRankUpdate splitRandomShareFromUpdate(const uint32_t srcTid, const uint32_t dstTid, PageRankUpdate& update) {
        // The update variable must be encrypted by the public key of dstTid
        // Ought not to be called
        PageRankUpdate update_share;
        return update_share;
    }

    static void mergeUpdateShare(const uint32_t dstTid, PageRankUpdate& update, PageRankUpdate& update_share) {
       // Ought not to be called
    }
};

void print(const PageRankUpdate& update) {
    // printf("PageRankUpdate---------\n");
    // printf("count %d\n", update.count.cnt());
    // printf(update.isEncrypted? "True":"False");
    // printf("---------\n");
}

/*
 * Algorithm kernel definition.
 */
template<typename GraphTileType>
class PageRankEdgeCentricAlgoKernel : public GraphGASLite::SSEdgeCentricAlgoKernel<GraphTileType> {
public:
    static Ptr<PageRankEdgeCentricAlgoKernel> instanceNew(const string& name,
            const double beta, const double tolerance) {
        return Ptr<PageRankEdgeCentricAlgoKernel>(new PageRankEdgeCentricAlgoKernel(name, beta, tolerance));
    }

protected:
    typedef typename GraphTileType::UpdateType UpdateType;
    typedef typename GraphTileType::VertexType VertexType;
    typedef typename GraphTileType::EdgeType::WeightType EdgeWeightType;

    std::pair<UpdateType, bool> scatter(const GraphGASLite::IterCount&, Ptr<VertexType>& src, EdgeWeightType& weight) const {
        std::pair<UpdateType, bool> ret;
        ret.first = PageRankUpdate(0);
        ret.second = true;
        return ret;
    }

    struct Task genScatterTask(const GraphGASLite::IterCount& iter, Ptr<VertexType>& src, EdgeWeightType& weight, const GraphGASLite::VertexIdx& dstId) const {
        struct Task task;
        return task;
    }

    struct Task genScatterTask(Ptr<VertexType> src, const ShareVec& vertexData, EdgeWeightType weight, uint64_t srcId, uint64_t srcTid, uint64_t dstId, uint64_t dstTid, bool isDummy=false) const {
        double odeg = 0;
        if (src != nullptr) {
            odeg = (double)src->outDeg().cnt();
            // printf("Odeg %lf\n", odeg);
        }
        
        Div<double>* div = new Div<double>();
        if (vertexData.shares.size() != 1) {
            printf("Wrong vertex data shares length for PageRank Scatter");
        }

        div->setCipherEntryPlainNum(0, 1);
        div->writeShareToOperand(vertexData.shares[0], 0, 0);

        div->operands[1].plain = odeg;
        if (div->operands[1].plain <= 0) div->operands[1].plain = 1;

        div->isOperandEncrypted[0] = true;
        div->isOperandEncrypted[1] = false;

        struct Task task = {DIV_DOUBLE, dstId, false, (void*)div, srcId, dstTid, srcTid};
        task.isDummy = isDummy;
        return task;
    }

    struct Task genDummyScatterTask(const Task& src_task) const {
        struct Task task;
        return task;
    }

    UpdateType getScatterTaskResult(const Task& task) const {
        if (!task.finished) {
            throw ResultNotReadyException("Task result not ready yet!\n");
        }
        Div<double>* div = (Div<double>*)task.buf;

        if (!div->isOperandEncrypted[0]) {
            PageRankUpdate result = PageRankUpdate(div->operands[0].plain);
            return result; 
        } else {
            PageRankUpdate result = PageRankUpdate(div->operands[0].enc);
            return result; 
        }       
    }

    bool gather(const GraphGASLite::IterCount&, Ptr<VertexType>& dst, const UpdateType& update) const {
        return true;
    }

    struct Task genGatherTask(const GraphGASLite::IterCount& iter, Ptr<VertexType>& dst, const UpdateType& update) const {
        struct Task task;
        return task;        
    }

    struct Task genGatherTask(const ShareVec& vertexData, const ShareVec& update, uint64_t dstId, uint64_t dstTid, bool isDummy=false) const {
        Add<double>* add = new Add<double>();

        add->setCipherEntryPlainNum(0, 1);
        add->setCipherEntryPlainNum(1, 1);
        add->writeShareToOperand(vertexData.shares[0], 0, 0);
        add->writeShareToOperand(update.shares[0], 1, 0);

        add->isOperandEncrypted[0] = true;
        add->isOperandEncrypted[1] = true;
        add->useMask = true;
        add->operandMask = true;

        struct Task task = {ADD_DOUBLE, dstId, false, (void*)add, (uint64_t)-1, dstTid, dstTid};
        task.isDummy = isDummy;
        return task;        
    }

    void writeGatherTaskResult(const Task& task, Ptr<VertexType>& dst) const {
    }

    void writeGatherTaskResult(const Task& task, ShareVec& dst) const {
        if (!task.finished) {
            throw ResultNotReadyException("Task result not ready yet!\n");
        }
        Add<double>* add = (Add<double>*)task.buf;

        dst.shares.clear();
        dst.shares.push_back(add->operands[0].enc.share[0]);
    }

    void onAlgoKernelStart(Ptr<GraphTileType>& graph) const {
        for (auto vertexIter = graph->vertexIter(); vertexIter != graph->vertexIterEnd(); ++vertexIter) {
            auto& v = vertexIter->second;
            v->data().rank = 1;
        }
    }

    void onAlgoKernelEnd(Ptr<GraphTileType>& graph) const {
    }

    uint32_t getPlainNumPerOperand() const {
        return 1;
    }

protected:
    PageRankEdgeCentricAlgoKernel(const string& name, const double beta, const double tolerance)
        : GraphGASLite::SSEdgeCentricAlgoKernel<GraphTileType>(name),
          beta_(beta), tolerance_(tolerance)
    {
        // Nothing else to do.
    }

private:
    const double beta_;
    const double tolerance_;
};

#endif // ALGO_KERNELS_EDGE_CENTRIC_PAGERANK_PAGERANK_H_

