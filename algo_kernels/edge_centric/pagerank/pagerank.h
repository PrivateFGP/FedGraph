#ifndef ALGO_KERNELS_EDGE_CENTRIC_PAGERANK_PAGERANK_H_
#define ALGO_KERNELS_EDGE_CENTRIC_PAGERANK_PAGERANK_H_

#include "graph.h"
#include "edge_centric_algo_kernel.h"
#include "task.h"

#include <queue>

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

/*
 * Graph types definitions.
 */
struct PageRankData {
    GraphGASLite::DegreeCount collected;
    union {
        double sum;
        CipherEntry enc_sum;
    };
    union {
        double rank;
        CipherEntry enc_rank;
    };
    bool isSumEncrypted;
    bool isRankEncrypted;

    PageRankData(const GraphGASLite::VertexIdx&)
        : collected(0), sum(0), rank(0), isSumEncrypted(false), isRankEncrypted(false)
    {
        // Nothing else to do.
    }
};

struct PageRankUpdate {
    GraphGASLite::DegreeCount count;
    double contribute;
    CipherEntry enc_contribute; // Cipher Text of concatenated contribute and count
    bool isEncrypted;

    PageRankUpdate(const double contribute_ = 0, const GraphGASLite::DegreeCount& count_ = 0)
        : contribute(contribute_), count(count_), isEncrypted(false)
    {
        // Nothing else to do.
    }

    PageRankUpdate(const CipherEntry& cipher, const GraphGASLite::DegreeCount& count_ = 0)
        : enc_contribute(cipher), count(count_), isEncrypted(true)
    {
        // Nothing else to do.
    }

    PageRankUpdate& operator+=(const PageRankUpdate& update) {
        contribute += update.contribute;
        count += update.count;
        return *this;
    }

    friend class boost::serialization::access;

    template <typename Archive> 
    void serialize(Archive &ar, const unsigned int version) 
    { 
        ar & count;
        ar & contribute; 
        ar & enc_contribute;
        ar & isEncrypted;
    }

    static Task genTask(const PageRankUpdate& update1, const PageRankUpdate& update2, uint64_t vertexIndex, uint64_t srcTid, uint64_t dstTid, bool isFinal=false, bool isDummy=false) {
        AddMixedPair<double, uint32_t>* ap = new AddMixedPair<double, uint32_t>();

        if (!update1.isEncrypted) {
            ap->operands[0].plain_a = update1.contribute;
        } else {
            ap->operands[0].enc_a = update1.enc_contribute;
        }

        if (!update2.isEncrypted) {
            ap->operands[1].plain_a = update2.contribute;
        } else {
            ap->operands[1].enc_a = update2.enc_contribute;
        }

        ap->operands[0].plain_b = update1.count.cnt();
        ap->operands[1].plain_b = update2.count.cnt();

        ap->isOperandEncrypted[0] = update1.isEncrypted;
        ap->isOperandEncrypted[1] = update2.isEncrypted;

        struct Task task = {ADD_MIXED_PAIR_DOUBLE_UINT, vertexIndex, false, (void*)ap, (uint64_t)-1, dstTid, srcTid, 0, isFinal, isDummy};
        return task;
    }   

    static PageRankUpdate getTaskResult(const struct Task& task) {
        if (!task.finished) {
            throw ResultNotReadyException("Task result not ready yet!\n");
        }
        AddMixedPair<double, uint32_t>* ap = (AddMixedPair<double, uint32_t>*)task.buf;
        if (!ap->isOperandEncrypted[0]) {
            PageRankUpdate result = PageRankUpdate(ap->operands[0].plain_a, GraphGASLite::DegreeCount(ap->operands[0].plain_b));
            return result;
        } else {
            PageRankUpdate result = PageRankUpdate(ap->operands[0].enc_a, GraphGASLite::DegreeCount(ap->operands[0].plain_b));
            return result;
        }
    }

    static PageRankUpdate splitRandomShareFromUpdate(const uint32_t srcTid, const uint32_t dstTid, PageRankUpdate& update) {
        // The update variable must be encrypted by the public key of dstTid
        PageRankUpdate update_share;
        CryptoUtil& cryptoUtil = CryptoUtil::getInstance();
        if (!update.isEncrypted) {
            update.enc_contribute = cryptoUtil.encryptToCipherEntry((double)update.contribute, srcTid);
            update.contribute = 0;
            update.isEncrypted = true;
            update_share.isEncrypted = false;
            update_share.contribute = 0;
        } else {
            int64_t new_share = cryptoUtil.splitRandomShareFromCipherEntry(update.enc_contribute);
            update_share.enc_contribute = update.enc_contribute; // Legal?
            update_share.isEncrypted = true;
            update.enc_contribute = cryptoUtil.encryptToCipherEntry(new_share, srcTid);
            update.isEncrypted = true;
        }

        return update_share;
    }

    static void mergeUpdateShare(const uint32_t dstTid, PageRankUpdate& update, PageRankUpdate& update_share) {
        CryptoUtil& cryptoUtil = CryptoUtil::getInstance();
        if (update_share.isEncrypted) {
            int64_t share = cryptoUtil.decryptCipherEntryToInt64(update_share.enc_contribute);
            cryptoUtil.mergeShareIntoCipherEntry(update.enc_contribute, share, update.enc_contribute.tid);
        }

        return;
    }
};

void print(const PageRankUpdate& update) {
    printf("PageRankUpdate---------\n");
    printf("count %d\n", update.count.cnt());
    printf(update.isEncrypted? "True":"False");
    printf("---------\n");
}

/*
 * Algorithm kernel definition.
 */
template<typename GraphTileType>
class PageRankEdgeCentricAlgoKernel : public GraphGASLite::EdgeCentricAlgoKernel<GraphTileType> {
public:
    static Ptr<PageRankEdgeCentricAlgoKernel> instanceNew(const string& name,
            const double beta, const double tolerance) {
        return Ptr<PageRankEdgeCentricAlgoKernel>(new PageRankEdgeCentricAlgoKernel(name, beta, tolerance));
    }

protected:
    typedef typename GraphTileType::UpdateType UpdateType;
    typedef typename GraphTileType::VertexType VertexType;
    typedef typename GraphTileType::EdgeType::WeightType EdgeWeightType;

    std::pair<UpdateType, bool> scatter(const GraphGASLite::IterCount&, Ptr<VertexType>& src, EdgeWeightType&) const {
        auto& data = src->data();
        auto odeg = src->outDeg();
        auto contribute = data.rank / odeg;
        std::pair<UpdateType, bool> ret;
        ret.first = PageRankUpdate(contribute, 1);
        ret.second = true;
        return ret;
    }

    struct Task genScatterTask(const GraphGASLite::IterCount& iter, Ptr<VertexType>& src, EdgeWeightType& weight, const GraphGASLite::VertexIdx& dstId) const {
        auto& data = src->data();
        auto odeg = src->outDeg();
        Div<double>* div = new Div<double>();

        if (!data.isRankEncrypted) {
            div->operands[0].plain = (double)data.rank;
        } else {
            div->operands[0].enc = data.enc_rank;
        }
        div->operands[1].plain = (double)odeg.cnt();

        div->isOperandEncrypted[0] = data.isRankEncrypted;
        div->isOperandEncrypted[1] = false;

        struct Task task = {DIV_DOUBLE, (uint64_t)dstId, false, (void*)div, (uint64_t)(src->vid()), this->getVertexTid(dstId), this->getVertexTid(src->vid())};
        return task;
    }

    struct Task genDummyScatterTask(const Task& src_task) const {
        Div<double>* div = new Div<double>();
        Div<double>* src_div = (Div<double>*)src_task.buf;
        *div = *src_div;
        struct Task task = src_task;
        task.buf = (void*)div;
        task.isDummy = true;
        return task;
    }

    UpdateType getScatterTaskResult(const Task& task) const {
        if (!task.finished) {
            throw ResultNotReadyException("Task result not ready yet!\n");
        }
        Div<double>* div = (Div<double>*)task.buf;

        if (!div->isOperandEncrypted[0]) {
            PageRankUpdate result = PageRankUpdate(div->operands[0].plain, GraphGASLite::DegreeCount(1));
            return result; 
        } else {
            PageRankUpdate result = PageRankUpdate(div->operands[0].enc, GraphGASLite::DegreeCount(1));
            return result; 
        }       
    }

    bool gather(const GraphGASLite::IterCount&, Ptr<VertexType>& dst, const UpdateType& update) const {
        auto& data = dst->data();
        data.sum += update.contribute;
        data.collected += update.count;
        if (data.collected == dst->inDeg()) {
            double newRank = beta_ * data.sum + (1 - beta_);
            bool converge = (std::abs(newRank - data.rank) <= tolerance_);
            data.rank = newRank;
            data.sum = 0;
            data.collected = 0;
            return converge;
        }
        // Convergency is unknown until all updates are collected.
        return true;
    }

    struct Task genGatherTask(const GraphGASLite::IterCount& iter, Ptr<VertexType>& dst, const UpdateType& update) const {
        auto& data = dst->data();
        AddMixedPair<double, uint32_t>* ap = new AddMixedPair<double, uint32_t>();

        if (!data.isSumEncrypted) {
            ap->operands[0].plain_a = data.sum;
        } else {
            ap->operands[0].enc_a = data.enc_sum;
        }
        
        if (!update.isEncrypted) {
            ap->operands[1].plain_a = update.contribute;
        } else {
            ap->operands[1].enc_a = update.enc_contribute;
        }
        
        ap->operands[0].plain_b = data.collected.cnt();
        ap->operands[1].plain_b = update.count.cnt();

        ap->isOperandEncrypted[0] = data.isSumEncrypted;
        ap->isOperandEncrypted[1] = update.isEncrypted;

        struct Task task = {ADD_PAIR_DOUBLE_UINT, dst->vid(), false, (void*)ap, (uint64_t)-1, this->getVertexTid(dst->vid()), this->getVertexTid(dst->vid())};
        // TODO: Add convergency handle
        //       What if there are multiple tasks in Gather? Dependency?
        return task;        
    }

    void writeGatherTaskResult(const Task& task, Ptr<VertexType>& dst) const {
        if (!task.finished) {
            throw ResultNotReadyException("Task result not ready yet!\n");
        }
        AddMixedPair<double, uint32_t>* ap = (AddMixedPair<double, uint32_t>*)task.buf;
        auto& data = dst->data();

        if (!ap->isOperandEncrypted[0]) {
            data.sum = ap->operands[0].plain_a;
        } else {
            data.enc_sum = ap->operands[0].enc_a;
        }
        data.collected = GraphGASLite::DegreeCount(ap->operands[0].plain_b);

        data.isSumEncrypted = ap->isOperandEncrypted[0];

        if (data.collected == dst->inDeg()) {
            if (!data.isSumEncrypted) {
                data.rank = data.sum;
            } else {
                data.enc_rank = data.enc_sum;
            }
            data.isRankEncrypted = data.isSumEncrypted;

            data.sum = 0;
            data.isSumEncrypted = false;
            data.collected = 0;
        }
        
        return;
    }

    Task genDecryptRankTask(const PageRankData& data, uint64_t vertexIndex) const {
        DecryptEntry<double>* dec = new DecryptEntry<double>();

        if (!data.isRankEncrypted) {
            // Exception
            throw InvalidArgumentException("Expect encrypted value!\n");
        }

        dec->operand.enc = data.enc_rank;
        dec->isOperandEncrypted[0] = true;

        struct Task task = {DEC_DOUBLE, vertexIndex, false, (void*)dec, (uint64_t)-1, this->getVertexTid(vertexIndex), this->getVertexTid(vertexIndex)};
        return task;
    }

    void writeDecryptRankTaskResult(const Task& task, Ptr<VertexType>& dst) const {
        if (!task.finished) {
            throw ResultNotReadyException("Task result not ready yet!\n");
        }
        DecryptEntry<double>* dec = (DecryptEntry<double>*)task.buf;
        auto& data = dst->data();

        if (dec->isOperandEncrypted[0]) {
            // Exception
            throw InvalidArgumentException("Expect plain value!\n");
        }

        data.isRankEncrypted = false;
        data.rank = dec->operand.plain;
        
        return;
    }

    void onAlgoKernelStart(Ptr<GraphTileType>& graph) const {
        for (auto vertexIter = graph->vertexIter(); vertexIter != graph->vertexIterEnd(); ++vertexIter) {
            auto& v = vertexIter->second;
            v->data().rank = 1;
        }
    }

    void onAlgoKernelEnd(Ptr<GraphTileType>& graph) const {
        // Decrypt vertex data (rank)

        std::queue<Task> taskq;
        // Task generation
        for (auto vertexIter = graph->vertexIter(); vertexIter != graph->vertexIterEnd(); ++vertexIter) {
            auto& v = vertexIter->second;
            if (v->data().isRankEncrypted) {
                taskq.push(genDecryptRankTask(v->data(), v->vid()));
            }
        }

        task_queue_handler(&taskq);

        while (!taskq.empty()) {
            auto ret = taskq.front();
            const auto dstId = GraphGASLite::VertexIdx(ret.vertexIndex);
            auto v = graph->vertex(dstId);
            writeDecryptRankTaskResult(ret, v);
            ret.delete_task_content_buf();
            taskq.pop();
        }

        // For those who have no in-edges, gather() will never be invoked. Set the rank to the teleport value.
        for (auto vertexIter = graph->vertexIter(); vertexIter != graph->vertexIterEnd(); ++vertexIter) {
            auto& v = vertexIter->second;
            if (v->inDeg() == 0) {
                assert(std::abs(v->data().rank - 1) <= tolerance_);  // initial value
                v->data().rank = 1 - beta_;
            }
        }
    }

protected:
    PageRankEdgeCentricAlgoKernel(const string& name, const double beta, const double tolerance)
        : GraphGASLite::EdgeCentricAlgoKernel<GraphTileType>(name),
          beta_(beta), tolerance_(tolerance)
    {
        // Nothing else to do.
    }

private:
    const double beta_;
    const double tolerance_;
};

#endif // ALGO_KERNELS_EDGE_CENTRIC_PAGERANK_PAGERANK_H_

