#ifndef EDGE_CENTRIC_ALGO_KERNEL_H_
#define EDGE_CENTRIC_ALGO_KERNEL_H_

#include "algo_kernel.h"
#include <thread>

namespace GraphGASLite {

template<typename GraphTileType>
class EdgeCentricAlgoKernel : public BaseAlgoKernel<GraphTileType> {
public:
    typedef typename GraphTileType::VertexType VertexType;
    typedef typename GraphTileType::EdgeType::WeightType EdgeWeightType;
    typedef typename GraphTileType::UpdateType UpdateType;
    typedef typename GraphTileType::MirrorVertexType MirrorVertexType;

public:
    AlgoKernelTag tag() const final {
        return AlgoKernelTag::EdgeCentric;
    }

protected:
    /**
     * Edge-centric scatter function.
     *
     * @param iter      Current iteration count.
     * @param src       Source vertex.
     * @param weight    Weight of the edge.
     *
     * @return          A pair consisting of the output update data, and a bool
     *                  denoting whether the update is valid.
     */
    virtual std::pair<UpdateType, bool>
    scatter(const IterCount& iter, Ptr<VertexType>& src, EdgeWeightType& weight) const = 0;

    /**
     * Edge-centric scatter task generation function.
     *
     * @param iter      Current iteration count.
     * @param src       Source vertex.
     * @param weight    Weight of the edge.
     *
     * @return          The Task object.
     */
    virtual struct Task
    genScatterTask(const IterCount& iter, Ptr<VertexType>& src, EdgeWeightType& weight, const GraphGASLite::VertexIdx& dstId) const = 0;

    virtual struct Task genDummyScatterTask(const Task& src_task) const = 0;

    /**
     * Edge-centric scatter task, get result to an UpdateType object.
     *
     * @param task      Task object
     *
     * @return          The UpdateType object.
     */
    virtual UpdateType
    getScatterTaskResult(const Task& task) const = 0;

    /**
     * Edge-centric gather function.
     *
     * @param iter      Current iteration count.
     * @param dst       Destination vertex.
     * @param update    Input update data.
     *
     * @return          Whether this vertex is converged.
     */
    virtual bool
    gather(const IterCount& iter, Ptr<VertexType>& dst, const UpdateType& update) const = 0;

    /**
     * Edge-centric gather task generation function.
     *
     * @param iter      Current iteration count.
     * @param dst       Destination vertex.
     * @param update    Input update data.
     *
     * @return          The Task object.
     */
    virtual struct Task
    genGatherTask(const IterCount& iter, Ptr<VertexType>& dst, const UpdateType& update) const = 0;

    /**
     * Edge-centric gather task, write result to destination point.
     *
     * @param task      Task object
     * @param dst       Destination vertex pointer
     *
     * @return          The UpdateType object.
     */
    virtual void
    writeGatherTaskResult(const Task& task, Ptr<VertexType>& dst) const = 0;

protected:
    using typename BaseAlgoKernel<GraphTileType>::CommSyncType;
    bool onIteration(Ptr<GraphTileType>& graph, CommSyncType& cs, const IterCount& iter) const;
    void runAlgoKernelServer(std::vector<std::thread>& threads) const;
    void closeAlgoKernelServer(std::vector<std::thread>& threads) const;

protected:
    EdgeCentricAlgoKernel(const string& name)
        : BaseAlgoKernel<GraphTileType>(name)
    {
        // Nothing else to do.
    }
};

template<typename GraphTileType>
bool EdgeCentricAlgoKernel<GraphTileType>::
onIteration(Ptr<GraphTileType>& graph, CommSyncType& cs, const IterCount& iter) const {

    printf("tid-> %lld, iteration-> %lld\n", graph->tid(), iter.cnt());
    const auto tid = graph->tid();

    // As there is no explicit barrier in each iteration, the earliest time
    // when the producer can know it is safe to reset comm. utility, is at the
    // beginning of the next iteration after the barrier b/w iterations.
    cs.keyValProdDelAll(tid);

    // Task queue.
    std::queue<Task> taskq;

    // Scatter.
    for (auto edgeIter = graph->edgeIter(); edgeIter != graph->edgeIterEnd(); ++edgeIter) {
        const auto srcId = edgeIter->srcId();
        const auto dstId = edgeIter->dstId();
        // Return reference to allow update to weight.
        auto& weight = edgeIter->weight();

        // Scatter.
        auto src = graph->vertex(srcId);
        auto ret = genScatterTask(iter, src, weight, dstId);
        taskq.push(ret);
    }

    std::cout<<tid<<" "<<"Finish Scatter task generation"<<std::endl;

    // Compute on Scatter task queue.
    this->taskQueueHandler(&taskq);

    std::cout<<tid<<" "<<"Finish Scatter task computation"<<std::endl;

    // if (iter.cnt() == 1) exit(0);

    // Write task results to update data
    while (!taskq.empty()) {
        // printf("1\n");
        auto& ret = taskq.front();
        const auto update = getScatterTaskResult(ret);
        const auto dstId = GraphGASLite::VertexIdx(ret.vertexIndex);
        if (graph->hasVertex(dstId)) {
            // Local destination.
            cs.keyValNew(tid, tid, dstId, update);
        } else {
            // Remote destination, use mirror vertex.
            auto mv = graph->mirrorVertex(dstId);
            // mv->updateNew(update);
            mv->appendUpdateToQueue(update);
        }
        ret.delete_task_content_buf();
        taskq.pop();
    }

    std::cout<<tid<<" "<<"Finish pushing updates to the queue of dst vertex / mirror vertex"<<std::endl;

    // Compute tasks in each mirror vertex's update queue
    std::queue<Ptr<MirrorVertexType>> mvq;
    for (auto mvIter = graph->mirrorVertexIter(); mvIter != graph->mirrorVertexIterEnd(); ++mvIter) {
        auto mv = mvIter->second;
        // if (!mv->hasUpdate()) {
        //     // Skip if no update accumulated.
        //     continue;
        // }
        if (!mv->isUpdateQueueEmpty()) {
            mvq.push(mv);
        }        
    }

    // printf("HERE2\n");

    // The reason for the outer iteration is that a update task in a mirror vertex relies on the previous task's result.
    // The number of the outer interation is linear with the degree of mirror vertices.
    // Future job might optimize it to log(degree)-level.
    int cur_cnt = 0;
    while (!mvq.empty()) {
        std::queue<Task>().swap(taskq);
        std::queue<Ptr<MirrorVertexType>> mvqNextRound;
        while (!mvq.empty()) {
            auto mv = mvq.front();
            mvq.pop();
            mv->getTasks(taskq);
            // taskq.push(mv->getTask());

            if (!mv->isUpdateQueueEmpty()) {
                mvqNextRound.push(mv);
            }
        }

        // printf("--%d\n", taskq.size());
        // Handle taskq using SGX
        // printf("Before taskq.size() = %d\n", taskq.size());
        if (!taskq.empty())
            this->taskQueueHandler(&taskq);
        // printf("After taskq.size() = %d\n", taskq.size());

        // printf("++%d\n", taskq.size());
        
        // Write back results to the mirror vertex.
        while (!taskq.empty()) {
            auto ret = taskq.front();
            auto mv = graph->mirrorVertex(GraphGASLite::VertexIdx(ret.vertexIndex));
            mv->writeTaskResult(ret);
            ret.delete_task_content_buf();
            taskq.pop();
        }

        std::swap(mvq, mvqNextRound);
        cur_cnt++;
    }

    std::cout<<tid<<" "<<"Mirror vertex update queue merge round number: "<<cur_cnt<<std::endl;

    std::cout<<tid<<" "<<"Finish pre-merging updates for mirror vertex"<<std::endl;

    // Send data.
    for (auto mvIter = graph->mirrorVertexIter(); mvIter != graph->mirrorVertexIterEnd(); ++mvIter) {
        auto mv = mvIter->second;
        const auto dstTileId = mv->masterTileId();
        const auto dstId = mv->vid();
        if (!mv->hasUpdate()) {
            // Skip if no update accumulated.
            continue;
        }
        const auto& accUpdate = mv->accUpdate();
        // cs.keyValNew(tid, dstTileId, dstId, accUpdate);
        cs.remoteKeyValNew(tid, dstTileId, dstId, accUpdate);
    }

    // for (uint32_t idx = 0; idx < cs.threadCount(); idx++) {
    //     cs.endTagNew(tid, idx);
    // }

    // End remoteKeyValNew
    for (uint32_t idx = 0; idx < cs.threadCount(); idx++)
        if (idx != tid)
            cs.endRemoteKeyValNew(tid, idx);

    std::cout<<tid<<" "<<"Finish sending accumulated update of mirror vertex"<<std::endl;

    // Remote receive data first
    for (uint32_t idx = 0; idx < cs.threadCount(); idx++) {
        if (idx != tid) {
            while (true) {
                if (cs.getKeyValNew(idx, tid) != 0)
                    break; 
            }
        }
    }

    std::cout<<tid<<" "<<"Finish receiving accumulated update of mirror vertex"<<std::endl;

    // End getKeyValNew
    for (uint32_t idx = 0; idx < cs.threadCount(); idx++)
        if (idx != tid)
            cs.endGetKeyValNew(idx, tid);

    // std::cout<<tid<<" "<<"Finish end recv"<<std::endl;

    // Ensure end getKeyValNew
    for (uint32_t idx = 0; idx < cs.threadCount(); idx++)
        if (idx != tid)
            cs.EnsureEndGetKeyValNew(tid, idx);

    // Clear updates in mirror vertex.
    for (auto mvIter = graph->mirrorVertexIter(); mvIter != graph->mirrorVertexIterEnd(); ++mvIter) {
        auto mv = mvIter->second;
        mv->updateDelAll();
    }

    std::cout<<tid<<" "<<"Begin Gather for each vertex."<<std::endl;

    // Receive data and gather.
    bool converged = false;

    // Partition based on producer tiles.
    auto updatePartitions = cs.keyValTiles(tid);

    for (const auto& prtn : updatePartitions) {
        // For each update ...
        for (const auto& u : prtn) {
            const auto dstId = u.key();
            const auto& update = u.val();

            // Append update to the update queue of the dst vertex
            auto dst = graph->vertex(dstId);
            dst->appendUpdateToQueue(update);
        }
    }

    std::queue<Ptr<VertexType>> vq;
    for (auto vIter = graph->vertexIter(); vIter != graph->vertexIterEnd(); ++vIter) {
        auto v = vIter->second;
        if (!v->hasUpdate()) {
            // Skip if no update accumulated.
            continue;
        }
        if (!v->isUpdateQueueEmpty()) {
            vq.push(v);
        }
    }

    std::cout<<tid<<" "<<"Begin Gather update queue merge."<<std::endl;

    cur_cnt = 0;

    while (!vq.empty()) {
        std::queue<Task>().swap(taskq);
        std::queue<Ptr<VertexType>> vqNextRound;
        while (!vq.empty()) {
            auto v = vq.front();
            vq.pop();
            v->getTasks(taskq);
            // taskq.push(v->getTask());

            if (!v->isUpdateQueueEmpty()) {
                vqNextRound.push(v);
            }
        }

        // printf("---%d\n", taskq.size());
        // Handle taskq using SGX
        // printf("Gather Before taskq.size() = %d\n", taskq.size());
        if (!taskq.empty())
            this->taskQueueHandler(&taskq);
        // printf("Gather After taskq.size() = %d\n", taskq.size());

        // printf("+++%d\n", taskq.size());
        
        // Write back results to the mirror vertex.
        while (!taskq.empty()) {
            auto ret = taskq.front();
            auto v = graph->vertex(GraphGASLite::VertexIdx(ret.vertexIndex));
            v->writeTaskResult(ret);
            ret.delete_task_content_buf();
            taskq.pop();
        }

        std::swap(vq, vqNextRound);

        cur_cnt += 1;
    }

    std::cout<<tid<<" "<<"Gather update queue merge round number: "<<cur_cnt<<std::endl;

    std::cout<<tid<<" "<<"Begin Gather taskq computation."<<std::endl;

    for (auto vIter = graph->vertexIter(); vIter != graph->vertexIterEnd(); ++vIter) {
        auto v = vIter->second;
        if (!v->hasUpdate()) {
            // Skip if no update accumulated.
            continue;
        }
        
        auto ret = genGatherTask(iter, v, v->accUpdate());
        taskq.push(ret);
    }

    // Execute task queue computation.
    this->taskQueueHandler(&taskq);

    // Write back task results to destination vertex.
    while (!taskq.empty()) {
        auto ret = taskq.front();
        const auto dstId = GraphGASLite::VertexIdx(ret.vertexIndex);
        auto v = graph->vertex(dstId);
        writeGatherTaskResult(ret, v);
        ret.delete_task_content_buf();
        taskq.pop();
    }

    std::cout<<tid<<" "<<"Finish Gather and Apply computation"<<std::endl;

    // exit(0);

    cs.keyValConsDelAll(tid);

    for (auto vIter = graph->vertexIter(); vIter != graph->vertexIterEnd(); ++vIter) {
        auto v = vIter->second;
        v->updateDelAll();
    }

    return converged;
}

template<typename GraphTileType>
void EdgeCentricAlgoKernel<GraphTileType>::runAlgoKernelServer(std::vector<std::thread>& threads) const {
    issue_server_recv_threads(threads);
}

template<typename GraphTileType>
void EdgeCentricAlgoKernel<GraphTileType>::closeAlgoKernelServer(std::vector<std::thread>& threads) const {
    TaskComm& clientTaskComm = TaskComm::getClientInstance();
    clientTaskComm.sendFinish();
    for (auto& thrd : threads)
        thrd.join();
}

}

#endif