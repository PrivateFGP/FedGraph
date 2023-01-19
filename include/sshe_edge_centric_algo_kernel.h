#ifndef SSHE_EDGE_CENTRIC_ALGO_KERNEL_H_
#define SSHE_EDGE_CENTRIC_ALGO_KERNEL_H_

#include "edge_centric_algo_kernel.h"
#include <thread>
#include <chrono>

namespace GraphGASLite {

template<typename GraphTileType>
class SSHEEdgeCentricAlgoKernel : public EdgeCentricAlgoKernel<GraphTileType> {
public:
    typedef typename GraphTileType::VertexType VertexType;
    typedef typename GraphTileType::EdgeType::WeightType EdgeWeightType;
    typedef typename GraphTileType::UpdateType UpdateType;
    typedef typename GraphTileType::MirrorVertexType MirrorVertexType;

protected:
    using typename BaseAlgoKernel<GraphTileType>::CommSyncType;
    bool onIteration(Ptr<GraphTileType>& graph, CommSyncType& cs, const IterCount& iter) const final;
    void runAlgoKernelServer(std::vector<std::thread>& threads) const final;
    void closeAlgoKernelServer(std::vector<std::thread>& threads) const final;
    using EdgeCentricAlgoKernel<GraphTileType>::scatter;
    using EdgeCentricAlgoKernel<GraphTileType>::genScatterTask;
    using EdgeCentricAlgoKernel<GraphTileType>::genDummyScatterTask;
    using EdgeCentricAlgoKernel<GraphTileType>::getScatterTaskResult;
    using EdgeCentricAlgoKernel<GraphTileType>::gather;
    using EdgeCentricAlgoKernel<GraphTileType>::genGatherTask;
    using EdgeCentricAlgoKernel<GraphTileType>::writeGatherTaskResult;

protected:
    SSHEEdgeCentricAlgoKernel(const string& name)
        : EdgeCentricAlgoKernel<GraphTileType>(name)
    {
        // Nothing else to do.
    }
};

template<typename GraphTileType>
bool SSHEEdgeCentricAlgoKernel<GraphTileType>::
onIteration(Ptr<GraphTileType>& graph, CommSyncType& cs, const IterCount& iter) const {
    printf("tid-> %lld, iteration-> %lld\n", graph->tid(), iter.cnt());    
    const auto tid = graph->tid();
    TaskComm& clientTaskComm = TaskComm::getClientInstance();
    size_t tileNum = clientTaskComm.getTileNum();

    // As there is no explicit barrier in each iteration, the earliest time
    // when the producer can know it is safe to reset comm. utility, is at the
    // beginning of the next iteration after the barrier b/w iterations.
    cs.keyValProdDelAll(tid);

    // Task queue.
    std::queue<Task> taskq;

    std::cout<<tid<<" "<<"Begin Scatter task generation"<<std::endl;
    // Scatter task generation
    for (auto edgeIter = graph->edgeIter(); edgeIter != graph->edgeIterEnd(); ++edgeIter) {
        const auto srcId = edgeIter->srcId();
        const auto dstId = edgeIter->dstId();
        // Return reference to allow update to weight.
        auto& weight = edgeIter->weight();

        // Scatter.
        auto src = graph->vertex(srcId);
        auto ret = genScatterTask(iter, src, weight, dstId);

        if (graph->hasVertex(dstId)) {
            // Local destination.
            auto dst = graph->vertex(dstId);
            dst->pushToTaskq(ret);
        } else {
            // Remote destination, use mirror vertex.
            auto mv = graph->mirrorVertex(dstId);
            mv->pushToTaskq(ret);
        }
    }

    std::vector<std::vector<Task>> taskvs;
    taskvs.resize(tileNum);
    std::vector<Task> plain_taskv;

    std::cout<<tid<<" "<<"Begin Scatter taskq construction"<<std::endl;
    std::cout<<tid<<" "<<"Begin Scatter taskq splitting"<<std::endl;
    // Split taskq for each vertex & append dummy task & prepare for counting sort
    std::vector<uint64_t> maxv_size;
    maxv_size.resize(tileNum);
    for (int i=0; i<tileNum; ++i) maxv_size[i] = 0;

    for (auto vIter = graph->vertexIter(); vIter != graph->vertexIterEnd(); ++vIter) {
        auto v = vIter->second;
        std::vector<Task>& cur_taskq = v->getTaskq();
        std::vector<std::vector<Task>>& cur_taskvs = v->getTaskvs();
        // std::vector<Task>& cur_plain_taskv = v->getPlainTaskv();
        if (cur_taskq.size() > 0) {
            split_task_queue_on_destination(cur_taskq, cur_taskvs, plain_taskv);
            for (int i=0; i<cur_taskvs.size(); ++i) {
                std::vector<Task>& cur_taskv = cur_taskvs[i];
                if (cur_taskv.size() > 0) {
                    uint64_t dummy_num = get_next_power_of_2(cur_taskv.size()) - cur_taskv.size();
                    for (int i=0; i<dummy_num; ++i) cur_taskv.push_back(genDummyScatterTask(cur_taskv[0]));
                }
                if (cur_taskv.size() > maxv_size[i]) maxv_size[i] = cur_taskv.size();
            }
            cur_taskq.clear();
        }
    }

    for (auto mvIter = graph->mirrorVertexIter(); mvIter != graph->mirrorVertexIterEnd(); ++mvIter) {
        auto mv = mvIter->second;
        std::vector<Task>& cur_taskq = mv->getTaskq();
        std::vector<std::vector<Task>>& cur_taskvs = mv->getTaskvs();
        // std::vector<Task>& cur_plain_taskv = mv->getPlainTaskv();
        if (cur_taskq.size() > 0) {
            split_task_queue_on_destination(cur_taskq, cur_taskvs, plain_taskv);
            for (int i=0; i<cur_taskvs.size(); ++i) {
                std::vector<Task>& cur_taskv = cur_taskvs[i];
                if (cur_taskv.size() > 0) {
                    uint64_t dummy_num = get_next_power_of_2(cur_taskv.size()) - cur_taskv.size();
                    for (int i=0; i<dummy_num; ++i) cur_taskv.push_back(genDummyScatterTask(cur_taskv[0]));
                }
                if (cur_taskv.size() > maxv_size[i]) maxv_size[i] = cur_taskv.size();
            }
            cur_taskq.clear();
        }
    }

    // Counting sort
    std::vector<std::vector<std::vector<uint64_t>>> countingVec;
    countingVec.resize(tileNum);
    for (int i=0; i<tileNum; ++i) {
        countingVec[i].resize(maxv_size[i]);
    }
    for (auto vIter = graph->vertexIter(); vIter != graph->vertexIterEnd(); ++vIter) {
        auto v = vIter->second;
        uint64_t vid = (uint64_t)vIter->second->vid();
        std::vector<std::vector<Task>>& cur_taskvs = v->getTaskvs();
        if (cur_taskvs.size() > 0) {
            for (int i=0; i<cur_taskvs.size(); ++i) {
                std::vector<Task>& cur_taskv = cur_taskvs[i];
                if (cur_taskv.size() > 0) {
                    countingVec[i][cur_taskv.size()-1].push_back(vid);
                }                
            }
        }
    }

    for (auto mvIter = graph->mirrorVertexIter(); mvIter != graph->mirrorVertexIterEnd(); ++mvIter) {
        auto mv = mvIter->second;
        uint64_t vid = (uint64_t)mvIter->second->vid();
        std::vector<std::vector<Task>>& cur_taskvs = mv->getTaskvs();
        if (cur_taskvs.size() > 0) {
            for (int i=0; i<cur_taskvs.size(); ++i) {
                std::vector<Task>& cur_taskv = cur_taskvs[i];
                if (cur_taskv.size() > 0) {
                    countingVec[i][cur_taskv.size()-1].push_back(vid);
                }                
            }
        }
    }

    // Taskvs construction based on counting sort result
    for (int i=0; i<tileNum; ++i) {
        for (int j=countingVec[i].size()-1; j>=0; --j) {
            if (countingVec[i][j].size() > 0) {
                for (int k=0; k<countingVec[i][j].size(); ++k) {
                    // Get vertex or mirror vertex
                    const auto dstId = countingVec[i][j][k];
                    if (graph->hasVertex(dstId)) {
                        auto dst = graph->vertex(dstId);
                        std::vector<std::vector<Task>>& cur_taskvs = dst->getTaskvs();
                        taskvs[i].insert(taskvs[i].end(), cur_taskvs[i].begin(), cur_taskvs[i].end());
                        cur_taskvs[i].clear();
                    } else {
                        auto mv = graph->mirrorVertex(dstId);
                        std::vector<std::vector<Task>>& cur_taskvs = mv->getTaskvs();
                        taskvs[i].insert(taskvs[i].end(), cur_taskvs[i].begin(), cur_taskvs[i].end());
                        cur_taskvs[i].clear();
                    }
                }
            }
        }
    }

    // Clean per-vertex taskvs
    for (auto vIter = graph->vertexIter(); vIter != graph->vertexIterEnd(); ++vIter) {
        auto v = vIter->second;
        uint64_t vid = (uint64_t)vIter->second->vid();
        std::vector<std::vector<Task>>& cur_taskvs = v->getTaskvs();
        if (cur_taskvs.size() > 0) cur_taskvs.clear();
    }

    for (auto mvIter = graph->mirrorVertexIter(); mvIter != graph->mirrorVertexIterEnd(); ++mvIter) {
        auto mv = mvIter->second;
        uint64_t vid = (uint64_t)mvIter->second->vid();
        std::vector<std::vector<Task>>& cur_taskvs = mv->getTaskvs();
        if (cur_taskvs.size() > 0) cur_taskvs.clear();
    }
    
    std::cout<<tid<<" "<<"Finish Scatter task generation"<<std::endl;

    // Scatter & Pre-merge for taskvs
    std::vector<std::thread> cipher_threads;
    for (int i=0; i<taskvs.size(); ++i) {
        if (i != tid) { // Enter even when 'taskvs[i].size() == 0'
            cipher_threads.emplace_back([this, tid, i, &taskvs]() {
                std::vector<Task>& taskv = taskvs[i];
                std::vector<Task> new_taskv;
                TaskComm& clientTaskComm = TaskComm::getClientInstance();
                TaskqHandlerConfig& thc = clientTaskComm.getTaskqHandlerConfig(i);
                uint64_t rotation_ub = get_rotation_upper_bound(taskv);

                // Compute Scatter task vector 
                std::cout<<tid<<" "<<"Begin Scatter task computation for taskvs["<<i<<"]"<<std::endl;
                thc.sendOperand = true;
                thc.mergeResult = false;
                thc.sendTaskqDigest = true;
                thc.rotation = 0;
                single_destination_task_vector_handler(taskv, i);
                std::cout<<tid<<" "<<"Finish Scatter task computation for taskvs["<<i<<"]"<<std::endl;
                   
                std::cout<<tid<<" "<<"Begin pre-merging for taskvs["<<i<<"]"<<std::endl;
                for (auto& ret : taskv) {
                    const auto update = this->getScatterTaskResult(ret);
                    const auto dstId = ret.vertexIndex;
                    const auto dstTid = ret.dstTid;
                    new_taskv.push_back(UpdateType::genTask(update, update, dstId, dstId, dstTid, true, ret.isDummy));
                    ret.delete_task_content_buf();
                }
                taskv.clear();
                std::swap(taskv, new_taskv);

                // Send rotation ub
                clientTaskComm.sendRotationUB(rotation_ub, i);
                printf("Client rotation upper bound %lld\n", rotation_ub);

                thc.sendOperand = false;
                thc.mergeResult = false;
                thc.sendTaskqDigest = false;
                thc.rotation = 1;
                while (thc.rotation < rotation_ub) {
                    auto t1 = std::chrono::high_resolution_clock::now();
                    client_task_vector_direct_handler(taskv, i);
                    thc.rotation = thc.rotation << 1;       
                    auto t2 = std::chrono::high_resolution_clock::now();
                    std::cout << "::client_task_vector_direct_handler "
                        << std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count() / 1000
                    << " seconds\n";
                }

                // thc.rotation == rotation_ub
                // Receive encrypted share in the last round of rotation
                thc.mergeResult = true;
                // single_destination_task_vector_handler(taskv, i);
                client_task_vector_direct_handler(taskv, i);
                client_direct_recv(taskv, i);
                std::cout<<tid<<" "<<"Finish pre-merging for taskvs["<<i<<"]"<<std::endl;
            });
        }
    }

    // Scatter & Pre-merge for plain_taskv
    std::thread plain_thread([this, tid, &plain_taskv]() {
        std::cout<<tid<<" "<<"Begin Scatter&Pre-merging task computation for plain_taskv"<<std::endl;

        std::vector<Task> new_taskv;
        plain_task_vector_handler(plain_taskv);

        bool has_merge = false;
        for (int i=0; i<plain_taskv.size();) {
            const auto dstId = plain_taskv[i].vertexIndex;
            const auto dstTid = plain_taskv[i].dstTid;
            const auto update1 = this->getScatterTaskResult(plain_taskv[i]);
            if (i+1 < plain_taskv.size() && plain_taskv[i].vertexIndex == plain_taskv[i+1].vertexIndex) {
                const auto update2 = this->getScatterTaskResult(plain_taskv[i+1]);
                new_taskv.push_back(UpdateType::genTask(update1, update2, dstId, dstId, dstTid, true));
                plain_taskv[i].delete_task_content_buf();
                plain_taskv[i+1].delete_task_content_buf();
                has_merge = true;
                i += 2;
            } else {
                new_taskv.push_back(UpdateType::genTask(update1, UpdateType(), dstId, dstId, dstTid, true));
                plain_taskv[i].delete_task_content_buf();
                i += 1;
            }
        }
        plain_taskv.clear();
        std::swap(plain_taskv, new_taskv);
        plain_task_vector_handler(plain_taskv);

        while (has_merge == true) {
            has_merge = false;
            for (int i=0; i<plain_taskv.size();) {
                const auto dstId = plain_taskv[i].vertexIndex;
                const auto dstTid = plain_taskv[i].dstTid;
                const auto update1 = UpdateType::getTaskResult(plain_taskv[i]);
                if (i+1 < plain_taskv.size() && plain_taskv[i].vertexIndex == plain_taskv[i+1].vertexIndex) {
                    const auto update2 = UpdateType::getTaskResult(plain_taskv[i+1]);
                    new_taskv.push_back(UpdateType::genTask(update1, update2, dstId, dstId, dstTid, true));
                    plain_taskv[i].delete_task_content_buf();
                    plain_taskv[i+1].delete_task_content_buf();
                    has_merge = true;
                    i += 2;
                } else {
                    new_taskv.push_back(UpdateType::genTask(update1, UpdateType(), dstId, dstId, dstTid, true));
                    plain_taskv[i].delete_task_content_buf();
                    i += 1;
                }
            }
            plain_taskv.clear();
            std::swap(plain_taskv, new_taskv);
            plain_task_vector_handler(plain_taskv);        
        }

        std::cout<<tid<<" "<<"Finish Scatter&Pre-merging task computation for plain_taskv"<<std::endl;
    });

    for (auto& thrd : cipher_threads)
        thrd.join();

    plain_thread.join();
    
    // Now extract the update variable from the task vectors, and append to corresponding vertex / mirror vertex
    std::cout<<tid<<" "<<"Begin pre-merging results extraction"<<std::endl;
    uint64_t cur_vid = -1;
    for (int j=0; j<plain_taskv.size(); ++j) {
        const auto dstId = plain_taskv[j].vertexIndex;
        if (dstId != cur_vid) {
            const auto update = UpdateType::getTaskResult(plain_taskv[j]);
            if (graph->hasVertex(dstId)) {
                auto v = graph->vertex(dstId);
                v->appendUpdateToQueue(update);
            } else {
                auto mv = graph->mirrorVertex(dstId);
                mv->appendUpdateToQueue(update);
            }
            cur_vid = dstId;
        }
        plain_taskv[j].delete_task_content_buf();
    }

    for (int i=0; i<taskvs.size(); ++i) {
        if (i != tid && taskvs[i].size() != 0) {
            cur_vid = -1;
            for (int j=0; j<taskvs[i].size(); ++j) {
                if (taskvs[i][j].isDummy) {
                    taskvs[i][j].delete_task_content_buf();
                    continue;
                }
                const auto dstId = taskvs[i][j].vertexIndex;
                if (dstId != cur_vid) {
                    const auto update = UpdateType::getTaskResult(taskvs[i][j]);
                    if (graph->hasVertex(dstId)) {
                        auto v = graph->vertex(dstId);
                        v->appendUpdateToQueue(update);
                    } else {
                        auto mv = graph->mirrorVertex(dstId);
                        mv->appendUpdateToQueue(update);
                    }
                    cur_vid = dstId;
                }
                taskvs[i][j].delete_task_content_buf();
            }
        }
    }

    std::cout<<tid<<" "<<"Finish pushing updates to the queue of dst vertex / mirror vertex"<<std::endl;

    // Send data.
    for (auto mvIter = graph->mirrorVertexIter(); mvIter != graph->mirrorVertexIterEnd(); ++mvIter) {
        auto mv = mvIter->second;
        if (!mv->hasUpdate()) {
            // Skip if no update accumulated.
            continue;
        }

        const auto dstTileId = mv->masterTileId();
        const auto dstId = mv->vid();
        std::queue<UpdateType>& updateQueue = mv->updateQueue();
        while (!updateQueue.empty()) {
            UpdateType cur_update = updateQueue.front();
            cs.remoteKeyValNew(tid, dstTileId, dstId, cur_update);
            updateQueue.pop();
        }
    }

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

    clientTaskComm.clean();
    std::vector<TaskqHandlerConfig>& thcs = clientTaskComm.getTaskqHandlerConfig();
    for (auto& thc : thcs) {
        thc.rotation = 0;
        thc.mergeResult = true;
        thc.sendOperand = true;
        thc.sendOperand = true;
    }

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

        if (!taskq.empty())
            this->taskQueueHandler(&taskq);
        
        // Write back results to the mirror vertex.
        while (!taskq.empty()) {
            auto ret = taskq.front();
            auto v = graph->vertex(GraphGASLite::VertexIdx(ret.vertexIndex));
            v->writeTaskResult(ret);
            ret.delete_task_content_buf();
            taskq.pop();
        }

        std::swap(vq, vqNextRound);

    }

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
    clientTaskComm.sendFinish();

    cs.keyValConsDelAll(tid);

    for (auto vIter = graph->vertexIter(); vIter != graph->vertexIterEnd(); ++vIter) {
        auto v = vIter->second;
        v->updateDelAll();
    }

    return converged;
}

template<typename GraphTileType>
void SSHEEdgeCentricAlgoKernel<GraphTileType>::runAlgoKernelServer(std::vector<std::thread>& threads) const {
	TaskComm& taskComm = TaskComm::getServerInstance();
	size_t tileNum = taskComm.getTileNum();
	size_t tileIndex = taskComm.getTileIndex();
    uint32_t mpcBasePort = taskComm.getMPCBasePort();

    uint64_t maxIters = this->maxIters().cnt();

	for (int i = 0; i < tileNum; i++) {
		if (i != tileIndex) {
			threads.emplace_back([this, i, &taskComm, tileIndex, tileNum, mpcBasePort, maxIters]() {
                uint64_t iter = 0;
                TaskqHandlerConfig& thc = taskComm.getTaskqHandlerConfig(i);
                std::vector<Task>& taskv = taskComm.getTaskv(i);
                std::vector<Task> new_taskv;
                while (iter < maxIters) { // On iteration
                    printf("Server iter %lld\n", iter);
                    int ret = 0;
                    printf("Server iter %lld Scatter\n", iter);
                    // Scatter
                    thc.sendOperand = true;
                    thc.mergeResult = false;
                    thc.sendTaskqDigest = true;
                    thc.rotation = 0;                
                    while (ret == 0) {
                        ret = taskComm.recvEntrance(i);
                        if (taskComm.getRemoteFinishqFlag(i)) {
                            // taskComm.clean(i);
                            break;
                        }
                    }

                    printf("Server iter %lld pre-merging\n", iter);
                    // Pre-merge
                    for (auto& ret : taskv) {
                        ret.finished = true;
                        const auto update = this->getScatterTaskResult(ret);
                        const auto dstId = ret.vertexIndex;
                        const auto dstTid = ret.dstTid;
                        new_taskv.push_back(UpdateType::genTask(update, update, dstId, dstId, dstTid, true));
                        ret.delete_task_content_buf();
                    }
                    taskv.clear();
                    std::swap(taskv, new_taskv);

                    // Receive rotation ub
                    taskComm.recvEntrance(i);
                    uint64_t rotation_ub = thc.rotation;
                    printf("Server iter %lld rotation upper bound %lld\n", iter, rotation_ub);

                    thc.sendOperand = false;
                    thc.mergeResult = false;
                    thc.sendTaskqDigest = false;
                    thc.rotation = 1;
                    while (thc.rotation < rotation_ub) {
                        server_task_vector_direct_handler(taskv, i);
                        thc.rotation = thc.rotation << 1;
                    }

                    // Send encrypted share in the last round of rotation
                    ret = 0;
                    thc.mergeResult = true;
                    server_task_vector_direct_handler(taskv, i);
                    taskComm.sendEncryptedTaskqShare(taskv, i);
                    taskComm.clean(i);

                    thc.sendOperand = true;
                    thc.mergeResult = true;
                    thc.sendTaskqDigest = true;
                    thc.rotation = 0;

                    printf("Server iter %lld Begin Gather\n", iter);
                    // Gather
                    while (ret == 0) {
                        ret = taskComm.recvEntrance(i);
                        if (taskComm.getRemoteFinishqFlag(i)) {
                            taskComm.clean(i);
                        }
                    }
                    // require(ret == -1)
                    printf("Server iter %lld End Gather\n", iter); 

                    iter++;
                }
                printf("Server finish iteration\n");
                printf("Server begin decryption\n");
                int ret = 0;
                while (ret == 0) {
                    ret = taskComm.recvEntrance(i);
                    if (taskComm.getRemoteFinishqFlag(i)) {
                        taskComm.clean(i);
                    }
                }
                printf("Server finish decryption\n");
			});
		}
	}
}

template<typename GraphTileType>
void SSHEEdgeCentricAlgoKernel<GraphTileType>::closeAlgoKernelServer(std::vector<std::thread>& threads) const {
    TaskComm& clientTaskComm = TaskComm::getClientInstance();
    clientTaskComm.sendFinish();
    for (auto& thrd : threads)
        thrd.join();
}

}

#endif