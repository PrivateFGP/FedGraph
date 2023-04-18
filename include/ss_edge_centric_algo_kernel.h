#ifndef SS_EDGE_CENTRIC_ALGO_KERNEL_H_
#define SS_EDGE_CENTRIC_ALGO_KERNEL_H_

#include "edge_centric_algo_kernel.h"
#include "ObliviousMapper.h"

#include <thread>
#include <chrono>

namespace GraphGASLite {

template<typename GraphTileType>
class SSEdgeCentricAlgoKernel : public EdgeCentricAlgoKernel<GraphTileType> {
public:
    typedef typename GraphTileType::VertexType VertexType;
    typedef typename GraphTileType::EdgeType::WeightType EdgeWeightType;
    typedef typename GraphTileType::UpdateType UpdateType;
    typedef typename GraphTileType::MirrorVertexType MirrorVertexType;

    struct GraphSummary {
        std::vector<Ptr<VertexType>> localVertexVec;
        std::vector<std::vector<Ptr<MirrorVertexType>>> mirrorVertexVecs;
        std::vector<std::vector<EdgeWeightType>> localEdgeWeightVecs;
        ShareVecVec localVertexSvv;
        std::vector<ShareVecVec> remoteVertexSvvs;
        std::vector<std::vector<bool>> isUpdateSrcVertexDummy;
        std::vector<std::vector<bool>> isGatherDstVertexDummy;
        std::vector<std::vector<uint64_t>> updateSrcVertexPos;
        std::vector<std::vector<uint64_t>> updateDstVertexPos;
        std::vector<uint64_t> localVertexPos;
        std::vector<std::vector<uint64_t>> mirrorVertexPos;
        std::vector<std::vector<uint64_t>> remoteMirrorVertexPos;
        std::vector<uint64_t> rotationUb;
        std::vector<ShareVecVec> remoteUpdateSvvs;
        std::vector<ShareVecVec> localUpdateSvvs;
    };
    typedef struct GraphSummary GraphSummary;

protected:
    using typename BaseAlgoKernel<GraphTileType>::CommSyncType;
    
    void getTwoPartyVertexDataVectorShare(GraphSummary& gs, ShareVecVec& vertexSvv0, ShareVecVec& vertexSvv1) const;
    void mergeTwoPartyVertexDataVectorShare(GraphSummary& gs, ShareVecVec& vertexSvv0, ShareVecVec& vertexSvv1) const;
    bool onIteration(Ptr<GraphTileType>& graph, CommSyncType& cs, GraphSummary& gs, const IterCount& iter) const;
    bool onIteration(Ptr<GraphTileType>& graph, CommSyncType& cs, const IterCount& iter) const {}
    void onPreprocessClient(Ptr<GraphTileType>& graph, CommSyncType& cs, GraphSummary& gs, bool doOMPreprocess = true) const;
    void onPreprocessServer(std::vector<std::thread>& threads, bool doOMPreprocess = true) const;
    void runAlgoKernelServer(std::vector<std::thread>& threads, CommSyncType& cs, GraphSummary& gs) const;
    void runAlgoKernelServer(std::vector<std::thread>& threads) const {}
    void closeAlgoKernelServer(std::vector<std::thread>& threads) const;
    void operator()(Ptr<GraphTileType>& graph, CommSyncType& cs) const override; // Override
    using EdgeCentricAlgoKernel<GraphTileType>::scatter;
    using EdgeCentricAlgoKernel<GraphTileType>::genScatterTask;
    using EdgeCentricAlgoKernel<GraphTileType>::genDummyScatterTask;
    using EdgeCentricAlgoKernel<GraphTileType>::getScatterTaskResult;
    using EdgeCentricAlgoKernel<GraphTileType>::gather;
    using EdgeCentricAlgoKernel<GraphTileType>::genGatherTask;
    using EdgeCentricAlgoKernel<GraphTileType>::writeGatherTaskResult;
    virtual uint32_t getPlainNumPerOperand() const = 0;
    virtual struct Task genScatterTask(Ptr<VertexType> src, const ShareVec& vertexData, EdgeWeightType weight, uint64_t srcId, uint64_t srcTid, uint64_t dstId, uint64_t dstTid, bool isDummy=false) const = 0;
    virtual void writeGatherTaskResult(const Task& task, ShareVec& dst) const = 0;
    virtual struct Task genGatherTask(const ShareVec& vertexData, const ShareVec& update, uint64_t dstId, uint64_t dstTid, bool isDummy=false) const = 0;
    void fromScatterTaskvResultToPreMergingTaskv(std::vector<Task>& taskv) const;

protected:
    SSEdgeCentricAlgoKernel(const string& name)
        : EdgeCentricAlgoKernel<GraphTileType>(name)
    {
        // Nothing else to do.
    }
};

template<typename GraphTileType>
void SSEdgeCentricAlgoKernel<GraphTileType>::
getTwoPartyVertexDataVectorShare(GraphSummary& gs, ShareVecVec& vertexSvv0, ShareVecVec& vertexSvv1) const {
    for (uint64_t i=0; i<gs.localVertexVec.size(); ++i) {
        const auto& v = gs.localVertexVec[i];
        auto& data = v->data();
        ShareVec cur_sv0;
        ShareVec cur_sv1;
        data.intoShareVec(cur_sv0, cur_sv1);
        vertexSvv0.shareVecs.push_back(cur_sv0);
        vertexSvv1.shareVecs.push_back(cur_sv1);
        if (cur_sv0.shares[0] + cur_sv1.shares[0] == 0) printf(">>>>>> %d\n", i);
    }    
}

template<typename GraphTileType>
void SSEdgeCentricAlgoKernel<GraphTileType>::
mergeTwoPartyVertexDataVectorShare(GraphSummary& gs, ShareVecVec& vertexSvv0, ShareVecVec& vertexSvv1) const {
    for (int i=0; i<gs.localVertexVec.size(); ++i) {
        auto& data = gs.localVertexVec[i]->data();
        data.fromShareVec(vertexSvv0.shareVecs[i], vertexSvv1.shareVecs[i]);
    }    
}

template<typename GraphTileType>
void SSEdgeCentricAlgoKernel<GraphTileType>::
operator()(Ptr<GraphTileType>& graph, CommSyncType& cs) const {
    TaskComm& clientTaskComm = TaskComm::getClientInstance();
    size_t tileNum = clientTaskComm.getTileNum();
    size_t tileIndex = clientTaskComm.getTileIndex();
    std::cout<<tileIndex<<" "<<"Initialize graph algo kernel"<<std::endl;
    
    this->onAlgoKernelStart(graph);

    std::cout<<tileIndex<<" "<<"Begin graph preprocessing"<<std::endl;
    
    GraphSummary gs;

    // Preprocessing
    std::vector<std::thread> preprocessServerThreads;

    bool doPreprocess = !clientTaskComm.getNoPreprocess();

    this->onPreprocessServer(preprocessServerThreads, doPreprocess);

    auto t_preprocess = std::chrono::high_resolution_clock::now();

    this->onPreprocessClient(graph, cs, gs, doPreprocess);
    
    print_duration(t_preprocess, "preprocess");

    // this->onPreprocessServer(preprocessServerThreads, false);
    // this->onPreprocessClient(graph, cs, gs, false);
    
    for (auto& thrd : preprocessServerThreads)
        thrd.join();

    std::cout<<tileIndex<<" "<<"Begin vertex data sharing"<<std::endl;
    // Share Vertex data
    ShareVecVec& localVertexSvv = gs.localVertexSvv;
    ShareVecVec remoteLocalVertexSvv;
    this->getTwoPartyVertexDataVectorShare(gs, localVertexSvv, remoteLocalVertexSvv);

    for (int i=0; i<tileNum; ++i) {
        if (i != tileIndex)
            clientTaskComm.sendShareVecVec(remoteLocalVertexSvv, i);
    }

    TaskComm& serverTaskComm = TaskComm::getServerInstance();
    std::vector<ShareVecVec>& remoteVertexSvvs = gs.remoteVertexSvvs;    
    remoteVertexSvvs.resize(tileNum);
    for (int i=0; i<tileNum; ++i) {
        if (i != tileIndex)
            serverTaskComm.recvShareVecVec(remoteVertexSvvs[i], i);
        std::cout<<tileIndex<<" Preprocess "<<remoteVertexSvvs[i].shareVecs.size()<<" "<<i<<std::endl;
    }

    std::cout<<tileIndex<<" "<<"Begin algo kernel iteration"<<std::endl;

    std::vector<std::thread> algo_kernel_server_threads;
    this->runAlgoKernelServer(algo_kernel_server_threads, cs, gs);
    
    IterCount iter(0);
    bool allConverged = false;
    while (!allConverged && iter < this->maxIters()) {
        auto t_iteration = std::chrono::high_resolution_clock::now();
        bool converged = this->onIteration(graph, cs, gs, iter);
        print_duration(t_iteration, "iteration");

        this->onIterationEnd(graph, iter);
        iter++;
    }

    this->closeAlgoKernelServer(algo_kernel_server_threads);

    std::cout<<graph->tid()<<" "<<"Finish all iterations and begin merging vertex data"<<std::endl;

    // Merge vertex data
    uint32_t dstTid = (tileNum + tileIndex - 1) % tileNum;
    uint32_t srcTid = (tileIndex + 1) % tileNum;
    // printf("Here0\n");
    serverTaskComm.sendShareVecVec(gs.remoteVertexSvvs[dstTid], dstTid);
    remoteLocalVertexSvv.shareVecs.clear();
    // printf("Here1\n");
    clientTaskComm.recvShareVecVec(remoteLocalVertexSvv, srcTid);
    // printf("Here2\n");
    this->mergeTwoPartyVertexDataVectorShare(gs, gs.localVertexSvv, remoteLocalVertexSvv);
    // printf("Here3\n");

    serverTaskComm.sendFinish();
    clientTaskComm.recvFinish();

    // this->onAlgoKernelEnd(graph);
    std::cout<<graph->tid()<<" "<<"Finish algo kernel"<<std::endl;
}

template<typename GraphTileType>
void SSEdgeCentricAlgoKernel<GraphTileType>::
onPreprocessClient(Ptr<GraphTileType>& graph, CommSyncType& cs, GraphSummary& gs, bool doOMPreprocess) const { 
    std::vector<Ptr<VertexType>>& localVertexVec = gs.localVertexVec;
    std::vector<std::vector<Ptr<MirrorVertexType>>>& mirrorVertexVecs = gs.mirrorVertexVecs;
    const auto tid = graph->tid();
    TaskComm& clientTaskComm = TaskComm::getClientInstance();
    size_t tileNum = clientTaskComm.getTileNum();
    size_t tileIndex = clientTaskComm.getTileIndex();
    uint64_t maxIters = this->maxIters().cnt();

    if (mirrorVertexVecs.size() == 0) mirrorVertexVecs.resize(tileNum);

    std::cout<<tid<<" "<<"Begin graph preprocessing"<<std::endl;

    // Record src vertices in each vertex / mirror vertex
    for (auto edgeIter = graph->edgeIter(); edgeIter != graph->edgeIterEnd(); ++edgeIter) {
        const auto srcId = edgeIter->srcId();
        const auto dstId = edgeIter->dstId();
        // std::cout<<(uint64_t)srcId<<" "<<(uint64_t)dstId<<std::endl;
        const auto& weight = edgeIter->weight();
        
        if (graph->hasVertex(dstId)) {
            // Local destination.
            auto v = graph->vertex(dstId);
            v->pushToSrcVertexv((uint64_t)srcId);
            v->pushToIncomingEdgev(weight);
            v->pushToIsSrcDummyv(false);
        } else {
            // Remote destination, use mirror vertex.
            auto mv = graph->mirrorVertex(dstId);
            mv->pushToSrcVertexv((uint64_t)srcId);
            mv->pushToIncomingEdgev(weight);
            mv->pushToIsSrcDummyv(false);
        }
    }

    // Construct local vertex Pos vec and update variable src vecs in Scatter
    std::vector<uint64_t>& localVertexPos = gs.localVertexPos;
    std::vector<std::vector<uint64_t>>& mirrorVertexPos = gs.mirrorVertexPos;
    std::vector<std::vector<uint64_t>>& updateSrcVertexPos = gs.updateSrcVertexPos;
    std::vector<std::vector<uint64_t>>& updateDstVertexPos = gs.updateDstVertexPos; // Duplicated
    std::vector<std::vector<bool>>& isUpdateSrcVertexDummy = gs.isUpdateSrcVertexDummy;
    std::vector<std::vector<bool>>& isGatherDstVertexDummy = gs.isGatherDstVertexDummy;
    std::vector<std::vector<EdgeWeightType>>& localEdgeWeightVecs = gs.localEdgeWeightVecs;
    std::vector<uint64_t>& rotationUb = gs.rotationUb;
    std::vector<ShareVecVec>& remoteUpdateSvvs = gs.remoteUpdateSvvs;
    std::vector<ShareVecVec>& localUpdateSvvs = gs.localUpdateSvvs;

    mirrorVertexPos.resize(tileNum);
    updateSrcVertexPos.resize(tileNum);
    updateDstVertexPos.resize(tileNum);
    isUpdateSrcVertexDummy.resize(tileNum);
    isGatherDstVertexDummy.resize(tileNum);
    localEdgeWeightVecs.resize(tileNum);
    rotationUb.resize(tileNum, 0);
    remoteUpdateSvvs.resize(tileNum);
    localUpdateSvvs.resize(tileNum);

    std::vector<uint64_t> maxSrcVSize;
    maxSrcVSize.resize(tileNum, 0);

    if (!clientTaskComm.getIsRotationBased()) {
        for (auto vIter = graph->vertexIter(); vIter != graph->vertexIterEnd(); ++vIter) {
            auto v = vIter->second;
            uint64_t vid = (uint64_t)v->vid();
            std::vector<uint64_t>& cur_srcvv = v->getSrcVertexv();
            std::vector<EdgeWeightType>& cur_incomingEdgev = v->getIncomingEdgev();
            std::vector<bool>& cur_isSrcDummyv = v->getIsSrcDummyv();

            // Add dummy src vertices
            uint64_t dummy_num = get_next_power_of_2(cur_srcvv.size()) - cur_srcvv.size();
            for (int i=0; i<dummy_num; ++i) {
                cur_srcvv.push_back(vid);
                cur_isSrcDummyv.push_back(true);
                cur_incomingEdgev.push_back(-1);
            }
            if (cur_srcvv.size() > maxSrcVSize[tileIndex]) maxSrcVSize[tileIndex] = cur_srcvv.size();
        }
    } else {
        for (auto vIter = graph->vertexIter(); vIter != graph->vertexIterEnd(); ++vIter) {
            auto v = vIter->second;
            uint64_t vid = (uint64_t)v->vid();
            std::vector<uint64_t>& cur_srcvv = v->getSrcVertexv();
            std::vector<EdgeWeightType>& cur_incomingEdgev = v->getIncomingEdgev();
            std::vector<bool>& cur_isSrcDummyv = v->getIsSrcDummyv();

            // Add a dummy src for zero incoming degree vertex
            if (cur_srcvv.size() == 0) {
                cur_srcvv.push_back(vid);
                cur_isSrcDummyv.push_back(true);
                cur_incomingEdgev.push_back(-1);
            }
            if (cur_srcvv.size() > maxSrcVSize[tileIndex]) maxSrcVSize[tileIndex] = cur_srcvv.size();
        }        
    }

    for (auto mvIter = graph->mirrorVertexIter(); mvIter != graph->mirrorVertexIterEnd(); ++mvIter) {
        auto mv = mvIter->second;
        uint64_t mvid = (uint64_t)mv->vid();
        uint64_t mvTid = this->getVertexTid(mvid);
        std::vector<uint64_t>& cur_srcvv = mv->getSrcVertexv();
        std::vector<EdgeWeightType>& cur_incomingEdgev = mv->getIncomingEdgev();
        std::vector<bool>& cur_isSrcDummyv = mv->getIsSrcDummyv();
        if (cur_srcvv.size() == 0) {
            printf("Unexpected mirror vertex with empty src vertex vec!\n");
            exit(-1);
        }     

        // Add dummy src vertex
        uint64_t dummy_num = get_next_power_of_2(cur_srcvv.size()) - cur_srcvv.size();
        for (int i=0; i<dummy_num; ++i) {
            cur_srcvv.push_back(cur_srcvv[0]);
            cur_isSrcDummyv.push_back(true);
            cur_incomingEdgev.push_back(-1);
        }
        if (cur_srcvv.size() > maxSrcVSize[mvTid]) maxSrcVSize[mvTid] = cur_srcvv.size();
    }

    // Counting sort
    std::vector<std::vector<std::vector<uint64_t>>> countingVec;
    countingVec.resize(tileNum);
    for (int i=0; i<tileNum; ++i) {
        countingVec[i].resize(maxSrcVSize[i]);
    }
    for (auto vIter = graph->vertexIter(); vIter != graph->vertexIterEnd(); ++vIter) {
        auto v = vIter->second;
        uint64_t vid = (uint64_t)vIter->second->vid();
        std::vector<uint64_t>& cur_srcvv = v->getSrcVertexv();
        if (cur_srcvv.size() > 0) {
            countingVec[tileIndex][cur_srcvv.size()-1].push_back(vid);
        }
    }
    for (auto mvIter = graph->mirrorVertexIter(); mvIter != graph->mirrorVertexIterEnd(); ++mvIter) {
        auto mv = mvIter->second;
        uint64_t mvid = (uint64_t)mvIter->second->vid();
        uint64_t mvTid = this->getVertexTid(mvid);
        std::vector<uint64_t>& cur_srcvv = mv->getSrcVertexv();
        if (cur_srcvv.size() > 0) {
            countingVec[mvTid][cur_srcvv.size()-1].push_back(mvid);
        }
    }

    // Update src vertex pos vec and isDummy vec construction based on counting sort result
    for (int i=0; i<tileNum; ++i) {
        for (int j=countingVec[i].size()-1; j>=0; --j) {
            if (countingVec[i][j].size() > 0) {
                if (j+1 > rotationUb[i]) rotationUb[i] = j+1; 
                for (int k=0; k<countingVec[i][j].size(); ++k) {
                    // Get vertex or mirror vertex
                    const auto dstId = countingVec[i][j][k];
                    if (tileIndex == i) {
                        auto v = graph->vertex(dstId);
                        v->setReorderedIndex(localVertexPos.size());
                        localVertexPos.push_back(dstId);
                        localVertexVec.push_back(v);
                        std::vector<uint64_t>& cur_srcvv = v->getSrcVertexv();
                        std::vector<bool>& cur_isSrcDummyv = v->getIsSrcDummyv();
                        std::vector<EdgeWeightType>& cur_incomingEdgev = v->getIncomingEdgev();
                        updateSrcVertexPos[i].insert(updateSrcVertexPos[i].end(), cur_srcvv.begin(), cur_srcvv.end());
                        updateDstVertexPos[i].insert(updateDstVertexPos[i].end(), cur_srcvv.size(), dstId);
                        isUpdateSrcVertexDummy[i].insert(isUpdateSrcVertexDummy[i].end(), cur_isSrcDummyv.begin(), cur_isSrcDummyv.end());
                        localEdgeWeightVecs[i].insert(localEdgeWeightVecs[i].end(), cur_incomingEdgev.begin(), cur_incomingEdgev.end());
                        isGatherDstVertexDummy[i].push_back(cur_isSrcDummyv[0]);
                    } else {
                        auto mv = graph->mirrorVertex(dstId);
                        uint64_t mvTid = this->getVertexTid(dstId);
                        mirrorVertexPos[mvTid].push_back(dstId);
                        mirrorVertexVecs[i].push_back(mv);
                        std::vector<uint64_t>& cur_srcvv = mv->getSrcVertexv();
                        std::vector<bool>& cur_isSrcDummyv = mv->getIsSrcDummyv();
                        std::vector<EdgeWeightType>& cur_incomingEdgev = mv->getIncomingEdgev();
                        updateSrcVertexPos[i].insert(updateSrcVertexPos[i].end(), cur_srcvv.begin(), cur_srcvv.end());
                        updateDstVertexPos[i].insert(updateDstVertexPos[i].end(), cur_srcvv.size(), dstId);
                        isUpdateSrcVertexDummy[i].insert(isUpdateSrcVertexDummy[i].end(), cur_isSrcDummyv.begin(), cur_isSrcDummyv.end());
                        localEdgeWeightVecs[i].insert(localEdgeWeightVecs[i].end(), cur_incomingEdgev.begin(), cur_incomingEdgev.end());
                    }
                }
            }
        }
    }

    // Send mirrorVertexPos to the target party (use cs to avoid channel conflicts in TaskComm)
    for (int i=0; i<tileNum; ++i) {
        if (i != tileIndex) {
            PosVec tmpPosVec;
            tmpPosVec.pos.swap(mirrorVertexPos[i]);
            cs.sendPosVec(tmpPosVec, tileIndex, i);
            tmpPosVec.pos.swap(mirrorVertexPos[i]);
        }
    }

    // Receive remote mirror vertex pos from the target party
    std::vector<std::vector<uint64_t>>& remoteMirrorVertexPos = gs.remoteMirrorVertexPos;
    remoteMirrorVertexPos.resize(tileNum);
    for (int i=0; i<tileNum; ++i) {
        if (i != tileIndex) {
            PosVec tmpPosVec;
            cs.recvPosVec(tmpPosVec, i, tileIndex);
            tmpPosVec.pos.swap(remoteMirrorVertexPos[i]);
            isGatherDstVertexDummy[i].resize(localVertexPos.size(), true);
            for (int j=0; j<remoteMirrorVertexPos[i].size(); ++j) {
                auto v = graph->vertex(remoteMirrorVertexPos[i][j]);
                isGatherDstVertexDummy[i][v->getReorderedIndex()] = false;
            }
        }
    }

    if (!doOMPreprocess) return;

    uint32_t plainNumPerOperand = getPlainNumPerOperand();

    auto t_preprocess_OM = std::chrono::high_resolution_clock::now();
    // Scatter & Pre-merge for taskvs
    std::cout<<tid<<" "<<"Begin preprocessing oblivious mapper"<<std::endl;
    std::vector<std::thread> threads;
    for (int i=0; i<tileNum; ++i) {
        if (i != tileIndex) {
            threads.emplace_back([this, tileIndex, tileNum, i, plainNumPerOperand, maxIters, &localVertexPos, &updateSrcVertexPos, &updateDstVertexPos, &mirrorVertexPos, &remoteMirrorVertexPos]() {
                uint32_t iter = 0;
                uint64_t batchSize = 0;
                for (iter=0; iter<maxIters; iter+=batchSize) {
                    std::cout<<tileIndex<<" "<<"Preprocessing oblivious mapper, iter "<<iter<<std::endl;
                    uint32_t preprocessId = 0;

                    // Oblivious Mapper
                    // srcPos: local vertex
                    // dstPos: update src vertex for local vertex
                    std::cout<<tileIndex<<" "<<"OM "<<0<<std::endl;
                    if (i == (tileIndex + 1) % tileNum) {
                        client_batch_oblivious_mapper_preprocess(localVertexPos, updateSrcVertexPos[tileIndex], plainNumPerOperand, iter, preprocessId, i);
                        preprocessId += 1;
                    }

                    // Oblivious Mapper
                    // srcPos: local vertex
                    // dstPos: update src vertex for mirror vertex
                    std::cout<<tileIndex<<" "<<"OM "<<1<<" "<<i<<std::endl;
                    // for (auto x : updateSrcVertexPos[i]) {
                    //     std::cout<<">>"<<x<<std::endl;
                    // }
                    batchSize = client_batch_oblivious_mapper_preprocess(localVertexPos, updateSrcVertexPos[i], plainNumPerOperand, iter, preprocessId, i);
                    preprocessId += 1;

                    // Oblivious Mapper
                    // srcPos: duplicated update dst vertex for local vertex
                    // dstPos: deduplicated update dst vertex for local vertex (same as local vertex)
                    std::cout<<tileIndex<<" "<<"OM "<<2<<" "<<i<<std::endl;
                    if (i == (tileIndex + 1) % tileNum) {
                        client_batch_oblivious_mapper_preprocess(updateDstVertexPos[tileIndex], localVertexPos, plainNumPerOperand, iter, preprocessId, i);
                        preprocessId += 1;
                    }

                    // Oblivious Mapper
                    // srcPos: duplicated update dst vertex for mirror vertex
                    // dstPos: deduplicated update dst vertex for mirror vertex (same as mirror vertex)
                    std::cout<<tileIndex<<" "<<"OM "<<3<<" "<<i<<std::endl;
                    client_batch_oblivious_mapper_preprocess(updateDstVertexPos[i], mirrorVertexPos[i], plainNumPerOperand, iter, preprocessId, i);
                    preprocessId += 1;

                    // Oblivious Mapper
                    // srcPos: remote mirror vertex
                    // dstPos: local vertex
                    std::cout<<tileIndex<<" "<<"OM "<<4<<" "<<i<<std::endl;
                    client_batch_oblivious_mapper_preprocess(remoteMirrorVertexPos[i], localVertexPos, plainNumPerOperand, iter, preprocessId, i, true);
                    preprocessId += 1;

                    std::cout<<tileIndex<<" "<<"OM "<<"end of iteration"<<" "<<i<<std::endl;
                }
            });
        }
    }

    for (auto& thrd : threads)
        thrd.join();
    
    print_duration(t_preprocess_OM, "preprocess_OM");
}

template<typename GraphTileType>
void SSEdgeCentricAlgoKernel<GraphTileType>::
onPreprocessServer(std::vector<std::thread>& threads, bool doOMPreprocess) const {
    if (!doOMPreprocess) return;
	TaskComm& serverTaskComm = TaskComm::getServerInstance();
	size_t tileNum = serverTaskComm.getTileNum();
	size_t tileIndex = serverTaskComm.getTileIndex();
    uint64_t maxIters = this->maxIters().cnt();

	for (int i = 0; i < tileNum; i++) {
		if (i != tileIndex) {
			threads.emplace_back([this, i, &serverTaskComm, tileIndex, tileNum, maxIters]() {
                uint32_t iter = 0;
                uint64_t batchSize = 0;
                for (iter=0; iter<maxIters; iter+=batchSize) {
                    uint32_t preprocessId = 0;

                    // Oblivious Mapper
                    // srcPos: local vertex
                    // dstPos: update src vertex for local vertex
                    if ((i + 1) % tileNum == tileIndex) {
                        server_batch_oblivious_mapper_preprocess(iter, preprocessId, i);
                        preprocessId += 1;
                    }

                    // Oblivious Mapper
                    // srcPos: local vertex
                    // dstPos: update src vertex for mirror vertex
                    batchSize = server_batch_oblivious_mapper_preprocess(iter, preprocessId, i);
                    preprocessId += 1;

                    // Oblivious Mapper
                    // srcPos: duplicated update dst vertex for local vertex
                    // dstPos: deduplicated update dst vertex for local vertex (same as local vertex)
                    if ((i + 1) % tileNum == tileIndex) {
                        server_batch_oblivious_mapper_preprocess(iter, preprocessId, i);
                        preprocessId += 1;
                    }

                    // Oblivious Mapper
                    // srcPos: duplicated update dst vertex for mirror vertex
                    // dstPos: deduplicated update dst vertex for mirror vertex (same as mirror vertex)
                    server_batch_oblivious_mapper_preprocess(iter, preprocessId, i);
                    preprocessId += 1;

                    // Oblivious Mapper
                    // srcPos: remote mirror vertex
                    // dstPos: local vertex
                    server_batch_oblivious_mapper_preprocess(iter, preprocessId, i);
                    preprocessId += 1;

                    // std::cout<<"Server final preprocessId "<<preprocessId<<" "<<i<<" "<<(tileIndex - 1) % tileNum<<std::endl;
                }            
			});
		}
	}
}

template<typename GraphTileType>
bool SSEdgeCentricAlgoKernel<GraphTileType>::
onIteration(Ptr<GraphTileType>& graph, CommSyncType& cs, GraphSummary& gs, const IterCount& iter) const {
    printf("tid-> %lld, iteration-> %lld\n", graph->tid(), iter.cnt());    
    const auto tid = graph->tid();
    TaskComm& clientTaskComm = TaskComm::getClientInstance();
    TaskComm& serverTaskComm = TaskComm::getServerInstance();
    size_t tileNum = clientTaskComm.getTileNum();
    size_t tileIndex = clientTaskComm.getTileIndex();

    uint32_t plainNumPerOperand = getPlainNumPerOperand();

    std::cout<<tid<<" "<<"Begin Scatter task generation"<<std::endl;
    
    std::vector<ShareVecVec> updateSrcs;
    updateSrcs.resize(tileNum);
    std::vector<std::thread> threads;
    bar_t barrier(tileNum - 1);
    for (int i=0; i<tileNum; ++i) {
        if (i != tileIndex) {
            threads.emplace_back([this, i, tileIndex, tileNum, plainNumPerOperand, &gs, &updateSrcs, &clientTaskComm, &serverTaskComm, &iter, &barrier, &graph](){
                uint32_t preprocessId = 0;

                auto t_Scatter_preparation = std::chrono::high_resolution_clock::now();
                std::cout<<tileIndex<<" "<<"Begin update src mapping with "<<i<<std::endl;

                if (i == (tileIndex + 1) % tileNum) {
                    client_oblivious_mapper_online(gs.localVertexPos, gs.updateSrcVertexPos[tileIndex],
                                        gs.localVertexSvv, updateSrcs[tileIndex], plainNumPerOperand, 
                                        iter.cnt(), preprocessId, i);
                    preprocessId += 1;
                }

                // printf("H1\n");

                client_oblivious_mapper_online(gs.localVertexPos, gs.updateSrcVertexPos[i],
                                    gs.localVertexSvv, updateSrcs[i], plainNumPerOperand, 
                                    iter.cnt(), preprocessId, i);
                preprocessId += 1;

                print_duration(t_Scatter_preparation, "Scatter_preparation");

                // if (i == 1) {
                //     for (int kk=0; kk<gs.localVertexPos.size(); ++kk) {
                //         if (gs.localVertexPos[kk] == 0) {
                //             std::cout<<"Zero pos "<<kk<<std::endl;
                //             auto& x = gs.localVertexSvv.shareVecs[kk];
                //             std::cout<<x.shares[0]<<" "<<x.shares[1]<<"  "<<std::endl;
                //         }
                //     }
                //     for (auto& x : gs.updateSrcVertexPos[i]) std::cout<<x<<" ";
                //     std::cout<<std::endl;
                //     for (auto& x : updateSrcs[i].shareVecs) std::cout<<x.shares[0]<<" "<<x.shares[1]<<"  ";
                //     std::cout<<std::endl;
                // }

                // Scatter and pre-merging for the mirror vertex of party i
                TaskqHandlerConfig& thc = clientTaskComm.getTaskqHandlerConfig(i);
                thc.sendOperand = false;
                thc.mergeResult = false;
                thc.sendTaskqDigest = false;
                std::vector<Task> taskv;

                auto clientComputeUpdate = [this, &gs, &thc, &clientTaskComm, &taskv, &graph, i, tileIndex](ShareVecVec& updateSrc, ShareVecVec& duplicatedUpdateSvv, const uint32_t dstTid) {
                    set_up_mpc_channel(true, i);
                    auto t_Scatter_computation = std::chrono::high_resolution_clock::now();

                    uint64_t scatterTaskNum = updateSrc.shareVecs.size();
                    if (scatterTaskNum != gs.localEdgeWeightVecs[dstTid].size()) {
                        printf("Unmatched src num and edge weight num during Scatter!\n");
                        exit(-1);
                    }

                    // if (i==1 && dstTid == 1) {
                    //     for (uint64_t j=0; j<scatterTaskNum; ++j) {
                    //         std::cout<<">> Client Scatter inputs"<<std::endl;
                    //         std::cout<<updateSrc.shareVecs[j].shares[0]<<" "<<updateSrc.shareVecs[j].shares[1]<<" "<<gs.localEdgeWeightVecs[dstTid][j]<<" "<<gs.updateSrcVertexPos[dstTid][j]<<std::endl;
                    //     }                        
                    // }

                    for (uint64_t j=0; j<scatterTaskNum; ++j) 
                        taskv.push_back(this->genScatterTask(graph->vertex(gs.updateSrcVertexPos[dstTid][j]), updateSrc.shareVecs[j], gs.localEdgeWeightVecs[dstTid][j], gs.updateSrcVertexPos[dstTid][j], tileIndex, gs.updateDstVertexPos[dstTid][j], dstTid, gs.isUpdateSrcVertexDummy[dstTid][j]));
                        
                    // Compute Scatter task vector 
                    thc.rotation = 0;
                    client_task_vector_direct_handler(taskv, i);

                    print_duration(t_Scatter_computation, "Scatter_computation");

                    auto t_premerging = std::chrono::high_resolution_clock::now();

                    // Pre-merging
                    this->fromScatterTaskvResultToPreMergingTaskv(taskv);

                    // Send rotation ub
                    // uint64_t rotation_ub = gs.rotationUb[dstTid];
                    uint64_t rotation_ub = taskv.size();
                    clientTaskComm.sendRotationUB(rotation_ub, i);

                    uint64_t clientPreprocessComm = 0;

                    thc.rotation = 1;
                    while (thc.rotation < rotation_ub) {
                        // if (i==1 && dstTid == 1) {
                        //     std::cout<<"Client pre-merging rotation before "<<thc.rotation<<" taskv size "<<taskv.size()<<std::endl;
                        //     for (Task& task : taskv) {
                        //         TaskPayload* tp = (TaskPayload*)(task.buf);
                        //         ShareVec operand0 = tp->getOperandShare(0);
                        //         ShareVec operand1 = tp->getOperandShare(1);
                        //         std::cout<<task.srcIndex<<" "<<task.vertexIndex<<" "<<task.isDummy<<" "<<operand0.shares[0]<<" "<<operand0.shares[1]<<" "<<operand1.shares[0]<<" "<<operand1.shares[1]<<std::endl;
                        //     }
                        // }
                        clientPreprocessComm += client_task_vector_direct_handler(taskv, i);

                        // if (i==1 && dstTid == 1) {
                        //     std::cout<<"Client pre-merging rotation after "<<thc.rotation<<" taskv size "<<taskv.size()<<std::endl;
                        //     for (Task& task : taskv) {
                        //         TaskPayload* tp = (TaskPayload*)(task.buf);
                        //         ShareVec operand0 = tp->getOperandShare(0);
                        //         ShareVec operand1 = tp->getOperandShare(1);
                        //         std::cout<<task.srcIndex<<" "<<task.vertexIndex<<" "<<task.isDummy<<" "<<operand0.shares[0]<<" "<<operand0.shares[1]<<" "<<operand1.shares[0]<<" "<<operand1.shares[1]<<std::endl;
                        //     }
                        // }

                        thc.rotation = thc.rotation << 1;
                    }

                    print_comm(clientPreprocessComm, "clientPreprocessComm");

                    if (rotation_ub <= 1) {
                        for (auto& task : taskv) task.finished = true;
                    }
                    
                    duplicatedUpdateSvv.shareVecs.clear();
                    duplicatedUpdateSvv.shareVecs.resize(scatterTaskNum);
                    for (uint64_t j=0; j<scatterTaskNum; ++j) {
                        UpdateType::getTaskResult(taskv[j], duplicatedUpdateSvv.shareVecs[j]);
                        taskv[j].delete_task_content_buf();
                    }

                    // if (i==1 && dstTid == 1) {
                    //     for (uint64_t j=0; j<duplicatedUpdateSvv.shareVecs.size(); ++j) {
                    //         std::cout<<"Scatter results"<<std::endl;
                    //         std::cout<<gs.updateDstVertexPos[1][j]<<" "<<duplicatedUpdateSvv.shareVecs[j].shares[0]<<" "<<duplicatedUpdateSvv.shareVecs[j].shares[1]<<std::endl;
                    //     }                        
                    // }
                    taskv.clear();

                    print_duration(t_premerging, "premerging");
                    close_mpc_channel(true, i);                   
                };

                std::cout<<tileIndex<<" "<<"Begin update extraction mapping with "<<i<<std::endl;
                ShareVecVec duplicatedUpdateSvv;
                if (i == (tileIndex + 1) % tileNum) {
                    clientComputeUpdate(updateSrcs[tileIndex], duplicatedUpdateSvv, tileIndex);
                    auto t_premerged_extraction = std::chrono::high_resolution_clock::now();
                    client_oblivious_mapper_online(gs.updateDstVertexPos[tileIndex], gs.localVertexPos,
                                        duplicatedUpdateSvv, gs.localUpdateSvvs[tileIndex], plainNumPerOperand, 
                                        iter.cnt(), preprocessId, i);
                    preprocessId += 1;
                    print_duration(t_premerged_extraction, "premerged_extraction");
                    duplicatedUpdateSvv.shareVecs.clear();
                }

                // printf("H2\n");

                clientComputeUpdate(updateSrcs[i], duplicatedUpdateSvv, i);
                auto t_premerged_extraction = std::chrono::high_resolution_clock::now();
                client_oblivious_mapper_online(gs.updateDstVertexPos[i], gs.mirrorVertexPos[i],
                                    duplicatedUpdateSvv, gs.remoteUpdateSvvs[i], plainNumPerOperand, 
                                    iter.cnt(), preprocessId, i);
                preprocessId += 1;
                print_duration(t_premerged_extraction, "premerged_extraction");

                printf("client here\n");
                Semaphore& local_update_ready_smp = clientTaskComm.getLocalUpdateReadySmp(i);
                Semaphore& remote_update_ready_smp = serverTaskComm.getRemoteUpdateReadySmp(i);
                remote_update_ready_smp.release();
                local_update_ready_smp.acquire();
                printf("client after\n");

                auto t_Gather_preparation = std::chrono::high_resolution_clock::now();

                std::cout<<tileIndex<<" "<<"Begin update extension mapping with "<<i<<std::endl;

                ShareVecVec tmpUpdateSvv;
                client_oblivious_mapper_online(gs.remoteMirrorVertexPos[i], gs.localVertexPos,
                                    gs.localUpdateSvvs[i], tmpUpdateSvv, plainNumPerOperand, 
                                    iter.cnt(), preprocessId, i, true);
                preprocessId += 1;

                gs.localUpdateSvvs[i].shareVecs.clear();
                gs.localUpdateSvvs[i].shareVecs.swap(tmpUpdateSvv.shareVecs);

                print_duration(t_Gather_preparation, "Gather_preparation");

                // Synchronize across all threads
                barrier.wait();

                auto t_Gather_computation = std::chrono::high_resolution_clock::now();

                std::cout<<tileIndex<<" "<<"Begin gather computation mapping with "<<i<<std::endl;
                // Gather
                if (i == (tileIndex + 1) % tileNum) {
                    set_up_mpc_channel(true, i);
                    for (int j=0; j<tileNum; ++j) {
                        uint64_t gatherTaskNum = gs.localVertexSvv.shareVecs.size();
                        if (gatherTaskNum != gs.localUpdateSvvs[j].shareVecs.size()) {
                            printf("Client Unmatched update num and vertex num during cooperation! %lld %d\n", gatherTaskNum, gs.localUpdateSvvs[j].shareVecs.size());
                            exit(-1);
                        }
                        for (int m=0; m<gatherTaskNum; ++m) {
                            taskv.push_back(this->genGatherTask(gs.localVertexSvv.shareVecs[m], gs.localUpdateSvvs[j].shareVecs[m], gs.localVertexPos[m], tileIndex, gs.isGatherDstVertexDummy[j][m]));
                        }

                        // if (j == 1 && i == 0)
                        //     for (uint64_t m=0; m<gatherTaskNum; ++m) {
                        //         std::cout<<"1Here "<<m<<" "<<gs.localVertexPos[m]<<" "<<gs.localUpdateSvvs[j].shareVecs[m].shares[0]<<" "<<gs.localUpdateSvvs[j].shareVecs[m].shares[1]<<" "<<gs.isGatherDstVertexDummy[j][m]<<std::endl;
                        //     }

                        thc.rotation = 0;
                        client_task_vector_direct_handler(taskv, i);
                        
                        for (int m=0; m<gatherTaskNum; ++m) {
                            this->writeGatherTaskResult(taskv[m], gs.localVertexSvv.shareVecs[m]);
                            // std::cout<<"gs.localVertexSvv.shareVecs[m].shares[0]: "<<gs.localVertexSvv.shareVecs[m].shares[0]<<std::endl;
                            taskv[m].delete_task_content_buf();
                        }
                        taskv.clear();

                        // if (j == 1 && i == 0)
                        //     for (uint64_t m=0; m<gatherTaskNum; ++m) {
                        //         std::cout<<"2Here "<<m<<" "<<gs.localVertexPos[m]<<" "<<gs.localVertexSvv.shareVecs[m].shares[0]<<" "<<gs.localVertexSvv.shareVecs[m].shares[1]<<" "<<gs.isGatherDstVertexDummy[j][m]<<std::endl;
                        //     }
                    }
                    close_mpc_channel(true, i);                    
                }

                print_duration(t_Gather_computation, "Gather_computation");
            });
        }
    }

    for (auto& thrd : threads)
        thrd.join();

    return false;
}

template<typename GraphTileType>
void SSEdgeCentricAlgoKernel<GraphTileType>::runAlgoKernelServer(std::vector<std::thread>& threads, CommSyncType& cs, GraphSummary& gs) const {
	TaskComm& clientTaskComm = TaskComm::getClientInstance();
    TaskComm& serverTaskComm = TaskComm::getServerInstance();
	size_t tileNum = serverTaskComm.getTileNum();
	size_t tileIndex = serverTaskComm.getTileIndex();

    uint64_t maxIters = this->maxIters().cnt();

    bar_t barrier(tileNum - 1);
	for (int i = 0; i < tileNum; i++) {
		if (i != tileIndex) {
			threads.emplace_back([this, i, &serverTaskComm, &clientTaskComm, &cs, &gs, &barrier, tileIndex, tileNum, maxIters]() {
                uint64_t iter = 0;
                TaskqHandlerConfig& thc = serverTaskComm.getTaskqHandlerConfig(i);
                thc.sendOperand = false;
                thc.mergeResult = false;
                thc.sendTaskqDigest = false;
                std::vector<Task>& taskv = serverTaskComm.getTaskv(i);
                while (iter < maxIters) { // On iteration
                    uint32_t preprocessId = 0;
                    printf("Server iter %lld\n", iter);
                    ShareVecVec coUpdateSrc;
                    ShareVecVec updateSrc;

                    if (tileIndex == (i + 1) % tileNum) {
                        server_oblivious_mapper_online(gs.remoteVertexSvvs[i], coUpdateSrc,
                                        iter, preprocessId, i);
                        preprocessId += 1;
                    }

                    server_oblivious_mapper_online(gs.remoteVertexSvvs[i], updateSrc,
                                    iter, preprocessId, i);
                    preprocessId += 1;
                    // if (i == 0) {
                    //     auto& x = gs.remoteVertexSvvs[i].shareVecs[11750];
                    //     std::cout<<x.shares[0]<<" "<<x.shares[1]<<"  "<<std::endl;;
                    //     for (auto& x : updateSrc.shareVecs) std::cout<<x.shares[0]<<" "<<x.shares[1]<<"  ";
                    //     std::cout<<std::endl;
                    // }

                    auto serverComputeUpdate = [this, &gs, &thc, &serverTaskComm, &taskv, i, tileIndex](ShareVecVec& updateSrc, ShareVecVec& duplicatedUpdateSvv) {
                        set_up_mpc_channel(false, i);
                        // Scatter
                        thc.rotation = 0;
                        uint64_t scatterTaskNum = updateSrc.shareVecs.size();

                        // if (i==0 && scatterTaskNum == 2) {
                        //     for (uint64_t j=0; j<scatterTaskNum; ++j) {
                        //         std::cout<<">> Server Scatter inputs"<<std::endl;
                        //         std::cout<<" "<<updateSrc.shareVecs[j].shares[0]<<" "<<updateSrc.shareVecs[j].shares[1]<<std::endl;
                        //     }                        
                        // }

                        for (uint64_t j=0; j<scatterTaskNum; ++j) {
                            taskv.push_back(this->genScatterTask(nullptr, updateSrc.shareVecs[j], -1, -1, -1, -1, -1));
                        }
                        server_task_vector_direct_handler(taskv, i);

                        // Pre-merge
                        this->fromScatterTaskvResultToPreMergingTaskv(taskv);

                        // Receive rotation ub
                        serverTaskComm.recvEntrance(i);
                        // uint64_t rotation_ub = thc.rotation;
                        uint64_t rotation_ub = taskv.size();

                        uint64_t serverPreprocessComm = 0;

                        thc.rotation = 1;
                        while (thc.rotation < rotation_ub) {
                            // if (i==0 && taskv.size()==7) {
                            //     std::cout<<"server pre-merging rotation before "<<thc.rotation<<" taskv size "<<taskv.size()<<std::endl;
                            //     for (Task& task : taskv) {
                            //         TaskPayload* tp = (TaskPayload*)(task.buf);
                            //         ShareVec operand0 = tp->getOperandShare(0);
                            //         ShareVec operand1 = tp->getOperandShare(1);
                            //         std::cout<<operand0.shares[0]<<" "<<operand0.shares[1]<<" "<<operand1.shares[0]<<" "<<operand1.shares[1]<<std::endl;
                            //     }
                            // }

                            serverPreprocessComm += server_task_vector_direct_handler(taskv, i);
                            thc.rotation = thc.rotation << 1;

                            // if (i==0 && taskv.size()==7) {
                            //     std::cout<<"server pre-merging rotation after "<<thc.rotation<<" taskv size "<<taskv.size()<<std::endl;
                            //     for (Task& task : taskv) {
                            //         TaskPayload* tp = (TaskPayload*)(task.buf);
                            //         ShareVec operand0 = tp->getOperandShare(0);
                            //         ShareVec operand1 = tp->getOperandShare(1);
                            //         std::cout<<operand0.shares[0]<<" "<<operand0.shares[1]<<" "<<operand1.shares[0]<<" "<<operand1.shares[1]<<std::endl;
                            //     }
                            // }
                        }

                        print_comm(serverPreprocessComm, "serverPreprocessComm");

                        if (rotation_ub <= 1) {
                            for (auto& task : taskv) task.finished = true;
                        }

                        duplicatedUpdateSvv.shareVecs.clear();
                        duplicatedUpdateSvv.shareVecs.resize(scatterTaskNum);
                        for (uint64_t j=0; j<scatterTaskNum; ++j) {
                            UpdateType::getTaskResult(taskv[j], duplicatedUpdateSvv.shareVecs[j]);
                            taskv[j].delete_task_content_buf();
                        }

                        // if (i==0 && duplicatedUpdateSvv.shareVecs.size() == 2) {
                        //     for (uint64_t j=0; j<duplicatedUpdateSvv.shareVecs.size(); ++j) {
                        //         std::cout<<"Server Scatter results"<<std::endl;
                        //         std::cout<<" "<<duplicatedUpdateSvv.shareVecs[j].shares[0]<<" "<<duplicatedUpdateSvv.shareVecs[j].shares[1]<<std::endl;
                        //     }                        
                        // }

                        taskv.clear();
                        close_mpc_channel(false, i);
                    };

                    ShareVecVec duplicatedUpdateSvv;
                    if (tileIndex == (i + 1) % tileNum) {
                        serverComputeUpdate(coUpdateSrc, duplicatedUpdateSvv);
                        // printf("++++%d\n", gs.remoteUpdateSvvs.size());
                        server_oblivious_mapper_online(duplicatedUpdateSvv, gs.remoteUpdateSvvs[tileIndex], 
                                            iter, preprocessId, i);
                        preprocessId += 1;
                        duplicatedUpdateSvv.shareVecs.clear();
                    }

                    serverComputeUpdate(updateSrc, duplicatedUpdateSvv);
                    server_oblivious_mapper_online(duplicatedUpdateSvv, gs.localUpdateSvvs[i], 
                                        iter, preprocessId, i);
                    preprocessId += 1;

                    printf("Server here\n");

                    Semaphore& local_update_ready_smp = clientTaskComm.getLocalUpdateReadySmp(i);
                    Semaphore& remote_update_ready_smp = serverTaskComm.getRemoteUpdateReadySmp(i);
                    local_update_ready_smp.release();
                    remote_update_ready_smp.acquire();

                    printf("server after\n");

                    ShareVecVec tmpUpdateSvv;
                    server_oblivious_mapper_online(gs.remoteUpdateSvvs[i], tmpUpdateSvv,
                                        iter, preprocessId, i);
                    preprocessId += 1;

                    gs.remoteUpdateSvvs[i].shareVecs.clear();
                    gs.remoteUpdateSvvs[i].shareVecs.swap(tmpUpdateSvv.shareVecs);

                    // Synchronize across all threads
                    barrier.wait();

                    std::cout<<"Barrier, "<<tileIndex<<" Server, "<<"iter: "<<iter<<" "<<i<<std::endl;

                    std::vector<ShareVecVec> remoteUpdateSvvs;

                    if (tileIndex != (i + 1) % tileNum) {
                        cs.sendShareVecVec(gs.remoteUpdateSvvs[i], tileIndex, (i + 1) % tileNum);
                    } else {
                        remoteUpdateSvvs.resize(tileNum);
                        remoteUpdateSvvs[tileIndex].shareVecs.swap(gs.remoteUpdateSvvs[i].shareVecs);
                        remoteUpdateSvvs[i].shareVecs.swap(gs.remoteUpdateSvvs[tileIndex].shareVecs);
                        for (int j=0; j<tileNum; ++j) {
                            if (j != tileIndex && j != i) {
                                cs.recvShareVecVec(remoteUpdateSvvs[j], j, tileIndex);
                            }
                        }
                    }

                    barrier.wait();

                    if (tileIndex != (i + 1) % tileNum) {
                        // serverTaskComm.sendShareVecVec(gs.remoteUpdateSvvs[i], (i + 1) % tileNum);
                        // std::cout<<tileIndex<<" Server delegates "<<gs.remoteUpdateSvvs[i].shareVecs.size()<<" to "<<(i + 1) % tileNum<<std::endl;
                        // std::cout<<gs.remoteUpdateSvvs[i].shareVecs[0].shares[0]<<" ++ "<<gs.remoteUpdateSvvs[i].shareVecs[0].shares[1]<<std::endl;
                        // clientTaskComm.recvShareVecVec(gs.remoteVertexSvvs[i], (i + 1) % tileNum);
                        // std::cout<<"Received Gather Taskv result, "<<tileIndex<<" Server, "<<"iter: "<<iter<<" "<<i<<std::endl;
                        cs.recvShareVecVec(gs.remoteVertexSvvs[i], (i + 1) % tileNum, tileIndex);
                    } else {
                        std::cout<<"Compute Gather Taskv, "<<tileIndex<<" Server, "<<"iter: "<<iter<<" "<<i<<std::endl;
                        set_up_mpc_channel(false, i);
                        for (int j=0; j<tileNum; ++j) {
                            uint64_t gatherTaskNum = remoteUpdateSvvs[j].shareVecs.size();
                            if (gatherTaskNum != gs.remoteVertexSvvs[i].shareVecs.size()) {
                                printf("%d Server Unmatched update num and vertex num during cooperation! %lld %d %d %d\n", tileIndex, gatherTaskNum, gs.remoteVertexSvvs[i].shareVecs.size(), i, j);
                                exit(-1);
                            }
                            for (int m=0; m<gatherTaskNum; ++m) {
                                taskv.push_back(this->genGatherTask(gs.remoteVertexSvvs[i].shareVecs[m], remoteUpdateSvvs[j].shareVecs[m], -1, i));
                            }

                            // if (j == 1 && i == 1) {
                            //     // uint64_t m=2718;
                            //     // std::cout<<"1THere "<<m<<" "<<remoteUpdateSvvs[j].shareVecs[m].shares[0]<<" "<<remoteUpdateSvvs[j].shareVecs[m].shares[1]<<std::endl;
                            //     // m=8796;
                            //     for (int m=0; m<gatherTaskNum; ++m)
                            //         std::cout<<"1THere "<<m<<" "<<remoteUpdateSvvs[j].shareVecs[m].shares[0]<<" "<<remoteUpdateSvvs[j].shareVecs[m].shares[1]<<std::endl;
                            // }                                                            

                            thc.rotation = 0;
                            server_task_vector_direct_handler(taskv, i);
                            
                            for (int m=0; m<gatherTaskNum; ++m) {
                                this->writeGatherTaskResult(taskv[m], gs.remoteVertexSvvs[i].shareVecs[m]);
                                taskv[m].delete_task_content_buf();
                            }
                            taskv.clear();

                            // if (j == 1 && i == 1) {
                            //     // uint64_t m=2718;
                            //     // std::cout<<"2THere "<<m<<" "<<gs.remoteVertexSvvs[i].shareVecs[m].shares[0]<<" "<<gs.remoteVertexSvvs[i].shareVecs[m].shares[1]<<std::endl;
                            //     // m=8796;
                            //     for (int m=0; m<gatherTaskNum; ++m)
                            //         std::cout<<"2THere "<<m<<" "<<gs.remoteVertexSvvs[i].shareVecs[m].shares[0]<<" "<<gs.remoteVertexSvvs[i].shareVecs[m].shares[1]<<std::endl;
                            // } 
                        }

                        // Send vertex data share to the other party
                        for (int j=0; j<tileNum; ++j) {
                            if (j != tileIndex && j != i) {
                                // serverTaskComm.sendShareVecVec(gs.remoteVertexSvvs[i], j);
                                cs.sendShareVecVec(gs.remoteVertexSvvs[i], tileIndex, j);
                            }
                        }
                        close_mpc_channel(false, i);
                    }

                    barrier.wait();

                    std::cout<<"Finish iteration, "<<tileIndex<<" Server, "<<"iter: "<<iter<<" "<<i<<std::endl;

                    // if (tileIndex != (i + 1) % tileNum) {
                    //     serverTaskComm.recvShareVecVec(gs.remoteVertexSvvs[i], (i + 1) % tileNum);
                    //     std::cout<<"Received Gather Taskv result, "<<tileIndex<<" Server, "<<"iter: "<<iter<<" "<<i<<std::endl;                        
                    // } else {
                    //     // Send vertex data share to the other party
                    //     for (int j=0; j<tileNum; ++j) {
                    //         if (j != tileIndex && j != i) {
                    //             clientTaskComm.sendShareVecVec(gs.remoteVertexSvvs[i], j);
                    //         }
                    //     }
                    // }

                    // barrier.wait();

                    iter++;
                }
			});
		}
	}
}

template<typename GraphTileType>
void SSEdgeCentricAlgoKernel<GraphTileType>::
fromScatterTaskvResultToPreMergingTaskv(std::vector<Task>& taskv) const {
    std::vector<Task> new_taskv;
    for (auto& ret : taskv) {
        const auto update = this->getScatterTaskResult(ret);
        const auto dstId = ret.vertexIndex;
        const auto dstTid = ret.dstTid;
        new_taskv.push_back(UpdateType::genTask(update, update, dstId, dstId, dstTid, true, ret.isDummy));
        ret.delete_task_content_buf();
    }
    taskv.clear();
    std::swap(taskv, new_taskv);    
}

template<typename GraphTileType>
void SSEdgeCentricAlgoKernel<GraphTileType>::closeAlgoKernelServer(std::vector<std::thread>& threads) const {
    // TaskComm& clientTaskComm = TaskComm::getClientInstance();
    // clientTaskComm.sendFinish();
    for (auto& thrd : threads)
        thrd.join();
}

}

#endif