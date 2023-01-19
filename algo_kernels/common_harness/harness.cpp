#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <thread>
#include <vector>
#include <condition_variable>
#include <chrono>
#include <ctime>
#include "engine.h"
#include "graph_io_util.h"
#include "graph.h"
#include "graph_common.h"

// Kernel harness header.
#include "kernel_harness.h"

#include "TaskqHandler.h"
#include "TaskUtil.h"

int main(int argc, char* argv[]) {
    /* Parse arguments. */

    size_t threadCount;
    size_t graphTileCount;
    size_t tileIndex; // The tile number is actually the thread number
    uint64_t maxIters;
    uint32_t numParts;
    bool undirected;
    bool noPreprocess;
    bool isCluster;
    bool isRotationBased;

    std::string edgelistFile;
    std::string partitionFile;
    std::string outputFile;
    std::string setting;

    AppArgs appArgs;

    int argRet = algoKernelArgs(argc, argv,
            threadCount, graphTileCount, tileIndex, maxIters, numParts, setting, noPreprocess, isCluster, isRotationBased, undirected,
            edgelistFile, partitionFile, outputFile, appArgs);

    if (argRet) {
        algoKernelArgsPrintHelp(appName, appArgs);
        return argRet;
    }

	// /* Load and initialize the signed enclave */
	// if (open_enclave(tileIndex) != 0) {
	// 	return -1;
	// }

    std::unordered_map<GraphGASLite::VertexIdx, GraphGASLite::TileIdx, std::hash<GraphGASLite::VertexIdx::Type> > tidMap;

    /* Make engine and load input. */

    GraphGASLite::Engine<Graph> engine;
    engine.graphTileIs(GraphGASLite::GraphIOUtil::graphTilesFromEdgeList<Graph>(
                threadCount, tileIndex, edgelistFile, partitionFile, 1, undirected, graphTileCount/threadCount, true, tidMap));
    engine.tileIndexIs(tileIndex);

    std::cout << "Graph loaded from " << edgelistFile <<
        (partitionFile.empty() ? "" : string(" and ") + partitionFile) <<
        " with " << graphTileCount << " graph tiles, " << 
        "into " << threadCount << " tiles." <<
        " Treated as " << (undirected ? "undirected" : "directed") << " graph." <<
        "Current tile is the No." << tileIndex << " tile." <<
        std::endl;

    /* Make algorithm kernel. */

    auto kernel = appArgs.algoKernel<Kernel>(appName);
    kernel->verboseIs(true);
    kernel->maxItersIs(maxIters);
    kernel->numPartsIs(numParts);
    kernel->tidMapIs(tidMap);
    kernel->curTidIs(tileIndex);
    engine.algoKernelNew(kernel);

    CryptoUtil& cryptoUtil = CryptoUtil::getInstance();
    cryptoUtil.tileIndexIs(tileIndex);
    cryptoUtil.setUpPaillierCipher();

    TaskComm& clientTaskComm = TaskComm::getClientInstance();
    clientTaskComm.tileNumIs(threadCount);
    clientTaskComm.tileIndexIs(tileIndex);
    clientTaskComm.settingIs(setting);
    clientTaskComm.noPreprocessIs(noPreprocess);
    clientTaskComm.isClusterIs(isCluster);
    clientTaskComm.isRotationBasedIs(isRotationBased);

    TaskComm& serverTaskComm = TaskComm::getServerInstance();
    serverTaskComm.tileNumIs(threadCount);
    serverTaskComm.tileIndexIs(tileIndex);
    serverTaskComm.settingIs(setting);
    serverTaskComm.noPreprocessIs(noPreprocess);
    serverTaskComm.isClusterIs(isCluster);
    serverTaskComm.isRotationBasedIs(isRotationBased);

    std::thread clientSetupThread([&clientTaskComm](){
        clientTaskComm.setUp(true);
    });

    std::thread serverSetupThread([&serverTaskComm](){
        serverTaskComm.setUp(false);
    });

    clientSetupThread.join();
    serverSetupThread.join();

    std::cout << "Algorithm kernel named " << appName <<
        " is " << algoKernelTagName(kernel->tag()) << ", " <<
        "with max iterations " << maxIters << " and number of partitions " << numParts << "." <<
        std::endl;

    std::cout << "Application parameters: " << appArgs << "." << std::endl;

    /* Run. */

    engine();

    /* Output. */
#ifdef VDATA
    if (!outputFile.empty()) {
        std::cout << "Output to " << outputFile << "." << std::endl;
        std::ofstream ofs(outputFile);
        auto g = engine.graphTile(tileIndex);
        for (auto vIter = g->vertexIter(); vIter != g->vertexIterEnd(); ++vIter) {
            auto v = vIter->second;
            ofs << v->vid() << "\t" << VDATA(v->data()) << std::endl;
        }
        ofs.close();
    }
#endif

    // close_enclave();
    clientTaskComm.closeChannels();

    // std::vector<std::condition_variable>& client_computeTaskq_cv = clientTaskComm.getClientComputeTaskqCv();
	// for (int i = 0; i < threadCount; i++) {
	// 	if (i != tileIndex) {
    //         client_computeTaskq_cv[i].notify_one();
    //     }
    // }

    // for (auto& thrd : client_computeTaskq_threads)
    //     thrd.join();

    serverTaskComm.closeChannels();

    return 0;
}

