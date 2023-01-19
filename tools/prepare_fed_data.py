#! /usr/bin/python

import os
import numpy as np
import argparse
import random

parser = argparse.ArgumentParser(description=\
        'Partition edge-list format graph with destination-mirror-only.')

# parser.add_argument('partitionFileName', help='Output partition file')
parser.add_argument('edgelistFileName', help='Input edge list file')

args = parser.parse_args()

# partitionFileName = os.path.expanduser(args.partitionFileName)
edgelistFileName = os.path.expanduser(args.edgelistFileName)


'''
Return a dict with mapping from vertex index to partition index.
'''
def extend_partition(inputEdgelistFileName, numParts, interEdgeRatio, scaler=0.2, mirrorDstDegree = 8):

    class VertexInfo:
        def __init__(self):
            self.dvSet = set()

        def addDstVertex(self, vid):
            self.dvSet.add(vid)

    # Vertex partition dict.
    vpartDict = {}
    # Vertex info dict.
    vinfoDict = {}
    # Edge count
    ne = 0

    outputEdgeList = []
    inputEdgeList = []
    inputVertexList = []
    maxVertexIndex = 0

    with open(inputEdgelistFileName, 'r') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue

            elems = line.split()
            assert len(elems) >= 2
            src = int(elems[0])
            dst = int(elems[1])

            # inputEdgeList.append([src, dst])

            if src > maxVertexIndex:
                maxVertexIndex = src
            if dst > maxVertexIndex:
                maxVertexIndex = dst

            vinfoDict.setdefault(src, VertexInfo()).addDstVertex(dst)
            vinfoDict.setdefault(dst, VertexInfo())
            vpartDict.setdefault(src, -1)
            vpartDict.setdefault(dst, -1)
            ne += 1

            if ne % 1000000 == 0:
                print('[partition] Read {} edges'.format(ne))


    print('[partition] {} edges, {} vertices'.format(ne, len(vinfoDict)))

    partitionGap = maxVertexIndex + 1

    if scaler != 1:
        if scaler < 1: # Shrink
            cnt = 0
            step = int(1/float(scaler))
            for (v, info) in dict(vinfoDict).items():
                if cnt%step != 0:
                    del vinfoDict[v]
                cnt+=1
        else: # Extend
            extendDict = {}
            for (v, info) in vinfoDict.items():
                for i in range(1, scaler):
                    vi = VertexInfo()
                    vi.dvSet = set([x+partitionGap*i for x in info.dvSet])
                    extendDict.setdefault(v+partitionGap*i, vi)
            vinfoDict.update(extendDict)
            partitionGap = scaler * partitionGap       

    for (v, info) in vinfoDict.items():
        inputVertexList.append(v)
        for dst in info.dvSet:
            inputEdgeList.append([v, dst])
            inputVertexList.append(dst)
    
    inputVertexList = list(set(inputVertexList))

    for i in range(numParts):
        offset = i * partitionGap
        # Duplicate edges by numParts
        for edge in inputEdgeList:
            outputEdgeList.append([v + offset for v in edge])
        # Assign partition index to vertices
        for v in inputVertexList:
            vpartDict[v + offset] = i

    for (v, part) in dict(vpartDict).items():
        if part == -1:
            del vpartDict[v]
    
    totalEdgeNum = len(outputEdgeList)
    interEdgeNum = int(totalEdgeNum * interEdgeRatio)
    onePartySendOutInterEdgeNum = int(interEdgeNum / (numParts * (numParts - 1)))
    onePartyDstVertexNum = int(onePartySendOutInterEdgeNum / mirrorDstDegree)
    print("onePartyDstVertexNum", onePartyDstVertexNum)
    print("totalEdgeNum", totalEdgeNum)
    # Add inter-edges
    for i in range(numParts):
        for j in range(numParts):
            if (i != j):
                for k in range(onePartyDstVertexNum):
                    for m in range(mirrorDstDegree):
                        outputEdgeList.append([inputVertexList[m]+i*partitionGap, inputVertexList[k]+j*partitionGap])

    return outputEdgeList, vpartDict

def writeHead(fh, numParts):
    fh.write('# Partitioned into {} partitions.\n'.format(numParts))
    fh.write('\n')

def writePartitionToFile(fileName, partDict, numParts = 2):
    with open(fileName, 'w') as fh:
        writeHead(fh, numParts)
        for (v, p) in partDict.items():
            fh.write('{}\t{}\n'.format(v, p))

def writeEdgeListToFile(fileName, edgeList):
    with open(fileName, 'w') as fh:
        for e in edgeList:
            fh.write('{} {}\n'.format(e[0], e[1]))

def my_makedir(path):
    if not os.path.isdir(path):
        os.makedirs(path)

defaultInterRatio = 0.4
defaultNumParts = 2

# Party num setting
numPartsList = [2, 3, 4, 5]
outputPathRoot = "data/party_num/"
my_makedir(outputPathRoot)
for num in numPartsList:
    outputEdgeList, vpartDict = extend_partition(edgelistFileName, num, defaultInterRatio)
    writeEdgeListToFile(outputPathRoot+str(num)+"p.dat", outputEdgeList)
    print("party_num " + str(num) + ":", "|E|=" + str(len(outputEdgeList)))
    writePartitionToFile(outputPathRoot+str(num)+"p.part", vpartDict, num)
    print("party_num " + str(num) + ":", "|V|=" + str(len(vpartDict)))

# Inter-edge ratio setting
interRatioList = [0.1, 0.2, 0.4, 0.8, 1]
outputPathRoot = "data/inter_ratio/"
my_makedir(outputPathRoot)
for ratio in interRatioList:
    outputEdgeList, vpartDict = extend_partition(edgelistFileName, defaultNumParts, ratio)
    writeEdgeListToFile(outputPathRoot+str(ratio)+"inter.dat", outputEdgeList)
    print("inter_ratio " + str(ratio) + ":", "|E|=" + str(len(outputEdgeList)))
    writePartitionToFile(outputPathRoot+str(ratio)+"inter.part", vpartDict)
    print("inter_ratio " + str(ratio) + ":", "|V|=" + str(len(vpartDict)))

# Graph scale
# scalerList = [0.2, 0.5, 1, 2, 5]
scalerList = [0.05, 0.2, 0.5, 1]
outputPathRoot = "data/graph_scale/"
my_makedir(outputPathRoot)
for scaler in scalerList:
    outputEdgeList, vpartDict = extend_partition(edgelistFileName, defaultNumParts, defaultInterRatio, scaler)
    writeEdgeListToFile(outputPathRoot+str(scaler)+"scale.dat", outputEdgeList)
    print("scaler " + str(scaler) + ":", "|E|=" + str(len(outputEdgeList)))
    writePartitionToFile(outputPathRoot+str(scaler)+"scale.part", vpartDict)
    print("scaler " + str(scaler) + ":", "|V|=" + str(len(vpartDict)))