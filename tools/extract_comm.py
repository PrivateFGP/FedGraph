# >>> x = [398, 752, 1784, 3520]
# >>> y = [376, 714, 1694, 3320]
# >>> apre = [int((x[i]-y[i])/14) for i in range(len(x))]
# >>> apre = [int((x[i]-y[i])/15) for i in range(len(x))]
# >>> print(apre)
# [1, 2, 6, 13]
# >>> apre = [(x[i]-y[i])/15 for i in range(len(x))]
# >>> print(apre)
# [1.4666666666666666, 2.533333333333333, 6.0, 13.333333333333334]
# >>> online = [(y[i])/2 for i in range(len(x))]
# >>> print(online)
# [188.0, 357.0, 847.0, 1660.0]

import numpy as np
from statistics import mean

commTags = ["clientPreprocessComm", "serverPreprocessComm"]
commHeader = ["Pre-merging Comm"]

log_root_path = "./log/"
comm_root_path = "./comm/"
# executable = "sssp-ss"

doPreprocess = False
defaultBandWidth = 400 #Mbps
defaultLatency = 1 #ms
defaultInterRatio = 0.4
defaultNumParts=2
defaultScaler=0.2

plainNumList = [2,1]

output_matrix = []

bandwidthList = [200, 400, 1000, 4000]

batchSize = 5

def add_table_split():
    table_tag = []
    for i in range(len(output_matrix[0])):
        table_tag.append("----")
    output_matrix.append(table_tag)

def print_matrix():
    for line in output_matrix:
        print("| " + " | ".join(line) + " |")

def getCommRow(fileName, isOnline=False):
    with open(fileName, "r", encoding='utf-8') as ifile:
        # print(fileName)
        commDict = {}
        commList = []
        for tag in commTags:
            commDict[tag] = []
        lines = ifile.readlines()
        
        for line in lines:
            line = line.strip("\n")
            line = line.strip(b'\xef\xbb\xbf'.decode('utf-8'))
            if line.startswith(">>"):
                for tag in commTags:
                    if line.startswith(">>" + tag + " "):
                        # print(line) 
                        commDict[tag].append(float(line.split(" ")[2]))
                        break
        
        for tag in commTags:
            commDict[tag] = float(mean(commDict[tag]))
        
        commDict["pre-merging"] = commDict["clientPreprocessComm"] + commDict["serverPreprocessComm"]

        # commList.append(str(commDict["pre-merging"]))
        cur_comm = commDict["pre-merging"]
        commList.append(f"{cur_comm:.2f}")

        return commList

def getTotalComm(fileName):
    with open(fileName, "r", encoding='utf-8') as ifile:
        lines = ifile.readlines()
        for line in lines:
            line = line.strip("\n")
            line = line.strip(b'\xef\xbb\xbf'.decode('utf-8'))
            if "vethA.peer" in line:
                segs = line.split()
                segs = [x for x in segs if x]
                cur_comm = segs[2]
                cur_comm = cur_comm.strip("MB")
                cur_comm = float(cur_comm)
                return cur_comm

def printTotalComm(executable_name):
    scalerList = [0.05, 0.2, 0.5, 1]
    noPreprocessList = []
    hasPreprocessList = []
    for scaler in scalerList:
        noPreprocessList.append(getTotalComm(comm_root_path + str("Falsepreprocess_") + str(scaler) + "scaler_" + executable_name + ".comm"))
        hasPreprocessList.append(getTotalComm(comm_root_path + str("Truepreprocess_") + str(scaler) + "scaler_" + executable_name + ".comm"))
    preprocessList = [(hasPreprocessList[i]-noPreprocessList[i])/5 for i in range(len(noPreprocessList))]
    print([f"{x:.2f}" for x in noPreprocessList])
    print([f"{x:.2f}" for x in preprocessList])
        
def print_round_halved_comm_table(executable_name):
    global output_matrix
    output_matrix = [[executable_name + " comm"] + commHeader]
    add_table_split()
    root_setting = "graph_scale"
    scalerList = [0.05, 0.2, 0.5, 1]
    log_path = log_root_path + "graph_scale/"
    for scaler in scalerList:
        i = 0
        row = []
        eval_setting = root_setting+"_"+executable_name+"_"+str(scaler)+"scale_"+str(defaultNumParts)+"p_"+str(defaultInterRatio)+"inter"+"_"+str(defaultBandWidth)+"b_"+str(defaultLatency)+"l"+"_online"
        log_f_name = log_path+eval_setting+"_"+str(i)+".log"
        durationRow = getCommRow(log_f_name, isOnline=True)
        row.append(durationRow[0])
        row = [str(scaler)] + row
        output_matrix.append(row)
    print_matrix()

def print_rotation_based_comm_table(executable_name):
    global output_matrix
    output_matrix = [[executable_name + " comm"] + commHeader]
    add_table_split()
    root_setting = "rotation_based"
    scalerList = [0.05, 0.2, 0.5, 1]
    log_path = log_root_path + "rotation_based/"
    for scaler in scalerList:
        i = 0
        row = []
        eval_setting = root_setting+"_"+"rotation_based_"+executable_name+"_"+str(scaler)+"scale_"+str(defaultNumParts)+"p_"+str(defaultInterRatio)+"inter"+"_"+str(defaultBandWidth)+"b_"+str(defaultLatency)+"l"+"_online"
        log_f_name = log_path+eval_setting+"_"+str(i)+".log"
        durationRow = getCommRow(log_f_name, isOnline=True)
        row.append(durationRow[0])
        row = [str(scaler)] + row
        output_matrix.append(row)
    print_matrix()


# Scale Table
# executableList = ["bfs-ss", "sssp-ss", "cc-ss", "pagerank-ss"]
executableList = ["bfs-ss", "sssp-ss", "cc-ss", "pagerank-ss"]
defaultExecutable = "sssp-ss"

print(">>>> Round-Halved")
print_round_halved_comm_table(defaultExecutable)

print(">>>> Shift-Based")
print_rotation_based_comm_table(defaultExecutable)

print(">>>> Total Comm")
printTotalComm(defaultExecutable)

