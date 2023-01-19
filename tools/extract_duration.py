import numpy as np
from statistics import mean

durationTags = ["preprocess", "iteration", "preprocess_OM", "Scatter" ,"Scatter_preparation", 
                    "Scatter_computation", "Gather", "premerging", "premerged_extraction", "Gather_preparation", "Gather_computation"]
breakDownTableHeader = ["Preprocess", "Online Iteration", "OM Preprocess", "Scatter" ,"Scatter Preparation", 
                    "Scatter Computation", "Gather", "Pre-merging", "Pre-merged Extraction", "Gather Preparation", "Gather Computation"]

log_root_path = "./log/"
# executable = "sssp-ss"

doPreprocess = False
defaultBandWidth = 400 #Mbps
defaultLatency = 1 #ms
defaultInterRatio = 0.4
defaultNumParts=2
defaultScaler=0.2

# plainNumList = [2,2,2,1]
plainNumList = [2,1]
keySize = 3072
# batchSize = [keySize/(x*65) for x in plainNumList]
batchSize = [5,5]


output_matrix = []

# bandwidthList = [50, 200, 400, 1000]
bandwidthList = [200, 400, 1000, 4000]

def add_table_split():
    table_tag = []
    for i in range(len(output_matrix[0])):
        table_tag.append("----")
    output_matrix.append(table_tag)

def print_matrix():
    for line in output_matrix:
        print("| " + " | ".join(line) + " |")

def getDurationRow(fileName, isOnline=False):
    # print(fileName)
    with open(fileName, "r", encoding='utf-8') as ifile:
        durationDict = {}
        durationList = []
        for tag in durationTags:
            durationDict[tag] = []
        lines = ifile.readlines()
        
        for line in lines:
            line = line.strip("\n")
            line = line.strip(b'\xef\xbb\xbf'.decode('utf-8'))
            if line.startswith("::"):
                # print(line)
                for tag in durationTags:
                    if line.startswith("::" + tag + " "): 
                        durationDict[tag].append(float(line.split(" ")[2]))
                        break
        
        for tag in durationTags:
            # print(tag, durationDict[tag])
            if tag.startswith("preprocess") and isOnline:
                durationDict[tag] = 0
                continue
            if tag in ["Scatter_computation", "premerging", "premerged_extraction"]:
                durationDict[tag] = float(mean(durationDict[tag]) * 2)
            elif tag not in ["Scatter", "Gather"]:
                durationDict[tag] = float(mean(durationDict[tag]))
        
        durationDict["Scatter"] = durationDict["Scatter_preparation"] + durationDict["Scatter_computation"]
        durationDict["Gather"] = durationDict["premerging"] + durationDict["premerged_extraction"] + durationDict["Gather_preparation"] + durationDict["Gather_computation"]

        # Fix me
        durationDict["preprocess"] = float(durationDict["preprocess"] / batchSize[-1]) + 1
        durationDict["preprocess_OM"] = float(durationDict["preprocess_OM"] / batchSize[-1]) + 1

        for tag in durationTags:
            # durationList.append(str(durationDict[tag]))
            durationList.append(f"{durationDict[tag]:.2f}")

        return durationList

# Graph scale
# scalerList = [0.2, 0.5, 1, 2, 5]
def print_scale_table(executable_name, cur_bandwidth, cur_latency):
    global output_matrix
    output_matrix = [[executable_name + " graph scaler"] + breakDownTableHeader]
    add_table_split()
    root_setting = "graph_scale"
    # scalerList = [1, 2, 5, 10]
    scalerList = [0.05, 0.2, 0.5, 1]
    log_path = log_root_path + "graph_scale/"
    for scaler in scalerList:
        i = 0
        eval_setting = root_setting+"_"+executable_name+"_"+str(scaler)+"scale_"+str(defaultNumParts)+"p_"+str(defaultInterRatio)+"inter"+"_"+str(cur_bandwidth)+"b_"+str(cur_latency)+"l"
        log_f_name = log_path+eval_setting+"_"+str(i)+".log"
        row = getDurationRow(log_f_name)
        row = [str(scaler)] + row
        output_matrix.append(row)
    print_matrix()

def print_scale_bandwidth_table(executable_name):
    global output_matrix
    output_matrix = [[executable_name + " scale bandwidth"] + [str(bandwidth) for bandwidth in bandwidthList]]
    add_table_split()
    root_setting = "graph_scale"
    scalerList = [0.05, 0.2, 0.5, 1]
    log_path = log_root_path + "graph_scale/"
    for scaler in scalerList:
        i = 0
        row = []
        for bandwidth in bandwidthList:
            eval_setting = root_setting+"_"+executable_name+"_"+str(scaler)+"scale_"+str(defaultNumParts)+"p_"+str(defaultInterRatio)+"inter"+"_"+str(bandwidth)+"b_"+str(defaultLatency)+"l"+"_online"
            log_f_name = log_path+eval_setting+"_"+str(i)+".log"
            durationRow = getDurationRow(log_f_name, isOnline=True)
            row.append(durationRow[durationTags.index("iteration")])
        row = [str(scaler)] + row
        output_matrix.append(row)
    print_matrix()

# Inter-edge ratio setting
# interRatioList = [0.05, 0.1, 0.2, 0.4, 0.8]
def print_inter_ratio_bandwidth(executable_name):
    global output_matrix
    output_matrix = [[executable_name + " inter-ratio bandwidth"] + [str(bandwidth) for bandwidth in bandwidthList]]
    add_table_split()
    root_setting = "inter_ratio"
    interRatioList = [0.1, 0.2, 0.4, 0.8, 1]
    log_path = log_root_path + "inter_ratio/"
    for ratio in interRatioList:
        i = 0
        row = []
        for bandwidth in bandwidthList:
            eval_setting = root_setting+"_"+executable_name+"_"+str(defaultScaler)+"scale_"+str(defaultNumParts)+"p_"+str(ratio)+"inter"+"_"+str(bandwidth)+"b_"+str(defaultLatency)+"l"+"_online"
            log_f_name = log_path+eval_setting+"_"+str(i)+".log"
            durationRow = getDurationRow(log_f_name, isOnline=True)
            row.append(durationRow[durationTags.index("iteration")])
        row = [str(ratio)] + row
        output_matrix.append(row)
    print_matrix()

# Party num setting
def print_party_num_bandwidth(executable_name):
    global output_matrix
    output_matrix = [[executable_name + " party num bandwidth"] + [str(bandwidth) for bandwidth in bandwidthList]]
    add_table_split()
    root_setting = "party_num"
    numPartsList = [2,3,4,5]
    log_path = log_root_path + "party_num/"
    for num in numPartsList:
        i = 0
        row = []
        for bandwidth in bandwidthList:
            eval_setting = root_setting+"_"+executable_name+"_"+str(defaultScaler)+"scale_"+str(num)+"p_"+str(defaultInterRatio)+"inter"+"_"+str(bandwidth)+"b_"+str(defaultLatency)+"l"+"_online"
            log_f_name = log_path+eval_setting+"_"+str(i)+".log"
            # print(log_f_name)
            durationRow = getDurationRow(log_f_name, isOnline=True)
            row.append(durationRow[durationTags.index("iteration")])
        row = [str(num)] + row
        output_matrix.append(row)
    print_matrix()

def print_latency_bandwidth(executable_name):
    global output_matrix
    output_matrix = [[executable_name + " latency bandwidth"] + [str(bandwidth) for bandwidth in bandwidthList]]
    add_table_split()
    root_setting = "latency"
    latencyList = [0.15, 1, 10, 20]
    log_path = log_root_path + "latency/"
    for latency in latencyList:
        i = 0
        row = []
        for bandwidth in bandwidthList:
            eval_setting = root_setting+"_"+executable_name+"_"+str(defaultScaler)+"scale_"+str(defaultNumParts)+"p_"+str(defaultInterRatio)+"inter"+"_"+str(bandwidth)+"b_"+str(latency)+"l"+"_online"
            log_f_name = log_path+eval_setting+"_"+str(i)+".log"
            # print(log_f_name)
            durationRow = getDurationRow(log_f_name, isOnline=True)
            row.append(durationRow[durationTags.index("iteration")])
        row = [str(latency)] + row
        output_matrix.append(row)
    print_matrix()

def print_round_halved_duration_table(executable_name):
    global output_matrix
    output_matrix = [[executable_name + " duration"] + ["pre-merging duration", "online duration"]]
    add_table_split()
    root_setting = "graph_scale"
    scalerList = [0.05, 0.2, 0.5, 1]
    log_path = log_root_path + "graph_scale/"
    for scaler in scalerList:
        i = 0
        row = []
        eval_setting = root_setting+"_"+executable_name+"_"+str(scaler)+"scale_"+str(defaultNumParts)+"p_"+str(defaultInterRatio)+"inter"+"_"+str(defaultBandWidth)+"b_"+str(defaultLatency)+"l"+"_online"
        log_f_name = log_path+eval_setting+"_"+str(i)+".log"
        durationRow = getDurationRow(log_f_name, isOnline=True)
        row.append(durationRow[durationTags.index("premerging")])
        row.append(durationRow[durationTags.index("iteration")])
        row = [str(scaler)] + row
        output_matrix.append(row)
    print_matrix()

def print_rotation_based_duration_table(executable_name):
    global output_matrix
    output_matrix = [[executable_name + " duration"] + ["pre-merging duration", "online duration"]]
    add_table_split()
    root_setting = "rotation_based"
    scalerList = [0.05, 0.2, 0.5, 1]
    log_path = log_root_path + "rotation_based/"
    for scaler in scalerList:
        i = 0
        row = []
        eval_setting = root_setting+"_"+"rotation_based_"+executable_name+"_"+str(scaler)+"scale_"+str(defaultNumParts)+"p_"+str(defaultInterRatio)+"inter"+"_"+str(defaultBandWidth)+"b_"+str(defaultLatency)+"l"+"_online"
        log_f_name = log_path+eval_setting+"_"+str(i)+".log"
        durationRow = getDurationRow(log_f_name, isOnline=True)
        row.append(durationRow[durationTags.index("premerging")])
        row.append(durationRow[durationTags.index("iteration")])
        row = [str(scaler)] + row
        output_matrix.append(row)
    print_matrix()

# Scale Table
# executableList = ["bfs-ss", "sssp-ss", "cc-ss", "pagerank-ss"]
executableList = ["sssp-ss", "pagerank-ss"]
defaultExecutable = "sssp-ss"

print(">>>> Scale Table")
for executable in executableList:
    # print_party_num(executable)
    # print_inter_edge_ratio(executable)
    print_scale_table(executable, defaultBandWidth, defaultLatency)

print(">>>> Scale Bandwidth Figure")
# Executbale Scale Bandwidth
for executable in executableList:
    print_scale_bandwidth_table(executable)

print(">>>> Inter-Ratio Bandwidth Figure")
# Executbale inter-ratio Bandwidth
for executable in executableList:
    print_inter_ratio_bandwidth(executable)

print(">>>> Party Num Figure")
# Executable party num bandwidth
for executable in executableList:
    print_party_num_bandwidth(executable)

# print(">>>> Latency Figure")
# # Executable party num bandwidth
# for executable in executableList:
#     print_latency_bandwidth(executable)

print(">>>> Round-Halved Duration Table")
print_round_halved_duration_table(defaultExecutable)

print(">>>> Shift-Based Duration Table")
print_rotation_based_duration_table(defaultExecutable)

