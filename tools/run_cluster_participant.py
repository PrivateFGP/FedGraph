import subprocess
import sys
import os
from threading import Thread
import psutil
import time
import pandas as pd

executable_root_path = "./../../bin/"
data_root_path = "./data/"
result_root_path = "./result/"
log_root_path = "./log/"
preprocess_root_path = "./preprocess/"
comm_root_path = "./comm/"
iterations = 1
defaultNumParts = 2
doPreprocess = False
defaultBandWidth = 400 #Mbps
defaultLatency = 1 #ms
isCluster = True
defaultInterRatio = 0.4
defaultScaler = 0.2
bandwidthList = [200, 400, 1000, 4000]
latencyList = [0.15, 1, 10, 20]

def my_makedir(path):
    if not os.path.isdir(path):
        os.makedirs(path)

# print(process.stdout.read().decode('utf-8'))

# python batch_Deckard.py > test_case/batch_similarity.txt

# ./sssp-ss -t 5 -g 5 -i 4 -m 2 -p 1 ./../data/test.dat ./../data/test5p.part ./../data/result_4.txt 0

processList = []

def setup_network(bandwidth, latency):
    # Clean first
    curProc = subprocess.Popen(["sudo", "bash", "./scripts/clean_network.sh"])
    curProc.wait()
    # Setup network
    curProc = subprocess.Popen(["sudo", "bash", "./scripts/setup_network.sh", str(bandwidth), str(latency)])
    curProc.wait()

def clean_network():
    curProc = subprocess.Popen(["sudo", "bash", "./scripts/clean_network.sh"])
    curProc.wait()

nsList = ["A", "B", "C", "D", "E"]

def get_size(bytes):
    """
    Returns size of bytes in a nice format
    """
    # for unit in ['', 'K', 'M', 'G', 'T', 'P']:
    #     if bytes < 1024:
    #         return f"{bytes:.2f}{unit}B"
    #     bytes /= 1024
    bytes /= 1024
    bytes /= 1024
    return f"{bytes:.2f}MB"

def commStart():
    io = psutil.net_io_counters(pernic=True)
    return io

def commEnd(io, doPreprocess, scaler, tag=""):
    # get the network I/O stats again per interface 
    io_2 = psutil.net_io_counters(pernic=True)
    # initialize the data to gather (a list of dicts)
    data = []
    for iface, iface_io in io.items():
        # new - old stats gets us the speed
        upload, download = io_2[iface].bytes_sent - iface_io.bytes_sent, io_2[iface].bytes_recv - iface_io.bytes_recv
        data.append({
            "iface": iface, "Download": get_size(download),
            "Upload": get_size(upload),
        })
    # construct a Pandas DataFrame to print stats in a cool tabular style
    df = pd.DataFrame(data)
    # sort values per column, feel free to change the column
    df.sort_values("Download", inplace=True, ascending=False)
    # clear the screen based on your OS
    # os.system("cls") if "nt" in os.name else os.system("clear")
    # print the stats
    # print(df.to_string())
    my_makedir(comm_root_path)
    with open(comm_root_path + str(doPreprocess) + "preprocess_" + str(scaler) + "scaler" + tag + ".comm", "w", encoding='utf-8') as ofile:
        ofile.write(df.to_string())


# Party num setting
def run_party_num(executable, executable_path, cur_bandwidth, cur_latency):
    root_setting = "party_num"
    numPartsList = [2,3,4,5]
    data_path = data_root_path + "party_num/"
    result_path = result_root_path + "party_num/"
    log_path = log_root_path + "party_num/"
    my_makedir(result_path)
    my_makedir(log_path)
    for num in numPartsList:
        preprocess_setting = root_setting+"_"+executable+"_"+str(defaultScaler)+"scale_"+str(num)+"p_"+str(defaultInterRatio)+"inter"
        eval_setting = preprocess_setting+"_"+str(cur_bandwidth)+"b_"+str(cur_latency)+"l"
        if not doPreprocess:
            eval_setting += "_online"
        preprocess_path = preprocess_root_path + preprocess_setting + "/"
        my_makedir(preprocess_path)
        processList = []
        for i in range(num):
            cmd = []
            if isCluster:
                cmd += ["sudo", "ip", "netns", "exec", nsList[i]]
            cmd += [executable_path, "-t", str(num), "-g", str(num), "-i", str(i), "-m", str(iterations), "-p", "1", "-s", preprocess_setting]
            if not doPreprocess:
                cmd += ["-n", "1"]
            if isCluster:
                cmd += ["-c", "1"]
            cmd += [data_path+str(num)+"p.dat", data_path+str(num)+"p.part",
                    result_path+eval_setting+"_"+str(i)+".result"]
            print(" ".join(cmd))
            log_f = open(log_path+eval_setting+"_"+str(i)+".log", 'w', encoding='utf-8')
            processList.append(subprocess.Popen(cmd, stdout=log_f))
        for process in processList:
            process.wait()
        processList = []

# Inter-edge ratio setting
# interRatioList = [0.1, 0.2, 0.4, 0.8, 1]
def run_inter_ratio(executable, executable_path, cur_bandwidth, cur_latency):
    root_setting = "inter_ratio"
    interRatioList = [0.1, 0.2, 0.4, 0.8, 1]
    data_path = data_root_path + "inter_ratio/"
    result_path = result_root_path + "inter_ratio/"
    log_path = log_root_path + "inter_ratio/"
    my_makedir(result_path)
    my_makedir(log_path)
    for ratio in interRatioList:
        preprocess_setting = root_setting+"_"+executable+"_"+str(defaultScaler)+"scale_"+str(defaultNumParts)+"p_"+str(ratio)+"inter"
        eval_setting = preprocess_setting+"_"+str(cur_bandwidth)+"b_"+str(cur_latency)+"l"
        if not doPreprocess:
            eval_setting += "_online"
        preprocess_path = preprocess_root_path + preprocess_setting + "/"
        my_makedir(preprocess_path)
        processList = []
        for i in range(defaultNumParts):
            cmd = []
            if isCluster:
                cmd = ["sudo", "ip", "netns", "exec", nsList[i]]
            cmd += [executable_path, "-t", str(defaultNumParts), "-g", str(defaultNumParts), "-i", str(i), "-m", str(iterations), "-p", "1", "-s", preprocess_setting]
            if not doPreprocess:
                cmd += ["-n", "1"]
            if isCluster:
                cmd += ["-c", "1"]
            cmd += [data_path+str(ratio)+"inter.dat", data_path+str(ratio)+"inter.part",
                    result_path+eval_setting+"_"+str(i)+".result"]
            print(" ".join(cmd))
            log_f = open(log_path+eval_setting+"_"+str(i)+".log", 'w', encoding='utf-8')
            processList.append(subprocess.Popen(cmd, stdout=log_f))
        for process in processList:
            process.wait()
        processList = []

# Graph scale
# scalerList = [0.1, 0.2, 0.4, 0.8, 1]
def run_graph_scale(executable, executable_path, cur_bandwidth, cur_latency):
    # Evaluation root setting
    root_setting = "graph_scale"
    scalerList = [0.05, 0.2, 0.5, 1]
    # scalerList = [0.1]
    data_path = data_root_path + "graph_scale/"
    result_path = result_root_path + "graph_scale/"
    log_path = log_root_path + "graph_scale/"
    my_makedir(result_path)
    my_makedir(log_path)
    for scaler in scalerList:
        preprocess_setting = root_setting+"_"+executable+"_"+str(scaler)+"scale_"+str(defaultNumParts)+"p_"+str(defaultInterRatio)+"inter"
        eval_setting = preprocess_setting+"_"+str(cur_bandwidth)+"b_"+str(cur_latency)+"l"
        if not doPreprocess:
            eval_setting += "_online"
        preprocess_path = preprocess_root_path + preprocess_setting + "/"
        my_makedir(preprocess_path)
        processList = []
        start_io = commStart()
        for i in range(defaultNumParts):
            cmd = []
            if isCluster:
                cmd = ["sudo", "ip", "netns", "exec", nsList[i]]
            cmd += [executable_path, "-t", str(defaultNumParts), "-g", str(defaultNumParts), "-i", str(i), "-m", str(iterations), "-p", "1", "-s", preprocess_setting]
            if not doPreprocess:
                cmd += ["-n", "1"]
            if isCluster:
                cmd += ["-c", "1"]
            cmd += [data_path+str(scaler)+"scale.dat", data_path+str(scaler)+"scale.part",
                    result_path+eval_setting+"_"+str(i)+".result"]
            print(" ".join(cmd))
            log_f = open(log_path+eval_setting+"_"+str(i)+".log", 'w', encoding='utf-8')
            processList.append(subprocess.Popen(cmd, stdout=log_f))
        for process in processList:
            process.wait()
        commEnd(start_io, doPreprocess, scaler, str("_") + executable)
        processList = []

# Graph scale
# scalerList = [0.1, 0.2, 0.4, 0.8, 1]
def run_rotation_based(executable, executable_path, cur_bandwidth, cur_latency):
    # Evaluation root setting
    root_setting = "rotation_based"
    scalerList = [0.05, 0.2, 0.5, 1]
    data_path = data_root_path + "graph_scale/"
    result_path = result_root_path + "rotation_based/"
    log_path = log_root_path + "rotation_based/"
    my_makedir(result_path)
    my_makedir(log_path)
    for scaler in scalerList:
        preprocess_setting = "rotation_based"+"_"+executable+"_"+str(scaler)+"scale_"+str(defaultNumParts)+"p_"+str(defaultInterRatio)+"inter"
        eval_setting = root_setting + "_" + preprocess_setting+"_"+str(cur_bandwidth)+"b_"+str(cur_latency)+"l"
        if not doPreprocess:
            eval_setting += "_online"
        preprocess_path = preprocess_root_path + preprocess_setting + "/"
        my_makedir(preprocess_path)
        processList = []
        start_io = commStart()
        for i in range(defaultNumParts):
            cmd = []
            if isCluster:
                cmd = ["sudo", "ip", "netns", "exec", nsList[i]]
            cmd += [executable_path, "-t", str(defaultNumParts), "-g", str(defaultNumParts), "-i", str(i), "-m", str(iterations), "-p", "1", "-s", preprocess_setting]
            if not doPreprocess:
                cmd += ["-n", "1"]
            if isCluster:
                cmd += ["-c", "1"]
            cmd += ["-r", "1"]
            cmd += [data_path+str(scaler)+"scale.dat", data_path+str(scaler)+"scale.part",
                    result_path+eval_setting+"_"+str(i)+".result"]
            print(" ".join(cmd))
            log_f = open(log_path+eval_setting+"_"+str(i)+".log", 'w', encoding='utf-8')
            processList.append(subprocess.Popen(cmd, stdout=log_f))
        for process in processList:
            process.wait()
        commEnd(start_io, doPreprocess, scaler, str("_") + executable + "_rotation_based")
        processList = []

# Latency
def run_latency(executable, executable_path, cur_bandwidth, cur_latency):
    # Evaluation root setting
    root_setting = "latency"
    data_path = data_root_path + "graph_scale/"
    result_path = result_root_path + "latency/"
    log_path = log_root_path + "latency/"
    my_makedir(result_path)
    my_makedir(log_path)
        
    preprocess_setting = root_setting+"_"+executable+"_"+str(defaultScaler)+"scale_"+str(defaultNumParts)+"p_"+str(defaultInterRatio)+"inter"
    eval_setting = preprocess_setting+"_"+str(cur_bandwidth)+"b_"+str(cur_latency)+"l"
    if not doPreprocess:
        eval_setting += "_online"
    preprocess_path = preprocess_root_path + preprocess_setting + "/"
    my_makedir(preprocess_path)
    processList = []
    for i in range(defaultNumParts):
        cmd = []
        if isCluster:
            cmd = ["sudo", "ip", "netns", "exec", nsList[i]]
        cmd += [executable_path, "-t", str(defaultNumParts), "-g", str(defaultNumParts), "-i", str(i), "-m", str(iterations), "-p", "1", "-s", preprocess_setting]
        if not doPreprocess:
            cmd += ["-n", "1"]
        if isCluster:
            cmd += ["-c", "1"]
        cmd += [data_path+str(defaultScaler)+"scale.dat", data_path+str(defaultScaler)+"scale.part",
                result_path+eval_setting+"_"+str(i)+".result"]
        print(" ".join(cmd))
        log_f = open(log_path+eval_setting+"_"+str(i)+".log", 'w', encoding='utf-8')
        processList.append(subprocess.Popen(cmd, stdout=log_f))
    for process in processList:
        process.wait()
    processList = []

# executableList = ["bfs-ss", "sssp-ss", "cc-ss", "pagerank-ss"]
# executableList = ["sssp-ss"]
defaultExecutable = "sssp-ss"
executableList = ["sssp-ss", "pagerank-ss"]

# print(">>>> Breakdown Table")
# # Table 1
# doPreprocess = True
# cur_bandwidth = defaultBandWidth
# cur_latency = defaultLatency
# if isCluster:
#     setup_network(cur_bandwidth, cur_latency)
#     # exit()
# for executable in executableList:
#     executable_path = executable_root_path + executable
#     run_graph_scale(executable, executable_path, cur_bandwidth, cur_latency)
# if isCluster:
#     clean_network()

# print(">>>> scale-bandwidth")
# # Figure 1 and Table 1
# doPreprocess = False
# cur_latency = defaultLatency
# for cur_bandwidth in bandwidthList:
#     if isCluster:
#         setup_network(cur_bandwidth, cur_latency)
#     for executable in executableList:
#         executable_path = executable_root_path + executable
#         run_graph_scale(executable, executable_path, cur_bandwidth, cur_latency)
#     if isCluster:
#         clean_network()

# print(">>>> inter-ratio-bandwidth")
# # Figure 2
# doPreprocess = True
# cur_bandwidth = defaultBandWidth
# cur_latency = defaultLatency
# if isCluster:
#     setup_network(cur_bandwidth, cur_latency)
# for executable in executableList:
#     executable_path = executable_root_path + executable
#     run_inter_ratio(executable, executable_path, cur_bandwidth, cur_latency)
# if isCluster:
#     clean_network()

# doPreprocess = False
# cur_latency = defaultLatency
# for cur_bandwidth in bandwidthList:
#     if isCluster:
#         setup_network(cur_bandwidth, cur_latency)
#     for executable in executableList:
#         executable_path = executable_root_path + executable
#         run_inter_ratio(executable, executable_path, cur_bandwidth, cur_latency)
#     if isCluster:
#         clean_network()

# print(">>>> party-num-bandwidth")
# # Figure 3
# doPreprocess = True
# cur_bandwidth = defaultBandWidth
# cur_latency = defaultLatency
# if isCluster:
#     setup_network(cur_bandwidth, cur_latency)
# for executable in executableList:
#     executable_path = executable_root_path + executable
#     run_party_num(executable, executable_path, cur_bandwidth, cur_latency)
# if isCluster:
#     clean_network()

# doPreprocess = False
# cur_latency = defaultLatency
# for cur_bandwidth in bandwidthList:
#     if isCluster:
#         setup_network(cur_bandwidth, cur_latency)
#     for executable in executableList:
#         executable_path = executable_root_path + executable
#         run_party_num(executable, executable_path, cur_bandwidth, cur_latency)
#     if isCluster:
#         clean_network()

# # print(">>>> latency-bandwidth")
# # # Figure 4
# # doPreprocess = True
# # cur_bandwidth = defaultBandWidth
# # cur_latency = defaultLatency
# # if isCluster:
# #     setup_network(cur_bandwidth, cur_latency)
# # for executable in executableList:
# #     executable_path = executable_root_path + executable
# #     run_latency(executable, executable_path, cur_bandwidth, cur_latency)
# # if isCluster:
# #     clean_network()

# # doPreprocess = False
# # # cur_latency = defaultLatency
# # for cur_latency in latencyList:
# #     for cur_bandwidth in bandwidthList:
# #         if isCluster:
# #             setup_network(cur_bandwidth, cur_latency)
# #         for executable in executableList:
# #             executable_path = executable_root_path + executable
# #             run_latency(executable, executable_path, cur_bandwidth, cur_latency)
# #         if isCluster:
# #             clean_network()

# print(">>>> Measure Communication")
# # Measure Communication
# # With preprocessing
# doPreprocess = True
# cur_bandwidth = defaultBandWidth
# cur_latency = defaultLatency
# if isCluster:
#     setup_network(cur_bandwidth, cur_latency)

# executable_path = executable_root_path + defaultExecutable
# run_graph_scale(defaultExecutable, executable_path, bandwidthList[-1], latencyList[0])

# if isCluster:
#     clean_network()

# # Without preprocessing
# doPreprocess = False
# cur_bandwidth = defaultBandWidth
# cur_latency = defaultLatency
# if isCluster:
#     setup_network(cur_bandwidth, cur_latency)

# executable_path = executable_root_path + defaultExecutable
# run_graph_scale(defaultExecutable, executable_path, bandwidthList[-1], latencyList[0])

# if isCluster:
#     clean_network()

print(">>>> diff rotation-based and round-halved")
# rotation-based vs round-halved

doPreprocess = True
cur_bandwidth = defaultBandWidth
cur_latency = defaultLatency
if isCluster:
    setup_network(cur_bandwidth, cur_latency)

executable_path = executable_root_path + defaultExecutable
run_rotation_based(defaultExecutable, executable_path, defaultBandWidth, defaultLatency)

if isCluster:
    clean_network()

doPreprocess = False
cur_bandwidth = defaultBandWidth
cur_latency = defaultLatency
if isCluster:
    setup_network(cur_bandwidth, cur_latency)

executable_path = executable_root_path + defaultExecutable
run_rotation_based(defaultExecutable, executable_path, defaultBandWidth, defaultLatency)

if isCluster:
    clean_network()

