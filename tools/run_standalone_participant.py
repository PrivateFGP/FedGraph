import subprocess
import sys
import os

executable_root_path = "./../../bin/"
data_root_path = "./data/"
result_root_path = "./result/"
log_root_path = "./log/"
preprocess_root_path = "./preprocess/"
iterations = 2
defaultNumParts = 2
defaultInterRatio = 0.4
doPreprocess = True

def my_makedir(path):
    if not os.path.isdir(path):
        os.makedirs(path)

# print(process.stdout.read().decode('utf-8'))

# python batch_Deckard.py > test_case/batch_similarity.txt

# ./sssp-ss -t 5 -g 5 -i 4 -m 2 -p 1 ./../data/test.dat ./../data/test5p.part ./../data/result_4.txt 0

processList = []

# # Party num setting
# numPartsList = [2,3,4,5]
# data_path = data_root_path + "party_num/"
# result_path = result_root_path + "party_num/"
# log_path = log_root_path + "party_num/"
# my_makedir(result_path)
# my_makedir(log_path)
# for num in numPartsList:
#     setting = executable+"_"+str(num)+"p"
#     preprocess_path = preprocess_root_path + setting + "/"
#     my_makedir(preprocess_path)
#     for i in range(num):
#         cmd = [executable_path, "-t", str(num), "-g", str(num), "-i", str(i), "-m", str(iterations), "-p", "1", "-s", setting]
#         if not doPreprocess:
#             cmd += ["-n", "1"]
#         cmd += [data_path+str(num)+"p.dat", data_path+str(num)+"p.part",
#                 result_path+setting+"_"+str(i)+".result"]
#         print(" ".join(cmd))
#         log_f = open(log_path+setting+"_"+str(i)+".log", 'w', encoding='utf-8')
#         processList.append(subprocess.Popen(cmd, stdout=log_f))
#     for process in processList:
#         process.wait()
#     processList = []

# # Inter-edge ratio setting
# # interRatioList = [0.05, 0.1, 0.2, 0.4, 0.8]
# interRatioList = [0.2, 0.4, 0.8]
# data_path = data_root_path + "inter_ratio/"
# result_path = result_root_path + "inter_ratio/"
# log_path = log_root_path + "inter_ratio/"
# my_makedir(result_path)
# my_makedir(log_path)
# for ratio in interRatioList:
#     setting = executable+"_"+str(ratio)+"inter"
#     preprocess_path = preprocess_root_path + setting + "/"
#     my_makedir(preprocess_path)
#     for i in range(defaultNumParts):
#         cmd = [executable_path, "-t", str(defaultNumParts), "-g", str(defaultNumParts), "-i", str(i), "-m", str(iterations), "-p", "1", "-s", setting]
#         if not doPreprocess:
#             cmd += ["-n", "1"]
#         cmd += [data_path+str(ratio)+"inter.dat", data_path+str(ratio)+"inter.part",
#                 result_path+setting+"_"+str(i)+".result"]
#         print(" ".join(cmd))
#         log_f = open(log_path+setting+"_"+str(i)+".log", 'w', encoding='utf-8')
#         processList.append(subprocess.Popen(cmd, stdout=log_f))
#     for process in processList:
#         process.wait()
#     processList = []

# Evaluation root setting
root_setting = "graph_scale"

# Graph scale
# scalerList = [0.2, 0.5, 1, 2, 5]
def run_graph_scale(executable, executable_path):
    scalerList = [1, 2, 5, 10]
    data_path = data_root_path + "graph_scale/"
    result_path = result_root_path + "graph_scale/"
    log_path = log_root_path + "graph_scale/"
    my_makedir(result_path)
    my_makedir(log_path)
    for scaler in scalerList:
        preprocess_setting = root_setting+"_"+executable+"_"+str(scaler)+"scale_"+str(defaultNumParts)+"p_"+str(defaultInterRatio)+"inter"
        eval_setting = preprocess_setting
        preprocess_path = preprocess_root_path + preprocess_setting + "/"
        my_makedir(preprocess_path)
        processList = []
        for i in range(defaultNumParts):
            cmd = [executable_path, "-t", str(defaultNumParts), "-g", str(defaultNumParts), "-i", str(i), "-m", str(iterations), "-p", "1", "-s", preprocess_setting]
            if not doPreprocess:
                cmd += ["-n", "1"]
            cmd += [data_path+str(scaler)+"scale.dat", data_path+str(scaler)+"scale.part",
                    result_path+eval_setting+"_"+str(i)+".result"]
            print(" ".join(cmd))
            log_f = open(log_path+eval_setting+"_"+str(i)+".log", 'w', encoding='utf-8')
            processList.append(subprocess.Popen(cmd, stdout=log_f))
        for process in processList:
            process.wait()
        processList = []

executableList = ["bfs-ss", "sssp-ss", "cc-ss", "pagerank-ss"]
for executable in executableList:
    executable_path = executable_root_path + executable
    run_graph_scale(executable, executable_path)