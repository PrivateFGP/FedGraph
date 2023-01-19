original_graph_file_path = "/home/zzh/project/FedGraph/MultipartyPSI/data/account_to_account"
formatted_graph_file_path = "/home/zzh/project/FedGraph/MultipartyPSI/data/large_test.dat"

edgeList = []

with open(original_graph_file_path, 'r', encoding='utf-8') as ifile:
    lines = ifile.readlines()
    for line in lines:
        line = line.strip("\n")
        segs = line.split(",")
        segs = [x for x in segs if x]
        edgeList.append([segs[0], segs[1]])

with open(formatted_graph_file_path, 'w', encoding='utf-8') as ofile:
    for edge in edgeList:
        ofile.write(edge[0] + " " + edge[1] + "\n")