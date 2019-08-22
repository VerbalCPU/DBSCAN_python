import sys
import math


def myDBscan(file_name, eps, minPts):
    with open(file_name, 'r') as data:
        x = []
        y = []
        result = []

        for line in data:
            p = line.split(";")
            x.append(float(p[0]))
            y.append(float(p[1]))
            result = list(zip(x,y))
        result.pop(0)
        return compute_clusters(result, eps, minPts)


def euclidean0_1(vector1, vector2):
    dist = [(a - b)**2 for a, b in zip(vector1, vector2)]
    dist = math.sqrt(sum(dist))
    return dist


def regionQuery(data,gene_data, eps):
    neighbors=[]

    for i in range(len(data)):
        if euclidean0_1(data[gene_data],data[i]) <= eps:
            neighbors.append(i)
    return neighbors


def expandCluster(data, gene_data, neighborPts, eps, minPts, C, label):
    label[gene_data] = C
    i = 0
    while i < len(neighborPts):
        Pt = neighborPts[i]

        if (label[Pt] == 0):
            label[Pt] = C
        if (label[Pt] == "unvisited"):
            label[Pt] = C
            new_neighborPts = list(regionQuery(data, Pt, eps))
            if (len(new_neighborPts) >= minPts):
                neighborPts = list((neighborPts + new_neighborPts))
        i = i + 1

def compute_clusters(data,eps,minPts):
    C=0
    label=list(data)
    for i in range(len(label)):
        label[i]="unvisited"
    for i in range(len(data)):
        gene_data=i
        if(label[i]=="unvisited" or label[i]==-1):
            neighborPts = list(regionQuery(data,gene_data, eps))
            if(len(neighborPts) < minPts):
                label[i]= 0
            else:
                C=C+1
                expandCluster(data,gene_data,neighborPts, eps, minPts,C, label)
    return ' '.join(map(str, label))


output = myDBscan(sys.argv[1],float(sys.argv[2]),int(sys.argv[3]))
print(output)



