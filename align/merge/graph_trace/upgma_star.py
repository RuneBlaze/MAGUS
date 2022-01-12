import heapq
from collections import deque 

from configuration import Configs

'''
This involves mostly copy-pasting of Vlad's code.
Hopefully this is within Vlad's permission to change it, since there is no explicit license to MAGUS.

If I cannot beat Vlad's code, my code will join Vlad's code.
'''

def upgmaSearch(graph, zeroWeight = False):
    Configs.log("Finding graph trace with upgma, or actually a modified version of region-growing search..")
    
    k = len(graph.context.subalignments)
    lowerBound = [graph.subsetMatrixIdx[i] for i in range(k)]
    upperBound = [graph.subsetMatrixIdx[i] + graph.subalignmentLengths[i] for i in range(k)] 
    graph.clusters = upgmaCluster(graph, lowerBound, upperBound, zeroWeight)

def upgmaCluster(graph, lowerBound, upperBound, zeroWeight = False):
    clusters = []
    clusterPointers = {}
    clusterPos = {}
    nodeClusters = {}
    weightMap = []
    absorbed = set()
    cantConnects = set()
    enforceTrace = True
    
    for s in range(len(lowerBound)):
        for a in range(lowerBound[s], upperBound[s]):
            clusters.append([a])
            idx = len(clusters)-1
            nodeClusters[a] = idx
            weightMap.append({})
            clusterPos[idx] = {s : a}
            clusterPointers[idx] = {s : (idx-1 if idx > lowerBound[s] else None, idx+1 if idx < upperBound[s]-1 else None)}
    
    
    heap = buildHeap(graph, nodeClusters, weightMap, lowerBound, upperBound)
    Configs.log("Built a heap of size {}..".format(len(heap)))
    crunchHeap(graph, heap, clusters, nodeClusters, clusterPos, clusterPointers, weightMap, cantConnects, absorbed, zeroWeight)
        
    # if enforceTrace:
    #     clusters = orderClusters(graph, clusters, nodeClusters, lowerBound, upperBound)
    return clusters

def buildHeap(graph, nodeClusters, weightMap, lowerBound, upperBound):
    heap = []
    for s in range(len(lowerBound)):
        for a in range(lowerBound[s], upperBound[s]):
            asub, apos = graph.matSubPosMap[a]
            i = nodeClusters[a]
            for b, value in graph.matrix[a].items():
                bsub, bpos = graph.matSubPosMap[b]
                if b <= a or asub == bsub or b < lowerBound[bsub] or b >= upperBound[bsub]:
                    continue
                j = nodeClusters[b]
                weightMap[j][i] = value
                weightMap[i][j] = value
                heapq.heappush(heap, (-1 * value, a, b))
    
    return heap

def crunchHeap(graph, heap, clusters, nodeClusters, clusterPos, clusterPointers, weightMap, cantConnects, absorbed, zeroWeight):
    enforceTrace = True
    while len(heap) > 0:
        value, a, b = heapq.heappop(heap)
        i, j = nodeClusters[a], nodeClusters[b]
        if i == j or orderPair(i,j) in cantConnects:
            continue
        if value != -1 * weightMap[i][j]:
            # this is a hack, but it should work. If the value is not up-to-date, we throw it out
            continue
        
        if not checkConnect(graph, i, j, clusters, clusterPos, enforceTrace):
            cantConnects.add(orderPair(i,j))
            continue

        size_i = len(clusters[i]) # the size of i before absorbing j
        size_j = len(clusters[j]) # the size of j before being absorbed
        
        absorbed.add(j)
        for e in clusters[j]:
            nodeClusters[e] = i
            clusters[i].append(e)
            asub, apos = graph.matSubPosMap[e]
            clusterPos[i][asub] = e
        clusters[j] = []

        if enforceTrace:
            for s in clusterPointers[j]:
                prev, nxt = clusterPointers[j][s]
                if prev is not None:
                    clusterPointers[prev][s] = (clusterPointers[prev][s][0], i) 
                if nxt is not None:
                    clusterPointers[nxt][s] = (i, clusterPointers[nxt][s][1])
                clusterPointers[i][s] = (prev, nxt)

            updateMergePointers(graph, i, clusterPointers, clusters, clusterPos)
         
        # print("Clusters left: {}".format(len(clusters) - len(absorbed)))
        for n in set(weightMap[i].keys()) | set(weightMap[j].keys()):
            if n in absorbed:
                continue
            # we do UPGMA star here
            if zeroWeight:
                weightMap[i][n] = (weightMap[i].get(n, 0) * size_i + weightMap[j].get(n, 0) * size_j)/(size_i + size_j)
            else:
                i_connected = n in weightMap[i]
                j_connected = n in weightMap[j]
                if i_connected and j_connected:
                    weightMap[i][n] = (weightMap[i][n] * size_i + weightMap[j][n] * size_j)/(size_i + size_j)
                elif i_connected:
                    weightMap[i][n] = weightMap[i][n]
                elif j_connected:
                    weightMap[i][n] = weightMap[j][n]
                else:
                    raise Exception("This should not happen. Absurd!")
            weightMap[n][i] = weightMap[i][n]
            heapq.heappush(heap, (-1 * weightMap[i][n], clusters[i][0], clusters[n][0]))

        # for n in weightMap[j]:
        #     if n in absorbed:
        #         continue
        #     weightMap[i][n] = weightMap[i].get(n, 0) + weightMap[j][n]
        #     weightMap[n][i] = weightMap[i][n]
        #     heapq.heappush(heap, (-1 * weightMap[i][n], clusters[i][0], clusters[n][0]))

def updateMergePointers(graph, i, clusterPointers, clusters, clusterPos):
    subsets = [graph.matSubPosMap[a][0] for a in clusters[i]]
    
    for s in subsets:
        queue = deque([i])
        visited = set([i])
        
        while len(queue) > 0:
            curNode = queue.popleft()
            
            if clusterPos[curNode].get(s, float('inf')) > clusterPos[i][s] or curNode == i:
                clusterPos[curNode][s] = clusterPos[i][s]
                
                for p in clusterPointers[curNode]:
                    prv, nxt = clusterPointers[curNode][p]
                    if prv not in visited and prv is not None:
                        queue.append(prv)
                        visited.add(prv)
           
            
def checkConnect(graph, i, j, clusters, clusterPos, enforceTrace):
    ci , cj = set([graph.matSubPosMap[a][0] for a in clusters[i]]), set([graph.matSubPosMap[a][0] for a in clusters[j]])
    for s in ci:
        if s in cj:
            return False
    
    if not enforceTrace:
        return True
    
    for s in ci:
        if clusterPos[j].get(s, float('inf')) <= clusterPos[i][s]:
            return False
    for s in cj:
        if clusterPos[i].get(s, float('inf')) <= clusterPos[j][s]:
            return False    
    return True

def orderClusters(graph, clusters, nodeClusters, lowerBound, upperBound):
    orderedClusters = []
    frontier = list(lowerBound)
    while True:
        foundGood = False
        for j in range(len(lowerBound)):
            good = True
            idx = frontier[j]
            if idx >= upperBound[j]:
                continue
            i = nodeClusters[idx]
            for b in clusters[i]:
                bsub, bpos = graph.matSubPosMap[b]
                if b > frontier[bsub]:
                    #print(bsub, b, frontier[bsub])
                    good = False
                    break
                
            if good:
                orderedClusters.append(clusters[i])
                for b in clusters[i]:
                    bsub, bpos = graph.matSubPosMap[b]
                    frontier[bsub] = b + 1
                foundGood = True
                break
        if not foundGood:
            break
    return orderedClusters

def orderPair(a, b):
    return (min(a, b), max(a, b))