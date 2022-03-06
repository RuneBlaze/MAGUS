'''
Created on May 28, 2020

@author: Vlad
'''

import os
import random
import time

from align.decompose import initial_tree, kmh
from helpers import treeutils, sequenceutils
from configuration import Configs

'''
Handles the different ways to decompose the dataset into subsets.
The main way is to estimate a guide tree, then use PASTA's centroid edge decomposition
on the guide tree. Can also decompose randomly (for high speed on huge datasets).
'''

def decomposeSequences(context):
    time1 = time.time()
    
    if len(context.subsetPaths) > 0:
        Configs.log("Subset paths already provided, skipping decomposition..")
    
    elif len(context.subalignmentPaths) > 0:
        context.subsetPaths = context.subalignmentPaths
        Configs.log("Subalignment paths already provided, skipping decomposition..")
    
    else:
        subsetsDir = os.path.join(context.workingDir, "decomposition")
        context.subsetPaths = []
        n = 1
        while True:
            filePath = os.path.join(subsetsDir, "subset_{}.txt".format(n))
            if not os.path.exists(filePath):
                break
            Configs.log("Detected existing subset file {}".format(filePath))
            context.subsetPaths.append(filePath)
            n = n + 1
        
        if len(context.subsetPaths) == 0:
            buildDecomposition(context, subsetsDir)
    
    time2 = time.time()  
    Configs.log("Decomposed {} into {} subsets in {} sec..".format(context.sequencesPath, len(context.subsetPaths), time2-time1))

def buildDecomposition(context, subsetsDir):  
    if not os.path.exists(subsetsDir):
        os.makedirs(subsetsDir)  
    if context.unalignedSequences is None:
        context.unalignedSequences = sequenceutils.readFromFasta(context.sequencesPath, removeDashes=True)    
    
    if (Configs.decompositionStrategy == "random" or context.guideTree == "random") and Configs.outputPath == context.outputFile:
        context.subsetPaths = randomDecomposition(subsetsDir, context.unalignedSequences, Configs.decompositionMaxNumSubsets)
        
    elif Configs.decompositionStrategy == "kmh":
        Configs.log("Decomposing {} with KMH..".format(context.sequencesPath))
        Configs.log("Targetting {} subsets..".format(Configs.decompositionMaxNumSubsets))
        context.subsetPaths = kmh.buildSubsetsKMH(context, subsetsDir)
    
    else:
        guideTreePath  = initial_tree.buildInitialTree(context, subsetsDir, context.guideTree)
        Configs.log("Using target subset size of {}, and maximum number of subsets {}..".format(Configs.decompositionMaxSubsetSize, Configs.decompositionMaxNumSubsets))
        weightSet = set([])
        if Configs.decompositionStrategy == "noodle":
            Configs.log("Total num. of full length sequences: {}".format(len(context.fullSequences)))
            weightSet = context.fullSequences
        context.subsetPaths = treeutils.decomposeGuideTree(subsetsDir, context.sequencesPath, guideTreePath, 
                                                   Configs.decompositionMaxSubsetSize, Configs.decompositionMaxNumSubsets, weightSet)        

def pasta_shim(n):
    return n.replace("_","").replace("/","").replace("-","").lower()

def guess_isfasta(p):
    # check if p is a file and it has a '>' in its first line
    if os.path.isfile(p):
        with open(p, 'r') as f:
            return f.readline().startswith('>')
    else:
        return False

def chooseSkeletonTaxa(sequences, skeletonSize, mode = "fulllength", maximalist = False, context = None):
    Configs.log("Choosing {} taxa for skeleton..".format(skeletonSize))
    allTaxa = list(sequences.keys())
    if Configs.skeletonSeqs:
        if guess_isfasta(Configs.skeletonSeqs):
            Configs.log("Using provided skeleton sequences.. from {}".format(Configs.skeletonSeqs))
            seqs = sequenceutils.readFromFasta(Configs.skeletonSeqs, removeDashes=True)
            rawSeqs = set(seqs.keys())
            skeletonTaxa = set([t for t in allTaxa if t in rawSeqs or pasta_shim(t) in rawSeqs])
            assert skeletonTaxa.issubset(allTaxa)
            assert len(skeletonTaxa) == len(rawSeqs)
            return list(skeletonTaxa), list(set(allTaxa) - skeletonTaxa)
        else:
            mode = Configs.skeletonSeqs
            Configs.log("Skeleton strategy: '{}'".format(mode))
    
    if mode == "fulllength" or mode == "median":
        seqLengths = [len(sequences[t].seq) for t in sequences]
        #topQuartile = numpy.quantile(seqLengths, 0.75)
        seqLengths.sort()
        threshold = 0.5 if mode == "median" else 0.75
        Configs.log(f"Using {mode} threshold of {threshold}")
        topQuartile = seqLengths[int(threshold*(len(seqLengths)-1))]
        
        fullLength = []
        notFullLength = []
        Configs.targetLength = topQuartile
        for t in allTaxa:
            if abs(len(sequences[t].seq) - topQuartile) < 0.25 * topQuartile:
                fullLength.append(t)
            else:
                notFullLength.append(t) 
        if context:
            #TODO: this will not work when MAGUS restarts
            Configs.log("Wrote full and fragged taxa set to context")
            context.fullSequences = set(fullLength)
            context.fragSequences = set(notFullLength)
        Configs.log(f"{len(fullLength)} full length taxa, {len(notFullLength)} fragmented taxa")
        random.shuffle(fullLength)
        random.shuffle(notFullLength)
        allTaxa = fullLength + notFullLength
        if maximalist:
            return fullLength[:skeletonSize], fullLength[skeletonSize:], notFullLength
    else:
        random.shuffle(allTaxa)
        
    skeletonTaxa = allTaxa[:skeletonSize]
    remainingTaxa = allTaxa[skeletonSize:]
    return skeletonTaxa, remainingTaxa

def randomDecomposition(subsetsDir, sequences, numSubsets):
    allTaxa = list(sequences.keys())
    random.shuffle(allTaxa)
    
    taxonSubsets = [allTaxa[i :: numSubsets] for i in range(numSubsets)]
    subsetPaths = []
    for n, subset in enumerate(taxonSubsets):
        subsetPath = os.path.join(subsetsDir, "subset_{}.txt".format(n+1))
        subsetPaths.append(subsetPath)                    
        sequenceutils.writeFasta(sequences, subsetPath, subset) 
    return subsetPaths
