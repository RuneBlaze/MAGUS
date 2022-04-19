'''
Created on Jul 28, 2021

@author: Vlad
'''

import os
import random
import heapq
import math

from helpers import sequenceutils
from tasks import task
from configuration import Configs
from tools import external_tools
import itertools


def requestMafftBackbones(context):
    missingBackboneFiles = {}
    numRuns = Configs.mafftRuns
    strategy = Configs.graphBuildStrategy.lower()
    if strategy.startswith("mst"):
        k = int(strategy[3])
        numRuns = (len(context.subsets) - 1) * k
    for n in range(numRuns):
        unalignedFile = os.path.join(context.graph.workingDir, "backbone_{}_unalign.txt".format(n+1))
        alignedFile = os.path.join(context.graph.workingDir, "backbone_{}_mafft.txt".format(n+1))
        if os.path.exists(alignedFile):
            Configs.log("Existing backbone file found: {}".format(alignedFile))            
            context.backbonePaths.append(alignedFile)
        else:
            missingBackboneFiles[unalignedFile] = alignedFile
            
        if Configs.graphBuildHmmExtend:
            context.backboneExtend.add(alignedFile)
    assignBackboneTaxa(context, missingBackboneFiles)
    for unalignedFile, alignedFile in missingBackboneFiles.items():
        if Configs.graphBuildMethod == "mafft":
            backboneTask = external_tools.buildMafftAlignment(unalignedFile, alignedFile, useFafft=Configs.fragAlignForBb)
        elif Configs.graphBuildMethod == "magus":
            backboneTask = external_tools.runMyself(unalignedFile, os.path.join(alignedFile + "_env"), alignedFile, None)
        else:
            assert False
        context.backboneTasks.append(backboneTask)
    
    if not Configs.graphBuildHmmExtend:
        for file in list(missingBackboneFiles.keys()) + context.backbonePaths:
            backbone = sequenceutils.readFromFasta(file)
            context.backboneTaxa.update(backbone)        
            
    task.submitTasks(context.backboneTasks)    
    
def assignBackboneTaxa(context, missingBackbones):
    if len(missingBackbones) == 0:
        return
    
    numTaxa = max(1, int(Configs.mafftSize/len(context.subsetPaths)))
    backbones = {file : {} for file in missingBackbones}

    strat = Configs.graphBuildStrategy.lower()
    if strat.startswith("mst"):
        divider = 2 if strat[-1] == "-" else 1
        buildBackbonesMST(context, backbones, numTaxa, divider)

    if Configs.graphBuildStrategy.lower() == "eligible":
        buildBackbonesEligible(context, backbones, numTaxa)
    
    def eligible_oracle2(b, s):
        if b >= 0.5:
            if s >= 0.9:
                return 'S'
        return 'L'
    def eligible_oracle3(b, s):
        return 'M'
    def eligible_oracle4(b, s):
        if b >= 0.5:
            return 'S'
        return 'L'
    def eligible_oracle5(b, s):
        if s >= 0.67:
            return 'S'
        return 'L'
    def eligible_oracle6(b, s):
        if b >= 0.5:
            if s >= 0.5:
                return 'S'
        return 'L'
    def eligible_oracle7(b, s):
        if b >= 0.5:
            if s >= 0.75:
                return 'S'
        return 'L'
    def eligible_oracle8(b, s):
        if b >= 0.0:
            return 'S'
        return 'L'
    def eligible_oracle9(b, s):
        if b >= 0.25:
            return 'S'
        return 'L'
    
    gbs = Configs.graphBuildStrategy.lower()
    if len(gbs) == len("eligibleX"):
        variant = int(gbs[-1])
        oracle = {
            2: eligible_oracle2,
            3: eligible_oracle3,
            4: eligible_oracle4,
            5: eligible_oracle5,
            6: eligible_oracle6,
            7: eligible_oracle7,
            8: eligible_oracle8,
            9: eligible_oracle9,
        }[variant]
        Configs.log(f"Using oracle {oracle}")
        buildBackbonesEligible(context, backbones, numTaxa, oracle)
    
    if Configs.graphBuildStrategy.lower() == "random":    
        buildBackbonesRandom(context, backbones, numTaxa)
        
    elif Configs.graphBuildStrategy.lower() == "longest":
        buildBackbonesLongest(context, backbones, numTaxa)
    
    elif Configs.graphBuildStrategy.lower() == "longestrandom":
        buildBackbonesLongestRandom(context, backbones, numTaxa)
    
    elif Configs.graphBuildStrategy.lower() == "longestrandom2":
        buildBackbonesLongestRandom2(context, backbones, numTaxa)
    
    elif Configs.graphBuildStrategy.lower() == "longestrandom3":
        buildBackbonesLongestRandom3(context, backbones, numTaxa)
    
    elif Configs.graphBuildStrategy.lower() == "longestrandom4":
        buildBackbonesLongestRandom4(context, backbones, numTaxa)
    
    elif Configs.graphBuildStrategy.lower() == "longestrandom5":
        buildBackbonesLongestRandom5(context, backbones, numTaxa)
    
    elif Configs.graphBuildStrategy.lower() == "longestrandom6":
        buildBackbonesLongestRandom6(context, backbones, numTaxa)
    
    elif Configs.graphBuildStrategy.lower() == "longestrandom7":
        buildBackbonesLongestRandom7(context, backbones, numTaxa)

    elif Configs.graphBuildStrategy.lower() == "coverage":
        buildBackbonesCoverage(context, backbones, numTaxa)  
    
    elif Configs.graphBuildStrategy.lower() == "coverage2":
        buildBackbonesCoverage2(context, backbones, numTaxa)                          
                    
    for file, backbone in backbones.items():
        sequenceutils.writeFasta(backbone, file)
    
    
def buildBackbonesRandom(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} RANDOM sequences per subset..".format(len(backbones), numTaxa))
    for file, backbone in backbones.items():
        if backbone:
            continue
        for subset in context.subsets:
            random.shuffle(subset)
            for taxon in subset[:numTaxa]:
                backbone[taxon] = context.unalignedSequences[taxon]       

def buildBackbonesMST(context, backbones, numTaxa, divider = 1):
    edges = []
    for n in context.MST.traverse_postorder():
        if n.is_root():
            continue
        if not n.label:
            assert False
        edges.append((int(n.label)-1, int(n.parent.label)-1))
    mst_num = len(backbones) // divider
    for i, (_, backbone) in enumerate(backbones.items()):
        if i >= mst_num:
            continue
        e_ix = i % len(edges) # which edge we are currently during
        u, v = edges[e_ix]
        m = Configs.mafftSize
        random.shuffle(context.subsets[u])
        random.shuffle(context.subsets[v])
        subset_a = context.subsets[u][:(m // 2)]
        subset_b = context.subsets[v][:(m // 2)]
        for taxon in subset_a:
            backbone[taxon] = context.unalignedSequences[taxon]
        for taxon in subset_b:
            backbone[taxon] = context.unalignedSequences[taxon]
    buildBackbonesRandom(context, backbones, numTaxa)
    

def buildBackbonesEligible(context, backbones, numTaxa, oracle = lambda a, b: 'L'):
    sequences = context.unalignedSequences
    Configs.log("Preparing {} backbones with {} ELIGIBLE sequences per subset..".format(len(backbones), numTaxa))
    seqLengths = [len(sequences[t].seq) for t in sequences]
    seqLengths.sort()
    threshold = 0.5 if Configs.skeletonSeqs == "median" else 0.75
    topQuartile = seqLengths[int(threshold*(len(seqLengths)-1))]
    Configs.targetLength = topQuartile
    Configs.log(f"Target length {topQuartile} for eligibility")
    bid = 0
    for file, backbone in backbones.items():
        bpercent = (bid + 1) / len(backbones)
        bb_log = ""
        perm = [y for y, _ in enumerate(context.subsets)]
        random.shuffle(perm)
        for sid, subset in enumerate(context.subsets):
            spercent = (perm[sid] + 1) / len(context.subsets)
            eligible = []
            ineligible = []
            for t in subset:
                if abs(len(sequences[t].seq) - Configs.targetLength) < 0.25 * Configs.targetLength:
                    eligible.append(t)
                else:
                    ineligible.append(t)
            if bid == 0:
                Configs.log(f"Total eligible sequences: {len(eligible)}, ineligible: {len(ineligible)}")
            random.shuffle(eligible)
            random.shuffle(ineligible)
            eligible_fst = eligible + ineligible
            ineligible_fst = ineligible + eligible
            mixed = [x for x in itertools.chain(*itertools.zip_longest(eligible, ineligible)) if x is not None]
            target = {
                'L': eligible_fst,
                'S': ineligible_fst,
                'M': mixed,
            }
            s = oracle(bpercent, spercent)
            bb_log += s
            for taxon in target[s][:numTaxa]:
                backbone[taxon] = context.unalignedSequences[taxon]
        Configs.log(f"Backbone composition: {bb_log}")
        bid += 1

def buildBackbonesLongest(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} LONGEST sequences per subset..".format(len(backbones), numTaxa))
    for subset in context.subsets:
        sortedByLength = sorted(subset, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        curIdx = 0
        for file, backbone in backbones.items():
            window = sortedByLength[curIdx : curIdx + numTaxa]
            if curIdx + numTaxa > len(sortedByLength):
                window.extend(sortedByLength[0 : min(curIdx, curIdx + numTaxa - len(sortedByLength))])
            for taxon in window:
                backbone[taxon] = context.unalignedSequences[taxon]
            curIdx = (curIdx + numTaxa) % len(sortedByLength)

def buildBackbonesLongestRandom(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} LONGEST (randomized) sequences per subset..".format(len(backbones), numTaxa))
    for subset in context.subsets:
        sortedByLength = sorted(subset, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        topSequences = sortedByLength[ : max(numTaxa, int(0.25*len(sortedByLength)))]
        for file, backbone in backbones.items():
            random.shuffle(topSequences)
            for taxon in topSequences[:numTaxa]:
                backbone[taxon] = context.unalignedSequences[taxon]

def buildBackbonesLongestRandom2(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} LONGEST (randomized-2) sequences per subset..".format(len(backbones), numTaxa))
    for subset in context.subsets:
        sortedByLength = sorted(subset, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        topSequences = sortedByLength[ : numTaxa * 2]
        for file, backbone in backbones.items():
            random.shuffle(topSequences)
            for taxon in topSequences[:numTaxa]:
                backbone[taxon] = context.unalignedSequences[taxon]

def buildBackbonesLongestRandom3(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} LONGEST (randomized-3) sequences per subset..".format(len(backbones), numTaxa))
    for subset in context.subsets:
        sortedByLength = sorted(subset, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        topSequences = sortedByLength[ : numTaxa * 3]
        for file, backbone in backbones.items():
            random.shuffle(topSequences)
            for taxon in topSequences[:numTaxa]:
                backbone[taxon] = context.unalignedSequences[taxon]

def buildBackbonesLongestRandom4(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} LONGEST (randomized-4) sequences per subset..".format(len(backbones), numTaxa))
    for subset in context.subsets:
        sortedByLength = sorted(subset, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        numTop = math.ceil(numTaxa/2)
        topSeq, remainder = sortedByLength[ : numTop], sortedByLength[numTop : ]
        for file, backbone in backbones.items():
            random.shuffle(remainder)
            for taxon in topSeq + remainder[ : numTaxa-numTop]:
                backbone[taxon] = context.unalignedSequences[taxon]

def buildBackbonesLongestRandom5(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} LONGEST (randomized-5) sequences per subset..".format(len(backbones), numTaxa))
    for subset in context.subsets:
        sortedByLength = sorted(subset, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        numTop = math.ceil(numTaxa/4)
        topSeq, remainder = sortedByLength[ : numTop], sortedByLength[numTop : ]
        for file, backbone in backbones.items():
            random.shuffle(remainder)
            for taxon in topSeq + remainder[ : numTaxa-numTop]:
                backbone[taxon] = context.unalignedSequences[taxon]

def buildBackbonesLongestRandom6(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} LONGEST (randomized-6) sequences per subset..".format(len(backbones), numTaxa))
    for subset in context.subsets:
        sortedByLength = sorted(subset, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        #topQuartileTaxon = sortedByLength[int(0.25*len(sortedByLength))]
        topLength = len(context.unalignedSequences[sortedByLength[0]].seq)
        fullLength = []
        notFullLength = []
        for t in sortedByLength:
            if abs(len(context.unalignedSequences[t].seq) - topLength) < 0.25 * topLength:
                fullLength.append(t)
            else:
                notFullLength.append(t) 
        
        Configs.log("Found {}/{} full-length sequences..".format(len(fullLength), len(subset)))        
        for file, backbone in backbones.items():
            random.shuffle(fullLength)
            random.shuffle(notFullLength)  
            allTaxa = fullLength + notFullLength  
            for taxon in allTaxa[:numTaxa]:
                backbone[taxon] = context.unalignedSequences[taxon]

def buildBackbonesLongestRandom7(context, backbones, numTaxa):
    Configs.log("Preparing {} backbones with {} LONGEST (randomized-7) sequences per subset..".format(len(backbones), numTaxa))
    for subset in context.subsets:
        sortedByLength = sorted(subset, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        #topQuartileTaxon = sortedByLength[int(0.25*len(sortedByLength))]
        topLength = len(context.unalignedSequences[sortedByLength[0]].seq)
        fullLength = []
        notFullLength = []
        for t in sortedByLength:
            if abs(len(context.unalignedSequences[t].seq) - topLength) < 0.5 * topLength:
                fullLength.append(t)
            else:
                notFullLength.append(t) 
        
        Configs.log("Found {}/{} full-length sequences..".format(len(fullLength), len(subset)))        
        for file, backbone in backbones.items():
            random.shuffle(fullLength)
            random.shuffle(notFullLength)  
            allTaxa = fullLength + notFullLength  
            for taxon in allTaxa[:numTaxa]:
                backbone[taxon] = context.unalignedSequences[taxon]

def buildBackbonesCoverage(context, backbones, numTaxa):
    context.awaitSubalignments()
    Configs.log("Preparing {} backbones with {} COVERAGE sequences per subset..".format(len(backbones), numTaxa))
        
    for subalignPath in context.subalignmentPaths:
        subalignment = sequenceutils.readFromFasta(subalignPath, removeDashes=False)
        taxa = list(subalignment.keys())
        sortedByLength = sorted(taxa, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        topSequences = sortedByLength[ : numTaxa*2]
        
        coverArray = None
        coverageDict = {}
        for taxon in sortedByLength:
            sequence = subalignment[taxon]
            if coverArray is None:
                coverArray = [[] for i in range(len(sequence.seq))]
            coverageDict[taxon] = []
            for i, c in enumerate(sequence.seq):
                if c not in ('-', '.', '_'):
                    coverArray[i].append(taxon)
                    coverageDict[taxon].append(i)
        
        
        for file, backbone in backbones.items():
            heap = []
            usedTaxons = set()
            coverage = [0] * len(coverArray)
            coverArrayPointers = [0] * len(coverArray)
            
            random.shuffle(topSequences)
            numSeq = int(numTaxa / 2)
            for taxon in topSequences[:numSeq]:
                usedTaxons.add(taxon)
                backbone[taxon] = context.unalignedSequences[taxon]
                for pos in coverageDict[taxon]:
                    coverage[pos] = coverage[pos] + 1
                
                
            for i in range(len(coverArray)):
                for p in range(len(coverArray[i])):
                    taxon = coverArray[i][p]
                    if not taxon in usedTaxons:
                        coverArrayPointers[i] = p
                        item = (coverage[i], -1 * len(coverageDict[taxon]), random.random(), i, taxon)
                        heapq.heappush(heap, item)
                        break
            
            while len(usedTaxons) < min(numTaxa, len(subalignment)):
                colCover, taxCover, trash, i, taxon = heapq.heappop(heap)
                if taxon in usedTaxons:
                    continue
                if colCover < coverage[i]:
                    item = (coverage[i], -1 * len(coverageDict[taxon]), random.random(), i, taxon)
                    heapq.heappush(heap, item)
                    continue
                
                usedTaxons.add(taxon)
                backbone[taxon] = context.unalignedSequences[taxon]                
                for pos in coverageDict[taxon]:
                    coverage[pos] = coverage[pos] + 1
                    
                    if taxon == coverArray[pos][coverArrayPointers[pos]]:
                        while coverArrayPointers[pos] < len(coverArray[pos]) and coverArray[pos][coverArrayPointers[pos]] in usedTaxons:
                            coverArrayPointers[pos] = coverArrayPointers[pos] + 1
                        if coverArrayPointers[pos] < len(coverArray[pos]):
                            newTaxon = coverArray[pos][coverArrayPointers[pos]]
                            item = (coverage[pos], -1 * len(coverageDict[newTaxon]), random.random(), pos, newTaxon)
                            heapq.heappush(heap, item) 

def buildBackbonesCoverage2(context, backbones, numTaxa):
    context.awaitSubalignments()
    Configs.log("Preparing {} backbones with {} COVERAGE2 sequences per subset..".format(len(backbones), numTaxa))
        
    for subalignPath in context.subalignmentPaths:
        subalignment = sequenceutils.readFromFasta(subalignPath, removeDashes=False)
        taxa = list(subalignment.keys())
        sortedByLength = sorted(taxa, key = lambda t : len(context.unalignedSequences[t].seq), reverse = True)
        
        coverArray = None
        coverageDict = {}
        for taxon in sortedByLength:
            sequence = subalignment[taxon]
            if coverArray is None:
                coverArray = [[] for i in range(len(sequence.seq))]
            coverageDict[taxon] = []
            for i, c in enumerate(sequence.seq):
                if c not in ('-', '.', '_'):
                    coverArray[i].append(taxon)
                    coverageDict[taxon].append(i)
        
        
        for file, backbone in backbones.items():
            heap = []
            heapTaxons = set()
            usedTaxons = set()
            usedSites = set()
             
            for i in range(len(coverArray)):
                for p in range(len(coverArray[i])):
                    taxon = coverArray[i][p]
                    if taxon not in heapTaxons:
                        item = (-1 * len(coverageDict[taxon]), random.random(), taxon)
                        heapTaxons.add(taxon)
                        heapq.heappush(heap, item)
                        break
            
            while len(usedTaxons) < min(numTaxa, len(subalignment)) and len(heap) > 0:
                taxCover, trash, taxon = heapq.heappop(heap)   
                for pos in coverageDict[taxon]:
                    if pos not in usedSites:
                        usedSites.add(pos)
                        usedTaxons.add(taxon)
                        backbone[taxon] = context.unalignedSequences[taxon]                
            Configs.log("Selected {} coverage sequences..".format(len(usedTaxons)))
            
            remainder = [t for t in taxa if t not in usedTaxons]
            random.shuffle(remainder)
            for taxon in remainder[ : numTaxa - len(usedTaxons)]:
                backbone[taxon] = context.unalignedSequences[taxon]