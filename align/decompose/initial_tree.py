'''
Created on May 29, 2020

@author: Vlad
'''

import os
import shutil
import time
import random

from align.decompose import decomposer
from helpers import sequenceutils, hmmutils, treeutils
from tasks import task
from tools import external_tools
from configuration import Configs

'''
Different options for estimating a guide tree.
The main ones are FastTree (for accuracy) and Clustal Omega's mbed (for speed).
'''

def buildInitialTree(context, workingDir, treeType):
    if treeType is not None and os.path.exists(treeType):
        Configs.log("Found user guide tree {}".format(treeType))
        return treeType
    
    tempDir = os.path.join(workingDir, "initial_tree")
    outputTreePath = os.path.join(tempDir, "initial_tree.tre")
    if os.path.exists(outputTreePath):
        Configs.log("Found existing initial tree {}".format(outputTreePath))
        return outputTreePath
    if os.path.exists(tempDir):
        shutil.rmtree(tempDir)
    os.makedirs(tempDir)
    
    time1 = time.time() 
    
    if treeType is None or treeType.lower() == "fasttree": 
        Configs.log("Building PASTA-style FastTree initial tree on {} with skeleton size {}..".format(context.sequencesPath, Configs.decompositionSkeletonSize))
        alignPath = os.path.join(tempDir, "initial_align.txt")
        buildInitialAlignment(context.unalignedSequences, tempDir, Configs.decompositionSkeletonSize, None, alignPath)
        if Configs.outputInitialAlignment:
            buildNaiveAlignment(context.unalignedSequences, tempDir, Configs.outputInitialAlignment)
        external_tools.runFastTree(alignPath, tempDir, outputTreePath, "fast").run()
    elif treeType is None or treeType.lower() == "epa-ng":
         mm = treeType[0] == 'E'
         Configs.log("Building epa-ng initial tree on {} with skeleton size {}..".format(context.sequencesPath, Configs.decompositionSkeletonSize))
         Configs.log("Maximalist mode: {}".format(mm))
         alignPath = os.path.join(tempDir, "initial_align.txt")
         epaNgPipeline(context.unalignedSequences, tempDir, Configs.decompositionSkeletonSize, None, alignPath, maximalist = mm)
    elif treeType is None or treeType.lower() == "fasttree-noml": 
        Configs.log("Building PASTA-style FastTree (NO ML) initial tree on {} with skeleton size {}..".format(context.sequencesPath, Configs.decompositionSkeletonSize))
        alignPath = os.path.join(tempDir, "initial_align.txt")
        buildInitialAlignment(context.unalignedSequences, tempDir, Configs.decompositionSkeletonSize, None, alignPath)
        external_tools.runFastTree(alignPath, tempDir, outputTreePath, "noml").run()
    elif treeType.lower() == "raxml":
        Configs.log("Building RAxML initial tree on {} with skeleton size {}..".format(context.sequencesPath, Configs.decompositionSkeletonSize))
        alignPath = os.path.join(tempDir, "initial_align.txt")
        buildInitialAlignment(context.unalignedSequences, tempDir, Configs.decompositionSkeletonSize, None, alignPath)
        external_tools.runRaxmlNg(alignPath, tempDir, outputTreePath, Configs.numCores).run()
    elif treeType.lower() == "parttree":
        Configs.log("Building MAFFT PartTree initial tree on {}..".format(context.sequencesPath))
        taxa = list(context.unalignedSequences.keys())
        external_tools.runMafftGuideTree(context.sequencesPath, tempDir, outputTreePath, Configs.numCores).run()
        treeutils.convertMafftGuideTree(outputTreePath, taxa)
    elif treeType.lower() == "clustal":
        Configs.log("Building Clustal Omega initial tree on {}..".format(context.sequencesPath))
        external_tools.runClustalOmegaGuideTree(context.sequencesPath, tempDir, outputTreePath, Configs.numCores).run()
    else:
        raise Exception("Guide tree {} not a file and not recognized..".format(treeType))

    time2 = time.time()
    Configs.log("Built initial tree on {} in {} sec..".format(context.sequencesPath, time2-time1))
    
    return outputTreePath

def buildInitialAlignment(sequences, tempDir, skeletonSize, initialAlignSize, outputAlignPath):
    skeletonPath = os.path.join(tempDir, "skeleton_sequences.txt")
    queriesPath = os.path.join(tempDir, "queries.txt") 
    hmmDir = os.path.join(tempDir, "skeleton_hmm")
    hmmPath = os.path.join(hmmDir, "hmm_model.txt")
    initialInsertPath = os.path.join(tempDir, "initial_insert_align.txt")
    if not os.path.exists(hmmDir):
        os.makedirs(hmmDir)
    
    if initialAlignSize is None or initialAlignSize > len(sequences):
        initialAlignSize = len(sequences)
    skeletonTaxa, remainingTaxa = decomposer.chooseSkeletonTaxa(sequences, skeletonSize)
    additional = initialAlignSize-skeletonSize
    random.shuffle(remainingTaxa)
    remainingTaxa, unusedTaxa = remainingTaxa[:additional], remainingTaxa[additional:]
    
    sequenceutils.writeFasta(sequences, skeletonPath, skeletonTaxa)
    external_tools.runMafft(skeletonPath, None, tempDir, outputAlignPath, Configs.numCores).run()
    
    if len(remainingTaxa) > 0:
        sequenceutils.writeFasta(sequences, queriesPath, remainingTaxa)    
        hmmutils.buildHmmOverAlignment(outputAlignPath, hmmPath).run()
        hmmTasks = hmmutils.hmmAlignQueries(hmmPath, queriesPath)
        task.submitTasks(hmmTasks)
        for hmmTask in task.asCompleted(hmmTasks):
            hmmutils.mergeHmmAlignments([hmmTask.outputFile], outputAlignPath, includeInsertions=False)
            if Configs.graphBuildMethod == "initial": # effectively NOP for now
                hmmutils.mergeHmmAlignments([hmmTask.outputFile], initialInsertPath, includeInsertions=True)

def epaNgPipeline(sequences, tempDir, skeletonSize, initialAlignSize, outputAlignPath, maximalist = False):
    skeletonPath = os.path.join(tempDir, "skeleton_sequences.txt")
    queriesPath = os.path.join(tempDir, "queries.txt") 
    hmmDir = os.path.join(tempDir, "skeleton_hmm")
    hmmPath = os.path.join(hmmDir, "hmm_model.txt")
    bb_tree = os.path.join(tempDir, "bb_tree.tre")
    rest_path = os.path.join(tempDir, "rest.fa")
    outputJplacePath = os.path.join(tempDir, "initial_jplace.jplace")
    outputTreePath = os.path.join(tempDir, "initial_tree.tre")
    initialInsertPath = os.path.join(tempDir, "initial_insert_align.txt")
    if not os.path.exists(hmmDir):
        os.makedirs(hmmDir)
    
    if initialAlignSize is None or initialAlignSize > len(sequences):
        initialAlignSize = len(sequences)
    # select the skeleton sequences
    if not maximalist:
        skeletonTaxa, remainingTaxa = decomposer.chooseSkeletonTaxa(sequences, skeletonSize)
        addTaxa = []
    else:
        skeletonTaxa, addTaxa, remainingTaxa = decomposer.chooseSkeletonTaxa(sequences, skeletonSize, maximalist = True)
    # additional = initialAlignSize - skeletonSize
    random.shuffle(remainingTaxa)
    random.shuffle(addTaxa)
    # remainingTaxa, unusedTaxa = remainingTaxa[:additional], remainingTaxa[additional:]
    
    sequenceutils.writeFasta(sequences, skeletonPath, skeletonTaxa)
    # build the skeleton alignment
    external_tools.runMafft(skeletonPath, None, tempDir, outputAlignPath, Configs.numCores).run()
    if len(addTaxa) > 0:
        addHmmPath = os.path.join(hmmDir, "add_hmm_model.txt")
        addQueriesPath = os.path.join(tempDir, "add_queries.txt")
        sequenceutils.writeFasta(sequences, addQueriesPath, addTaxa)
        hmmutils.buildHmmOverAlignment(outputAlignPath, addHmmPath).run()
        addHmmTasks = hmmutils.hmmAlignQueries(addHmmPath, addQueriesPath)
        task.submitTasks(addHmmTasks)
        for addHmmTask in task.asCompleted(addHmmTasks):
            hmmutils.mergeHmmAlignments([addHmmTask.outputFile], outputAlignPath, includeInsertions=False)
    # build an initial tree on the skeleton alignment, called the bb_tree
    external_tools.runFastTree(outputAlignPath, tempDir, bb_tree).run()
    if len(remainingTaxa) > 0:
        sequenceutils.writeFasta(sequences, queriesPath, remainingTaxa)    
        if len(addTaxa) <= 0:
            hmmutils.buildHmmOverAlignment(outputAlignPath, hmmPath).run()
        else:
            hmmPath = addHmmPath
        hmmTasks = hmmutils.hmmAlignQueries(hmmPath, queriesPath)
        task.submitTasks(hmmTasks)
        for hmmTask in task.asCompleted(hmmTasks):
            hmmutils.mergeHmmAlignments([hmmTask.outputFile], rest_path, includeInsertions=False)
        external_tools.runEpaNg(outputAlignPath, bb_tree, rest_path, tempDir, outputJplacePath).run()
        external_tools.runGappaGraft(outputJplacePath, tempDir, outputTreePath).run()
    else:
        shutil.copy(bb_tree, outputTreePath)
    # now we have an alignment over, run epa-ng

def buildNaiveAlignment(sequences, tempDir, outputAlignPath):
    # assuming that buildInitialAlignment has already been called
    queriesPath = os.path.join(tempDir, "queries_all.txt")
    hmmDir = os.path.join(tempDir, "skeleton_hmm")
    sequenceutils.writeFasta(sequences, queriesPath)
    hmmPath = os.path.join(hmmDir, "hmm_model.txt")
    hmmTasks = hmmutils.hmmAlignQueries(hmmPath, queriesPath)
    task.submitTasks(hmmTasks)
    for hmmTask in task.asCompleted(hmmTasks):
        hmmutils.mergeHmmAlignments([hmmTask.outputFile], outputAlignPath, includeInsertions=True)