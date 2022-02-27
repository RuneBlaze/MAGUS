'''
Created on Apr 14, 2020

@author: Vlad
'''

import os
import time
from sys import platform

from helpers import sequenceutils
from shutil import which as witch

PREFER_SELF_INSTALLED = True

def which(path):
    if PREFER_SELF_INSTALLED:
        return witch(path)
    else:
        return None

def relative_to_me(path):
    if platform == "darwin" or platform == "win32":
        return None
    r = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
    assert os.path.isfile(r)
    return r

class Configs:
    
    workingDir = None
    sequencesPath = None
    subsetPaths = None
    subalignmentPaths = None
    backbonePaths = None
    guideTree = "fasttree"
    outputPath = None
    dataType = None

    targetLength = None
    
    decompositionMaxNumSubsets = 25
    decompositionMaxSubsetSize = 50
    decompositionStrategy = "pastastyle"
    decompositionSkeletonSize = 300
    #decompositionKmhIterations = 1
    
    graphBuildMethod = "mafft"
    graphBuildHmmExtend = False
    graphBuildRestrict = False
    graphTransformMethod = None
    graphBuildStrategy = "random"
    
    graphClusterMethod = "mcl" 
    graphTraceMethod = "minclusters"
    graphTraceOptimize = False

    # useJulia = False
    
    mafftRuns = 10
    mafftSize = 200
    mclInflationFactor = 4

    upgmaKeepOrder = True
    
    constrain = True
    onlyGuideTree = False
    recurse = True
    recurseGuideTree = "fasttree"
    recurseThreshold = 200
    emulatePasta = False
    
    clustalPath = relative_to_me("tools/clustal/clustalo")
    mafftPath = which("mafft") or relative_to_me("tools/mafft/mafft")
    mclPath = which("mcl") or relative_to_me("tools/mcl/bin/mcl")
    mlrmclPath = relative_to_me("tools/mlrmcl/mlrmcl")
    hmmalignPath = which("hmmalign") or relative_to_me("tools/hmmer/hmmalign")
    hmmbuildPath = which("hmmbuild") or relative_to_me("tools/hmmer/hmmbuild")
    hmmsearchPath = which("hmmsearch") or relative_to_me("tools/hmmer/hmmsearch")
    fasttreePath = which("FastTreeMP") or relative_to_me("tools/fasttree/FastTreeMP") or which("FastTree")
    raxmlPath = relative_to_me("tools/raxmlng/raxml-ng")
    
    logPath = None
    errorPath = None
    debugPath = None
    
    numCores = 1
    searchHeapLimit = 5000
    alignmentSizeLimit = 100
    allowLossyCompression = True

    outputInitialAlignment = None
    onlyInitialAln = False
    
    @staticmethod
    def log(msg, path = None):
        print(msg)
        path = Configs.logPath if path is None else path
        Configs.writeMsg(msg, path)
    
    @staticmethod
    def error(msg, path = None):
        Configs.log(msg)
        path = Configs.errorPath if path is None else path
        Configs.writeMsg(msg, path)
    
    @staticmethod
    def debug(msg, path = None):
        path = Configs.debugPath if path is None else path
        Configs.writeMsg(msg, path)
    
    @staticmethod
    def writeMsg(msg, path):
        if path is not None:
            with open(path, 'a') as logFile:
                logFile.write("{}    {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S"), msg))
    
    @staticmethod
    def inferDataType(sequencesFile):
        if Configs.dataType is None:
            Configs.dataType = sequenceutils.inferDataType(sequencesFile)
            Configs.log("Data type wasn't specified. Inferred data type {} from {}".format(Configs.dataType.upper(), sequencesFile))
        return Configs.dataType 

def buildConfigs(args):
    Configs.outputPath = os.path.abspath(args.output)
    
    if args.directory is not None:
        Configs.workingDir = os.path.abspath(args.directory) 
    else:
        Configs.workingDir = os.path.join(os.path.dirname(Configs.outputPath), "magus_working_dir")
    if not os.path.exists(Configs.workingDir):
        os.makedirs(Configs.workingDir)
    
    Configs.sequencesPath = os.path.abspath(args.sequences) if args.sequences is not None else Configs.sequencesPath
    
    Configs.guideTree = os.path.abspath(args.guidetree) if args.guidetree is not None else Configs.guideTree
    if args.guidetree is not None:
        Configs.guideTree = os.path.abspath(args.guidetree) if os.path.exists(os.path.abspath(args.guidetree)) else args.guidetree
    
    Configs.subalignmentPaths = []
    for p in args.subalignments:
        path = os.path.abspath(p)
        if os.path.isdir(path):
            for filename in os.listdir(path):
                Configs.subalignmentPaths.append(os.path.join(path, filename))
        else:
            with open(path) as fh:
                l = fh.readline()
            if ">" not in l and "#" not in l:
                Configs.subalignmentPaths += l.split()
            else:
                Configs.subalignmentPaths.append(path)
    
    Configs.backbonePaths = []
    for p in args.backbones:
        path = os.path.abspath(p)
        if os.path.isdir(path):
            for filename in os.listdir(path):
                Configs.backbonePaths.append(os.path.join(path, filename))
        else:
            Configs.backbonePaths.append(path)

    if args.numprocs > 0:
        if os.cpu_count() < Configs.numCores:
            raise RuntimeError("Number of cores specified ({}) is greater than the number of cores available ({})".format(Configs.numCores, os.cpu_count()))
        else:
            Configs.numCores = args.numprocs
    else:
        Configs.numCores = os.cpu_count()

    Configs.decompositionMaxSubsetSize = args.maxsubsetsize
    Configs.decompositionMaxNumSubsets = args.maxnumsubsets
    Configs.decompositionStrategy = args.decompstrategy
    Configs.decompositionSkeletonSize = args.decompskeletonsize or 300
    Configs.dataType = args.datatype
    
    Configs.graphBuildMethod = args.graphbuildmethod
    Configs.graphBuildHmmExtend = args.graphbuildhmmextend.lower() == "true"
    # print(Configs.graphBuildHmmExtend)
    Configs.graphBuildRestrict = args.graphbuildrestrict.lower() == "true"
    Configs.graphBuildStrategy = args.graphbuildstrategy
    Configs.graphClusterMethod = args.graphclustermethod
    Configs.graphTraceMethod = args.graphtracemethod
    Configs.graphTraceOptimize = args.graphtraceoptimize.lower() == "true"

    Configs.mafftRuns = args.mafftruns
    Configs.mafftSize = args.mafftsize
    Configs.mclInflationFactor = args.inflationfactor
    Configs.mclReweightIts = args.reweightits
    
    Configs.constrain = args.constrain.lower() == "true"
    Configs.onlyGuideTree = args.onlyguidetree.lower() == "true"
    Configs.recurse = args.recurse.lower() == "true"
    Configs.recurseGuideTree = args.recurseguidetree
    Configs.recurseThreshold = args.recursethreshold
    
    Configs.logPath = os.path.join(Configs.workingDir, "log.txt")    
    Configs.errorPath = os.path.join(Configs.workingDir, "log_errors.txt")
    Configs.debugPath = os.path.join(Configs.workingDir, "log_debug.txt")

    # Configs.useJulia = args.julia
    Configs.upgmaKeepOrder = args.keepOrder
    Configs.upgmaZeroWeight = args.zeroWeight
    Configs.exp = args.exp
    # Configs.randomSamples = args.randomSamples
    Configs.skeletonSeqs = args.skeleton
    # Configs.upgmaNoNormalize = args.noNormalize
    Configs.graphTransformMethod = args.transform

    Configs.emulatePasta = args.p0

    if Configs.emulatePasta:
        Configs.log("Emulating PASTA: setting up PASTA-specific configs")
        Configs.decompositionSkeletonSize = args.decompskeletonsize or 100
        Configs.skeletonSeqs = "random"
    
    Configs.alignmentSizeLimit = args.alignsizelimit
    Configs.allowLossyCompression = args.allowlossycompression.lower() == "true"
    Configs.outputInitialAlignment = args.initialAln
    Configs.onlyInitialAln = args.onlyInitialAln