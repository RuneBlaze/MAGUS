'''
Created on Apr 14, 2020

@author: Vlad
'''

import subprocess
import os
import random
import shutil
from configuration import Configs
from tasks.task import Task
from helpers.translate_raxml import translate_raxml

def runCommand(**kwargs):
    command = kwargs["command"]
    Configs.log("Running an external tool, command: {}".format(command))
    runner = subprocess.run(command, shell = True, cwd = kwargs["workingDir"], universal_newlines = True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    try:    
        runner.check_returncode()
    except:
        Configs.error("Command encountered error: {}".format(command))
        Configs.error("Exit code: {}".format(runner.returncode))
        Configs.error("Output: {}".format(runner.stdout))
        raise
    for srcPath, destPath in kwargs.get("fileCopyMap", {}).items():
        if os.path.isfile(destPath):
            Configs.log("Removing existing file: {}".format(destPath))
            os.remove(destPath)
        shutil.move(srcPath, destPath)

def runClustalOmegaGuideTree(fastaPath, workingDir, outputPath, threads = 1):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    args = [Configs.clustalPath]
    args.extend(["-i", fastaPath, "--max-hmm-iterations=-1", "--guidetree-out={}".format(tempPath)])
    args.extend(["--threads={}".format(threads)])
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {tempPath : outputPath}, "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)

def generateMafftFilePathMap(inputPaths, outputDir):
    mafftMap = {inputPath : os.path.join(outputDir, "mafft_{}".format(os.path.basename(inputPath))) for inputPath in inputPaths}
    return mafftMap

def buildMafftAlignments(inputOutputPathMap):
    tasks = [buildMafftAlignment(inputPath, outputPath) for inputPath, outputPath in inputOutputPathMap.items()]
    return tasks
    
def buildMafftAlignment(inputPath, outputPath, subtablePath = None, useFafft = False):
    if not useFafft:
        return runMafft(inputPath, subtablePath, Configs.workingDir, outputPath, Configs.numCores)      
    else:
        return runFafft(inputPath, subtablePath, Configs.workingDir, outputPath, Configs.numCores)          

def runMafft(fastaPath, subtablePath, workingDir, outputPath, threads = 1):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    ep = "0.123"
    if Configs.newMafft:
        ep = "0"
    args = [Configs.mafftPath, "--localpair", "--maxiterate", "1000", "--ep", ep, 
            "--quiet", "--thread", str(threads), "--anysymbol"]
    if subtablePath is not None:
        args.extend(["--merge", subtablePath])
    args.extend([fastaPath, ">", tempPath])
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {tempPath : outputPath}, "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)

def runMyself(fastaPath, workingDir, outputPath, guideTree, threads = 1):
    # tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    args = {"workingDir" : workingDir, "outputFile" : outputPath,
            "subalignmentPaths" : Configs.subalignmentPaths, "sequencesPath" : fastaPath,
            "backbonePaths" : Configs.backbonePaths, "guideTree" : guideTree}
    task = createAlignmentTask(args)

def runFafft(fastaPath, subtablePath, workingDir, outputPath, threads = 1):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    if not Configs.targetLength:
        Configs.error("Target length not set")
        exit()
    assert not subtablePath
    args = ["python3", Configs.fafftPath]
    args.extend(["-i", fastaPath, "-w", workingDir, "-t", str(threads)])
    args.extend(["-o", tempPath, "-l", str(Configs.targetLength)])
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {tempPath : outputPath}, "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)

def runMafftAdd(existingAlnP, newSeqP, workingDir, outputPath, threads = 1, frag = False):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    args = [Configs.mafftPath]
    if frag:
        args.extend(["--auto", "--addfragments", newSeqP])
    else:
        args.extend(["--add", newSeqP])
    args.extend(["--keeplength", "--quiet","--anysymbol"])
    args.extend(["--thread", str(threads)])
    args.extend([existingAlnP, ">", tempPath])
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {tempPath : outputPath}, "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)

def splitAlignment(existingAlnP, combinedAlnP, workingDir, outputRefPath, outputQueryPath):
    args = ['epa-ng']
    args.extend(["--split", existingAlnP, combinedAlnP])
    args.extend(["-w", workingDir])
    refPath = os.path.join(workingDir, "reference.fasta")
    queryPath = os.path.join(workingDir, "query.fasta")
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {refPath : outputRefPath, queryPath: outputQueryPath}, "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputQueryPath, taskArgs = taskArgs)

def runMafftGuideTree(fastaPath, workingDir, outputPath, threads = 1):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    treeFile = os.path.join(os.path.dirname(fastaPath),  "{}.tree".format(os.path.basename(fastaPath)))
    args = [Configs.mafftPath, "--retree", "0", "--treeout", "--parttree",
            "--quiet", "--thread", str(threads), "--anysymbol"]
    args.extend(["--partsize", "1000"])
    args.extend([fastaPath, ">", tempPath])
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {treeFile : outputPath}, "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)

def runMcl(matrixPath, inflation, workingDir, outputPath):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    args = [Configs.mclPath, matrixPath, "--abc", "-o", tempPath]
    if inflation is not None:
        args.extend(["-I", str(inflation)])
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {tempPath : outputPath}, "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)

def runMlrMcl(matrixPath, granularity, balance, inflation, workingDir, outputPath):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    args = [Configs.mlrmclPath, matrixPath, "-o", tempPath]
    if granularity is not None:
        args.extend(["-c", str(granularity)])
    if balance is not None:
        args.extend(["-b", str(balance)])    
    if inflation is not None:
        args.extend(["-i", str(inflation)])
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {tempPath : outputPath}, "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)

def runFastTree(fastaFilePath, workingDir, outputPath, mode = "normal", intree = None):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    
    args = [Configs.fasttreePath]
    if Configs.inferDataType(fastaFilePath) == "protein":
        if Configs.emulatePasta and not Configs.exp:
            args.extend(["-wag", "-gamma"])
        else:
            args.extend(["-lg"])
    else:
        args.extend(["-nt", "-gtr"])
    
    if intree is not None:
        args.extend(["-intree", intree])
    
    # FIXME: this is a poor way to detect if we are using VeryFastTree
    if "Very" not in Configs.fasttreePath:
        if mode == "fast":
            args.extend(["-fastest", "-nosupport"]) 
        elif mode == "faster":
            args.extend(["-fastest", "-nosupport", "-mlnni", "4" ]) 
        elif mode == "noml":
            args.extend(["-fastest", "-nosupport", "-noml"])
    else:
        args.extend(["-nosupport", "-threads", str(Configs.numCores), '-double-precision'])
    
    args.extend([fastaFilePath, ">", tempPath])
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {tempPath : outputPath}, "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)

def runRaxmlNg(fastaFilePath, workingDir, outputPath, threads = Configs.numCores):
    # raxml-ng --msa prim.phy --model GTR+G --prefix T4 --threads 2 --seed 2 --tree pars{25},rand{25}
    baseName = os.path.basename(outputPath).replace(".","")
    raxmlFile = os.path.join(workingDir, "{}.raxml.bestTree".format(baseName))
    seed = random.randint(1, 1000000)
    args = [Configs.raxmlPath,
            "--msa", fastaFilePath,
            "--prefix", baseName,
            "--threads", str(threads),
            "--seed", str(seed)]
    
    if Configs.inferDataType(fastaFilePath) == "protein":
        args.extend(["--model", "LG+G"])
    else:
        args.extend(["--model", "GTR+G"])
        
    args.extend(["--tree", "pars{{{}}}".format(1)])
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {raxmlFile : outputPath}, "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)

def runRaxmlEvaluate(msaPath, workingDir, treePath, outputPath):
    baseName = os.path.basename(outputPath).replace(".","")
    raxmlFile = os.path.join(workingDir, "{}.raxml.bestTree".format(baseName))
    bestModelFile = os.path.join(workingDir, "{}.raxml.bestModel".format(baseName))
    seed = random.randint(1, 1000000)
    args = [Configs.raxmlPath,
            "--evaluate",
            "--msa", msaPath,
            "--prefix", baseName,
            "--tree", treePath,
            "--threads", "auto",
            "--seed", str(seed)]
    if Configs.inferDataType(msaPath) == "protein":
        args.extend(["--model", "LG+G"])
    else:
        args.extend(["--model", "GTR+G"])
    taskArgs = {"command" : subprocess.list2cmdline(args), 
        "fileCopyMap" : {raxmlFile : outputPath, bestModelFile: outputPath + ".model"}, 
        "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)

def runOldRmEvaluate(msaPath, workingDir, treePath, outputPath):
    # FIXME: this only handles nucleotide data
    args = ["raxmlHPC", "-f", "e", "-t" , treePath, "-m", "GTRGAMMA", 
    "-s", msaPath, "-n", "oldraxml", "w", workingDir]
    outTreeP = os.path.join(workingDir, "RAxML_result.oldraxml")
    infoP = os.path.join(workingDir, "RAxML_info.oldraxml")
    taskArgs = {"command" : subprocess.list2cmdline(args), 
        "fileCopyMap" : {outTreeP : outputPath, infoP: outputPath + ".info"}, 
        "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)

def runEpaNg(refMsaP, refTreeP, queryMsaP, workingDir, outputPath, threads = Configs.numCores):
    args = ['epa-ng']
    baseName = os.path.basename(outputPath).replace(".","")
    epaNgFile = os.path.join(workingDir, "epa_result.jplace")
    with open(refTreeP + ".model") as fh:
        modelArgs = fh.read().strip().split(",")[0]
    args.extend(["--ref-msa", refMsaP, "--tree", refTreeP, "--query", queryMsaP, "-w", workingDir])
    Configs.log("Model args for EPA-ng: {}".format(modelArgs))
    args.extend(["--model", modelArgs])
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {epaNgFile : outputPath}, "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)

def runPplacer(refMsaP, refTreeP, queryMsaP, workingDir, outputPath, threads = Configs.numCores):
    args = ['pplacer']
    baseName = os.path.basename(outputPath).replace(".","")
    pplacerFile = os.path.join(workingDir, "{}.jplace".format(baseName))
    translate_raxml(refTreeP + ".info", refTreeP + ".oldinfo")
    args.extend(["-r", refMsaP, "-t", refTreeP, "-s", refTreeP + ".oldinfo", queryMsaP])
    args.extend(["-j", str(threads), "-o", pplacerFile])
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {pplacerFile : outputPath}, "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)

def runGappaGraft(jplaceFile, workingDir, outputPath, fullyResolve = True):
    args = ['gappa', 'examine', 'graft']
    if fullyResolve:
        args.extend(['--fully-resolve'])
    args.extend(['--jplace-path', jplaceFile])
    args.extend(['--out-dir', workingDir])
    graftRes = os.path.join(workingDir, os.path.basename(jplaceFile).replace(".jplace", ".newick"))
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {graftRes : outputPath}, "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)

def runHmmBuild(alignmentPath, workingDir, outputPath):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    args = [Configs.hmmbuildPath,'--ere', '0.59', "--cpu", "1"]
    args.extend(["--symfrac", "0.0", "--informat", "afa", tempPath, alignmentPath])
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {tempPath : outputPath}, "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)

def runHmmAlign(hmmModelPath, fragPath, workingDir, outputPath):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    args = [Configs.hmmalignPath, "-o", tempPath]
    args.extend([hmmModelPath, fragPath])
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {tempPath : outputPath}, "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)

def runHmmSearch(hmmModelPath, fragPath, workingDir, outputPath):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    args = [Configs.hmmsearchPath,"--noali", "--cpu", "1", "-o", tempPath, "-E", "99999999", "--max"]
    args.extend([hmmModelPath, fragPath])
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {tempPath : outputPath}, "workingDir" : workingDir}
    return Task(taskType = "runCommand", outputFile = outputPath, taskArgs = taskArgs)