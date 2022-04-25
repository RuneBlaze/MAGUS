from sys import path
from os.path import join
from configuration import Configs
import os
try:
    os.environ["JULIA_NUM_THREADS"] = str(Configs.numCores)
    from julia import Julia
    jl = Julia()
    p = join(path[0], 'MagusNight')
    jl.eval(f'push!(LOAD_PATH, "{p}")')
    from julia import MagusNight
except ModuleNotFoundError:
    import imp
    MagusNight = None