from sys import path
from os.path import join
from julia import Julia
jl = Julia()
p = join(path[0], 'MagusNight')
jl.eval(f'push!(LOAD_PATH, "{p}")')
from julia import MagusNight

class JuliaInterface:
    magus_night = None
    jl = None

    @staticmethod
    def init_interface():
        if not JuliaInterface.magus_night:
            jl = Julia()
            JuliaInterface.jl = jl
            p = join(path[0], 'MagusNight')
            jl.eval(f'push!(LOAD_PATH, "{p}")') 
            from julia import MagusNight
            JuliaInterface.magus_night = MagusNight

def night():
    JuliaInterface.init_interface()
    return JuliaInterface.magus_night