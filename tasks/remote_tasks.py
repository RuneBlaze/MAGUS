from tools import external_tools
import os
from huey import SqliteHuey

huey = SqliteHuey(filename="/home/baqiaol2/scratch/huey.db")

@huey.task()
def remote_mafft_linsi(input_path, output_path):
    t = external_tools.runMafft(input_path, None, os.path.dirname(output_path), output_path, 4)
    t.run()
    return t

@huey.task()
def remote_gcm(input_path, output_path):
    t = external_tools.runGcmC(input_path, None, os.path.dirname(output_path), output_path, 4)
    t.run()
    return t