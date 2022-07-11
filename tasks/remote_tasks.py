from tools import external_tools
import os
from huey import RedisHuey
huey = RedisHuey(host=open('/projects/tallis/baqiaol2/multi/hostname.txt').read().strip(), port=6379, password='MAGUS')

@huey.task(retries=42, retry_delay=10)
def remote_mafft_linsi(input_path, output_path):
    t = external_tools.runMafft(input_path, None, os.path.dirname(output_path), output_path, 4)
    t.run()
    return t

@huey.task(retries=42, retry_delay=10)
def remote_gcm(input_path, support_size, output_path):
    t = external_tools.runGcmC(input_path, support_size, os.path.dirname(output_path), output_path, 4)
    t.run()
    return t