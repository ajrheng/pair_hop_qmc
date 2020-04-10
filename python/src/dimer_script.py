import os
import numpy as np
from shutil import copy

j2 = 0
beta = 10
l = 10

cwd = os.getcwd()
for j in np.arange(0,2,0.2):
    folder_name = "J="+str(j)
    folder_path = os.path.join(cwd,folder_name)
    driver = os.path.join(cwd,'sse_driver.py')
    script = os.path.join(cwd,'script.sh')
    os.mkdir(folder_path)
    copy(driver,folder_path)
    copy(script,folder_path)
    os.chdir(folder_path)
    
    file = open('params.txt','a')
    file.write(str(j)+'\n')
    file.write(str(j2)+'\n')
    file.write(str(beta)+'\n')
    file.write(str(l))
    file.close()

    os.system('qsub script.sh')

    os.chdir(cwd)
