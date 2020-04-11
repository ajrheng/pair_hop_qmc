import numpy as np
import sys

#sys.path.append('/home/users/ntu/aheng004/dimer/python/src')
sys.path.append('/home/alvin/Desktop/RA/python/src')

from mc_sse_dimer import *

if __name__ == "__main__":

    param = np.loadtxt('params.txt')
    if len(param) != 4:
        sys.exit("insufficient parameters")

    j1 = float(param[0])
    j2 = float(param[1])
    beta = float(param[2])
    L = int(param[3])


    test = mc_sse_dimer(j1,j2,beta,L)
    test.main()
    