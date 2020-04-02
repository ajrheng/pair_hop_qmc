from mc_sse_dimer import *

if __name__ == "__main__":
    test = mc_sse_dimer(1,0,10,4)
    test.lattice()
    test.init_matrix_ele()
    test.pvect0()
    test.vxweight()
    test.initvrtx()