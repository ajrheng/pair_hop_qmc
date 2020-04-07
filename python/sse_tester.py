from mc_sse_dimer import *

if __name__ == "__main__":
    test = mc_sse_dimer(1,0,10,4)
    test.init_state()
    test.lattice()
    test.init_matrix_ele()
    test.pvect0()
    test.vxweight()
    test.initvrtx()
    test.diagonal_update()
    test.linked_list()
    test.t_loop_update()
    test.d_loop_update()
    test.aft_loop_update()