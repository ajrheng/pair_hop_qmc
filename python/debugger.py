from dimer_debug_class import *

if __name__ == "__main__":

    f = open('log.txt','w') 
    f.close()
    m = open('error.txt','w')
    m.close()

    test = mc_sse_dimer(1.3,0,10,10,equil_steps = 1, mc_steps = 10000, num_runs = 1)

    test.init_state()
    test.lattice()
    test.init_matrix_ele()
    test.pvect0()
    test.vxweight()
    test.initvrtx()

    for i in range(test.mc_steps):
        print('\n STARTING MC STEP ', i)
        print(test.state)
        if (i+1)%1000 == 0:
            test.adjust_trun_cutoff()
        test.one_mc_step(i,0)