import numpy as np
import random
import time

def get_num_from_ns(ns):
    """
    given a len 4 array of states return the 4-bit num representation
    """
    num = 0
    for i in range(4):
        num += ns[i]*(4**i)
    return num

class mc_sse_dimer:

    #constants in the class that won't change
    MAX_BOND_NUM = 4**2
    MAX_VX_NUM = 4**4
    OP_NUM = 6 # operator numbers labelled 0-5
    Z = 4 #coordination number for square lattice

    #simulation variables
    wgt = np.zeros(MAX_BOND_NUM,dtype= np.float64)
    awgt = np.zeros(MAX_BOND_NUM,dtype= np.float64)
    dwgt = np.zeros(MAX_BOND_NUM,dtype= np.float64)
    vx_num_from_int = np.zeros(MAX_VX_NUM, dtype=np.int64) #ivx, +1 so index tot_vex is not illegal
    vx_num_from_int[:] = -1 #-1 to indicate it is an invalid vertex
    int_from_vx_num = np.zeros(MAX_VX_NUM,dtype=np.int64) #vxi
    op = np.zeros((OP_NUM,MAX_VX_NUM,),dtype=np.int64)
    oper_from_vx_num = np.zeros(MAX_VX_NUM,dtype=np.int64) #vxoper
    vx_num_aft_op = np.zeros((OP_NUM,MAX_BOND_NUM),dtype=np.int64) #vxcode
    vx_leg = np.zeros((4,MAX_VX_NUM),dtype=np.int64) 
    vx_new = np.zeros((4,4,4,MAX_VX_NUM),dtype=np.int64)
    t_worm_prob = np.zeros((4,4,4,MAX_VX_NUM),dtype=np.float64)
    d_worm_prob = np.zeros((4,4,4,MAX_VX_NUM),dtype=np.float64)
    vx_matrix_ele = np.zeros(MAX_VX_NUM, dtype=np.float64)

    l = 20 #initial length of opstring and related arrays
    opstring = np.zeros(l,dtype=np.int64) 
    vert = np.zeros(l,dtype=np.int64)
    link = np.zeros(4*l,dtype=np.int64)
    num_op = 0
    t_loop_len = 50
    d_loop_len = 50
    num_opers_t = 0
    num_opers_d = 0
    num_loops_t = 0
    num_loops_d = 0
    ns_to_iq = np.zeros((4,4),dtype=np.int64)
    iq_to_ns = np.zeros((2,MAX_BOND_NUM),dtype=np.int64)

    #observables variables
    num_op_for_energy = 0
    max_wgt = 0

    tp_wgt = np.zeros(4,dtype=np.float64)
    tm_wgt = np.zeros(4,dtype=np.float64)
    tz_wgt = np.zeros(4,dtype=np.float64)
    dp_wgt = np.zeros(4,dtype=np.float64)
    dm_wgt = np.zeros(4,dtype=np.float64)
    dz_wgt = np.zeros(4,dtype=np.float64)
    t2_wgt = np.zeros(4,dtype=np.float64)

    act_tp = np.zeros(4,dtype=np.int64); act_tp[:] = -1
    act_tm = np.zeros(4,dtype=np.int64); act_tm[:] = -1
    act_tz = np.zeros(4,dtype=np.int64); act_tz[:] = -1
    act_dp = np.zeros(4,dtype=np.int64); act_dp[:] = -1
    act_dm = np.zeros(4,dtype=np.int64); act_dm[:] = -1
    act_dz = np.zeros(4,dtype=np.int64); act_dz[:] = -1
    act_t2 = np.zeros(4,dtype=np.int64); act_t2[:] = -1

    nvx = 0 #counter for num of vertices
    passed = False 

    def __init__(self,j1,j2,beta,L,equil_steps = 50000, mc_steps = 10000, num_runs = 10):

        # Hamiltonian parameters
        self.j1 = j1
        self.j2 = j2
        self.beta = beta

        #MC parameters
        self.equil_steps = equil_steps
        self.mc_steps = mc_steps
        self.num_runs = num_runs

        # lattice parameters
        self.nx = L
        self.ny = L #square lattice
        self.nn = self.nx * self.ny
        self.nb = 2 * self.nn

        # arrays that depend on lattice size, hence put into constructor
        self.state = np.zeros(self.nn,dtype=np.int8) #create 1D array
        self.cords_of_site = np.zeros((2,self.nn),dtype=np.int64) #given site number, what is the coordinates 
        self.site_of_cords = np.zeros((self.nx,self.ny),dtype=np.int64) #given coordinates, what is the site number 
        self.bond = np.zeros((2,self.nb+1),dtype=np.int64) #index bond from 1 to nb (inclusive)
        self.first = np.zeros(self.nn,dtype=np.int64)
        self.last = np.zeros(self.nn, dtype=np.int64)

        random.seed(int(time.time())) #set random seed when constructor called

        with open('log.txt','a') as file:
            file.write("Parameters for this run: J = {0}, J2 = {1}, beta = {2}, L = {3}\n"\
                .format(self.j1, self.j2, self.beta, self.nx)) 

    def init_state(self):
        for i in range(len(self.state)):
            self.state[i] = random.choice([0,1,2,3]) #initializing initial state of lattice

    def lattice(self):
        i = 0
        for y in range(self.ny):
            for x in range(self.nx):
                self.cords_of_site[0,i] = x
                self.cords_of_site[1,i] = y
                self.site_of_cords[x,y] = i
                i = i+1

        q = 0
        for y in range(self.ny): 
            for x in range(self.nx):
                q = q+1
                x1 = x; y1 = y 
                x2 = (x1+1) % self.nx ; y2 = y1 
                self.bond[0,q] = self.site_of_cords[x1,y1]
                self.bond[1,q] = self.site_of_cords[x2,y2]

        for x in range(self.nx): 
            for y in range(self.ny):
                q = q+1
                x1 = x; y1 = y 
                x2 = x1 ; y2 = (y1+1)% self.ny 
                self.bond[0,q] = self.site_of_cords[x1,y1]
                self.bond[1,q] = self.site_of_cords[x2,y2]

        for iq in range(self.MAX_BOND_NUM):
            iiq = iq
            ns0 = iiq % 4; iiq = iiq // 4
            ns1 = iiq % 4; iiq = iiq // 4
            self.ns_to_iq[ns0,ns1] = iq
            self.iq_to_ns[0,iq] = ns0
            self.iq_to_ns[1,iq] = ns1

    def init_matrix_ele(self):

        self.dz_wgt[0] = 1
        self.dp_wgt[0] = np.sqrt(2)
        self.dm_wgt[0] = np.sqrt(2)

        self.tz_wgt[1] = -1
        self.tp_wgt[1] = np.sqrt(2)
        self.dp_wgt[1] = np.sqrt(2)
        self.t2_wgt[1] = 2

        self.tp_wgt[2] = np.sqrt(2)
        self.tm_wgt[2] = np.sqrt(2)
        self.dz_wgt[2] = 1
        self.t2_wgt[2] = 2

        self.tz_wgt[3] = 1
        self.tm_wgt[3] = np.sqrt(2)
        self.dm_wgt[3] = np.sqrt(2)
        self.t2_wgt[3] = 2

        self.act_dz[0] = 2
        self.act_dp[0] = 3
        self.act_dm[0] = 1

        self.act_tz[1] = 1
        self.act_tp[1] = 2
        self.act_dp[1] = 0
        self.act_t2[1] = 1

        self.act_tp[2] = 3
        self.act_tm[2] = 1
        self.act_dz[2] = 0
        self.act_t2[2] = 2

        self.act_tz[3] = 3
        self.act_tm[3] = 2
        self.act_dm[3] = 0
        self.act_t2[3] = 3

    def pvect0(self):

        for iq in range(self.MAX_BOND_NUM):
            s1 = self.iq_to_ns[0,iq]
            s2 = self.iq_to_ns[1,iq]
            self.wgt[iq] = 0.5 * (self.j1 + self.j2) * self.tz_wgt[s1] * self.tz_wgt[s2] 
            self.wgt[iq] += ( (0.5 * self.t2_wgt[s1] - 3/4) + (0.5*self.t2_wgt[s2] - 3/4) )/self.Z    
            if self.wgt[iq] > self.max_wgt:
                self.max_wgt = self.wgt[iq]

        #max_wgt += 1
        # self.wgt = np.add(self.wgt,0.5*(self.j1+self.j2))
        # self.awgt[:] = self.wgt[:]

        for iq in range(self.MAX_BOND_NUM):
            self.awgt[iq] = self.max_wgt - self.wgt[iq]
            if self.awgt[iq] > 1e-6:
                self.dwgt[iq] = 1/self.awgt[iq]
            else:
                self.dwgt[iq] = 1e6 

    def vxweight(self):

        self.vx_num_from_int[:] = -1
        self.vx_leg[:,:] = -1
        self.vx_num_aft_op[:,:] = -1
        self.nvx = 0
        ns = np.zeros(4,dtype=np.int64)

        #diagonal vertices
        opnum = 0
        for iq in range(self.MAX_BOND_NUM):
            if self.awgt[iq] != 0:
                ns[0] = ns[2] = self.iq_to_ns[0,iq]
                ns[1] = ns[3] = self.iq_to_ns[1,iq]
                iiv = get_num_from_ns(ns)
                self.nvx += 1
                self.vx_num_from_int[iiv] = self.nvx
                self.int_from_vx_num[self.nvx] = iiv
                self.op[opnum,iq] = iq
                self.oper_from_vx_num[self.nvx] = opnum
                self.vx_num_aft_op[opnum,iq] = self.nvx
                for i in range(4):
                    self.vx_leg[i,self.nvx] = ns[i]
                self.vx_matrix_ele[iiv] = self.awgt[iq]

        for iq in range(self.MAX_BOND_NUM):
            ns[:] = 0
            ns[0] = self.iq_to_ns[0,iq]
            ns[1] = self.iq_to_ns[1,iq]

            if (self.act_tp[ns[0]] != -1 and self.act_tm[ns[1]]!= -1):
                #T+ T- 
                opnum = 1
                ns[2] = self.act_tp[ns[0]]
                ns[3] = self.act_tm[ns[1]]
                iiv = get_num_from_ns(ns)
                self.nvx += 1
                jq = 0
                for i in range(2):
                    jq += ns[i+2]*(4**i)
                self.vx_num_from_int[iiv] = self.nvx
                self.int_from_vx_num[self.nvx] = iiv
                self.op[opnum,iq] = jq
                self.oper_from_vx_num[self.nvx] = opnum
                self.vx_num_aft_op[opnum,iq] = self.nvx
                for i in range(4):
                    self.vx_leg[i,self.nvx] = ns[i]
                self.vx_matrix_ele[iiv] = 0.5*(self.j1+self.j2)

            if (self.act_tm[ns[0]] != -1 and self.act_tp[ns[1]]!= -1):
                #T- T+ 
                opnum = 2
                ns[2] = self.act_tm[ns[0]]
                ns[3] = self.act_tp[ns[1]]
                iiv = get_num_from_ns(ns)
                self.nvx += 1
                jq = 0
                for i in range(2):
                    jq += ns[i+2]*(4**i)
                self.vx_num_from_int[iiv] = self.nvx
                self.int_from_vx_num[self.nvx] = iiv
                self.op[opnum,iq] = jq
                self.oper_from_vx_num[self.nvx] = opnum
                self.vx_num_aft_op[opnum,iq] = self.nvx
                for i in range(4):
                    self.vx_leg[i,self.nvx] = ns[i]
                self.vx_matrix_ele[iiv] = 0.5*(self.j1+self.j2)

        for iq in range(self.MAX_BOND_NUM):
            ns[:] = 0
            ns[0] = self.iq_to_ns[0,iq]
            ns[1] = self.iq_to_ns[1,iq]

            if (self.act_tz[ns[0]] != -1 and self.act_dz[ns[1]]!= -1):
                #Tz Dz
                opnum = 3
                ns[2] = self.act_tz[ns[0]]
                ns[3] = self.act_dz[ns[1]]
                iiv = get_num_from_ns(ns)
                self.nvx += 1
                jq = 0
                for i in range(2):
                    jq += ns[i+2]*(4**i)
                self.vx_num_from_int[iiv] = self.nvx
                self.int_from_vx_num[self.nvx] = iiv
                self.op[opnum,iq] = jq
                self.oper_from_vx_num[self.nvx] = opnum
                self.vx_num_aft_op[opnum,iq] = self.nvx
                for i in range(4):
                    self.vx_leg[i,self.nvx] = ns[i]
                self.vx_matrix_ele[iiv] = 0.5*abs(self.j1-self.j2)

            if (self.act_tp[ns[0]] != -1 and self.act_dm[ns[1]]!= -1):
                #T+ D-
                opnum = 4
                ns[2] = self.act_tp[ns[0]]
                ns[3] = self.act_dm[ns[1]]
                iiv = get_num_from_ns(ns)
                self.nvx += 1
                jq = 0
                for i in range(2):
                    jq += ns[i+2]*(4**i)
                self.vx_num_from_int[iiv] = self.nvx
                self.int_from_vx_num[self.nvx] = iiv
                self.op[opnum,iq] = jq
                self.oper_from_vx_num[self.nvx] = opnum
                self.vx_num_aft_op[opnum,iq] = self.nvx
                for i in range(4):
                    self.vx_leg[i,self.nvx] = ns[i]
                self.vx_matrix_ele[iiv] = 0.5*abs(self.j1-self.j2)

            if (self.act_tm[ns[0]] != -1 and self.act_dp[ns[1]]!= -1):
                #T- D+
                opnum = 5
                ns[2] = self.act_tm[ns[0]]
                ns[3] = self.act_dp[ns[1]]
                iiv = get_num_from_ns(ns)
                self.nvx += 1
                jq = 0
                for i in range(2):
                    jq += ns[i+2]*(4**i)
                self.vx_num_from_int[iiv] = self.nvx
                self.int_from_vx_num[self.nvx] = iiv
                self.op[opnum,iq] = jq
                self.oper_from_vx_num[self.nvx] = opnum
                self.vx_num_aft_op[opnum,iq] = self.nvx
                for i in range(4):
                    self.vx_leg[i,self.nvx] = ns[i]
                self.vx_matrix_ele[iiv] = 0.5*abs(self.j1-self.j2)

    def initvrtx(self):
        ns = np.zeros(4,dtype=np.int64)
        self.t_worm_prob[:,:,:,:] = 0
        self.d_worm_prob[:,:,:,:] = 0
        self.vx_new[:,:,:,:] = 0

        for i in range(1,self.nvx+1):
            iq = self.int_from_vx_num[i]
            for k in range(4):
                ns[k] = iq % 4; iq = iq//4

            for ic in range(4):
                instate = ns[ic]
                #====================== START OF T+ T- WORM============================#
                if (self.act_tp[instate]!=-1):
                    ns1 = np.copy(ns)
                    ns1[ic] = self.act_tp[instate]

                    for oc in range(4):
                        ns2 = np.copy(ns1)
                        outstate= ns2[oc]

                        if (self.act_tp[outstate]!=-1):
                            ns2[oc] = self.act_tp[outstate]
                            iiv = get_num_from_ns(ns2)
                            if (self.vx_num_from_int[iiv]!=-1):
                                vertex_num = self.vx_num_from_int[iiv]
                                self.vx_new[ic,oc,ns1[ic],i] = vertex_num
                                self.t_worm_prob[ic,oc,ns1[ic],i] = self.vx_matrix_ele[iiv]

                        if (self.act_tm[outstate]!=-1):
                            ns2[oc] = self.act_tm[outstate]
                            iiv = get_num_from_ns(ns2)
                            if (self.vx_num_from_int[iiv]!=-1):
                                vertex_num = self.vx_num_from_int[iiv]
                                self.vx_new[ic,oc,ns1[ic],i] = vertex_num
                                self.t_worm_prob[ic,oc,ns1[ic],i] = self.vx_matrix_ele[iiv]                        

                if (self.act_tm[instate]!=-1):
                    ns1 = np.copy(ns)
                    ns1[ic] = self.act_tm[instate]

                    for oc in range(4):
                        ns2 = np.copy(ns1)
                        outstate= ns2[oc]

                        if (self.act_tp[outstate]!=-1):
                            ns2[oc] = self.act_tp[outstate]
                            iiv = get_num_from_ns(ns2)
                            if (self.vx_num_from_int[iiv]!=-1):
                                vertex_num = self.vx_num_from_int[iiv]
                                self.vx_new[ic,oc,ns1[ic],i] = vertex_num
                                self.t_worm_prob[ic,oc,ns1[ic],i] = self.vx_matrix_ele[iiv]

                        if (self.act_tm[outstate]!=-1):
                            ns2[oc] = self.act_tm[outstate]
                            iiv = get_num_from_ns(ns2)
                            if (self.vx_num_from_int[iiv]!=-1):
                                vertex_num = self.vx_num_from_int[iiv]
                                self.vx_new[ic,oc,ns1[ic],i] = vertex_num
                                self.t_worm_prob[ic,oc,ns1[ic],i] = self.vx_matrix_ele[iiv]

                #====================== START OF D+ D- Dz WORM============================#
                if (self.act_dp[instate]!=-1):
                    ns1 = np.copy(ns)
                    ns1[ic] = self.act_dp[instate]

                    for oc in range(4):
                        ns2 = np.copy(ns1)
                        outstate= ns2[oc]

                        if (self.act_dp[outstate]!=-1):
                            ns2[oc] = self.act_dp[outstate]
                            iiv = get_num_from_ns(ns2)
                            if (self.vx_num_from_int[iiv]!=-1):
                                vertex_num = self.vx_num_from_int[iiv]
                                self.vx_new[ic,oc,ns1[ic],i] = vertex_num
                                self.d_worm_prob[ic,oc,ns1[ic],i] = self.vx_matrix_ele[iiv]

                        if (self.act_dm[outstate]!=-1):
                            ns2[oc] = self.act_dm[outstate]
                            iiv = get_num_from_ns(ns2)
                            if (self.vx_num_from_int[iiv]!=-1):
                                vertex_num = self.vx_num_from_int[iiv]
                                self.vx_new[ic,oc,ns1[ic],i] = vertex_num
                                self.d_worm_prob[ic,oc,ns1[ic],i] = self.vx_matrix_ele[iiv]   

                        if (self.act_dz[outstate]!=-1):
                            ns2[oc] = self.act_dz[outstate]
                            iiv = get_num_from_ns(ns2)
                            if (self.vx_num_from_int[iiv]!=-1):
                                vertex_num = self.vx_num_from_int[iiv]
                                self.vx_new[ic,oc,ns1[ic],i] = vertex_num
                                self.d_worm_prob[ic,oc,ns1[ic],i] = self.vx_matrix_ele[iiv]

                if (self.act_dm[instate]!=-1):
                    ns1 = np.copy(ns)
                    ns1[ic] = self.act_dm[instate]

                    for oc in range(4):
                        ns2 = np.copy(ns1)
                        outstate= ns2[oc]

                        if (self.act_dp[outstate]!=-1):
                            ns2[oc] = self.act_dp[outstate]
                            iiv = get_num_from_ns(ns2)
                            if (self.vx_num_from_int[iiv]!=-1):
                                vertex_num = self.vx_num_from_int[iiv]
                                self.vx_new[ic,oc,ns1[ic],i] = vertex_num
                                self.d_worm_prob[ic,oc,ns1[ic],i] = self.vx_matrix_ele[iiv]

                        if (self.act_dm[outstate]!=-1):
                            ns2[oc] = self.act_dm[outstate]
                            iiv = get_num_from_ns(ns2)
                            if (self.vx_num_from_int[iiv]!=-1):
                                vertex_num = self.vx_num_from_int[iiv]
                                self.vx_new[ic,oc,ns1[ic],i] = vertex_num
                                self.d_worm_prob[ic,oc,ns1[ic],i] = self.vx_matrix_ele[iiv]   

                        if (self.act_dz[outstate]!=-1):
                            ns2[oc] = self.act_dz[outstate]
                            iiv = get_num_from_ns(ns2)
                            if (self.vx_num_from_int[iiv]!=-1):
                                vertex_num = self.vx_num_from_int[iiv]
                                self.vx_new[ic,oc,ns1[ic],i] = vertex_num
                                self.d_worm_prob[ic,oc,ns1[ic],i] = self.vx_matrix_ele[iiv]


                if (self.act_dz[instate]!=-1):
                    ns1 = np.copy(ns)
                    ns1[ic] = self.act_dz[instate]

                    for oc in range(4):
                        ns2 = np.copy(ns1)
                        outstate= ns2[oc]

                        if (self.act_dp[outstate]!=-1):
                            ns2[oc] = self.act_dp[outstate]
                            iiv = get_num_from_ns(ns2)
                            if (self.vx_num_from_int[iiv]!=-1):
                                vertex_num = self.vx_num_from_int[iiv]
                                self.vx_new[ic,oc,ns1[ic],i] = vertex_num
                                self.d_worm_prob[ic,oc,ns1[ic],i] = self.vx_matrix_ele[iiv]

                        if (self.act_dm[outstate]!=-1):
                            ns2[oc] = self.act_dm[outstate]
                            iiv = get_num_from_ns(ns2)
                            if (self.vx_num_from_int[iiv]!=-1):
                                vertex_num = self.vx_num_from_int[iiv]
                                self.vx_new[ic,oc,ns1[ic],i] = vertex_num
                                self.d_worm_prob[ic,oc,ns1[ic],i] = self.vx_matrix_ele[iiv]   

                        if (self.act_dz[outstate]!=-1):
                            ns2[oc] = self.act_dz[outstate]
                            iiv = get_num_from_ns(ns2)
                            if (self.vx_num_from_int[iiv]!=-1):
                                vertex_num = self.vx_num_from_int[iiv]
                                self.vx_new[ic,oc,ns1[ic],i] = vertex_num
                                self.d_worm_prob[ic,oc,ns1[ic],i] = self.vx_matrix_ele[iiv]

        for i in range(1,self.nvx+1):
            for ic in range(4):
                for instate in range(4):
                    for oc in range(1,4):
                        self.t_worm_prob[ic,oc,instate,i] = self.t_worm_prob[ic,oc,instate,i] + self.t_worm_prob[ic,oc-1,instate,i]
                        self.d_worm_prob[ic,oc,instate,i] = self.d_worm_prob[ic,oc,instate,i] + self.d_worm_prob[ic,oc-1,instate,i]

        for i in range(1,self.nvx+1):
            for ic in range(4):
                for instate in range(4):
                    for oc in range(4):
                        if (self.t_worm_prob[ic,3,instate,i]!= 0):
                            self.t_worm_prob[ic,oc,instate,i] /= self.t_worm_prob[ic,3,instate,i]
                        if (self.d_worm_prob[ic,3,instate,i]!= 0):
                            self.d_worm_prob[ic,oc,instate,i] /= self.d_worm_prob[ic,3,instate,i]
                        if self.t_worm_prob[ic,oc,instate,i] < 1e-6:
                            self.t_worm_prob[ic,oc,instate,i] = -1
                        if self.d_worm_prob[ic,oc,instate,i] < 1e-6:
                            self.d_worm_prob[ic,oc,instate,i] = -1

    def diagonal_update(self):
        
        for i in range(self.l):
            ii = self.opstring[i]
            if ii == 0:
                b = random.randrange(1,self.nb+1) #note the range, randrange is not inclusive of end
                ns0 = self.state[self.bond[0,b]]
                ns1 = self.state[self.bond[1,b]]
                iq = self.ns_to_iq[ns0,ns1]
                accept_prob = self.awgt[iq]*self.beta*self.nb/(self.l-self.num_op)
                if accept_prob >= 1 or random.random() < accept_prob:
                    self.opstring[i] = 6 * b
                    self.num_op += 1

            elif ii%6 == 0:
                b = ii//6
                ns0 = self.state[self.bond[0,b]]
                ns1 = self.state[self.bond[1,b]]
                iq = self.ns_to_iq[ns0,ns1]
                accept_prob = self.dwgt[iq]*(self.l - self.num_op + 1)/(self.beta*self.nb)
                if accept_prob >= 1 or random.random() < accept_prob:
                    self.opstring[i] = 0
                    self.num_op -= 1

            else:
                b = ii//6
                o = ii%6
                ns0 = self.state[self.bond[0,b]]
                ns1 = self.state[self.bond[1,b]]
                iq = self.ns_to_iq[ns0,ns1]
                jq = self.op[o,iq]
                self.state[self.bond[0,b]] = self.iq_to_ns[0,jq]
                self.state[self.bond[1,b]] = self.iq_to_ns[1,jq]

        self.num_op_for_energy += self.num_op

    def linked_list(self):

        i = 0; i0 = 0; i1 = 1
        self.first[:] = -1; self.last[:] = -1; self.link[:] = -1; self.vert[:] = -1

        for j in range(self.l):
            ii = self.opstring[j]
            if ii != 0:
                o = ii%6
                b = ii//6
                s0 = self.bond[0,b]; s1 = self.bond[1,b]
                ns0 = self.state[s0]; ns1 = self.state[s1]
                iq = self.ns_to_iq[ns0,ns1]
                self.vert[i] = self.vx_num_aft_op[o,iq]
                jq = self.op[o,iq]
                self.state[s0] = self.iq_to_ns[0,jq]
                self.state[s1] = self.iq_to_ns[1,jq]
                if self.vert[i] == -1:
                    with open('error.txt','a') as file:
                        file.write('\nerror here, vert is -1\n')
                        file.write('s0: {0}, s1: {1}, ns0: {2}, ns1: {3}\n'.format(s0,s1,ns0,ns1))
                        file.write('opstring: {4}, i: {0}, o: {1}, iq: {2}, b: {3}\n'.format(i,o,iq,b,ii))
                        file.write('state aft op: {0}, {1}\n\n'.format(self.state[s0],self.state[s1]))
                p0 = self.last[s0]
                p1 = self.last[s1]
                if p0 != -1:
                    self.link[p0] = i0
                    self.link[i0] = p0
                else:
                    self.first[s0] = i0
                if p1 != -1:
                    self.link[p1] = i1
                    self.link[i1] = p1
                else:
                    self.first[s1] = i1
                self.last[s0] = i0 + 2
                self.last[s1] = i1 + 2
                i += 1; i0 += 4; i1 += 4

        for s in range(self.nn):
            i = self.first[s]
            if i != -1:
                p0 = self.last[s]
                self.link[p0] = i
                self.link[i] = p0

        # print(self.opstring)

        # for i in range(len(self.link)):
        #     print("p: {0}, link[p]: {1}, vert num: {2}".format(i,self.link[i],self.vert[i//4]))

    def t_loop_update(self):

        ml = 50*self.l

        for _ in range(self.t_loop_len):
            nv = 0
            vert_num0 = -1
            init_state = -1

            while init_state == -1 or vert_num0 == -1:
                init_p = random.randrange(0,4*self.num_op)
                vp0 = init_p//4
                vert_num0 = self.vert[vp0]
                in_leg0 = init_p%4
                init_state = self.vx_leg[in_leg0, vert_num0]
                if (init_state == -1 or vert_num0 == -1):
                    print(init_p, vp0, vert_num0, in_leg0, init_state)
            #print('out of while')

            bef_init_p = self.link[init_p]
            vx = self.vert[bef_init_p//4]
            in_leg = bef_init_p%4
            bef_init_state = self.vx_leg[in_leg, vx]

            p1 = init_p
            self.passed = False
        
            for i in range(ml):
                vp = p1//4
                vx = self.vert[vp]
                in_leg = p1%4
                in_state = self.vx_leg[in_leg, vx]

                if i==0:
                    if init_state == 0:
                        break
                    elif init_state == 2:
                        instate_aft_flip = random.choice([self.act_tp[init_state], self.act_tm[init_state]])
                        # if random.random() <= 0.5:
                        #     instate_aft_flip = self.act_tp[init_state]
                        # else:
                        #     instate_aft_flip = self.act_tm[init_state]
                    elif init_state == 1:
                        instate_aft_flip = self.act_tp[init_state]
                    else:
                        instate_aft_flip = self.act_tm[init_state]
                else:
                    instate_aft_flip = outstate_aft_flip

                if p1 == init_p:
                    init_state = instate_aft_flip
                elif p1 == bef_init_p:
                    bef_init_state = instate_aft_flip
            
                r = random.random()
                for out_leg in range(4):
                    #print(self.t_worm_prob[in_leg, out_leg, instate_aft_flip, vx])
                    if r <= self.t_worm_prob[in_leg, out_leg, instate_aft_flip, vx]:
                        new_vx = self.vx_new[in_leg, out_leg, instate_aft_flip, vx]
                        outstate_aft_flip = self.vx_leg[out_leg, new_vx]
                        self.vert[vp] = new_vx
                        break
                    if out_leg == 3:
                        with open('error.txt','a') as file:
                            file.write("didn't find a suitable outleg for t worm!\n") #should have reached break statement previously
                
                p1 = 4*vp + out_leg
                nv += 1

                if p1 == init_p:
                    init_state = outstate_aft_flip
                elif p1 == bef_init_p:
                    bef_init_state = outstate_aft_flip
                
                if (p1 == init_p and init_state == bef_init_state) or \
                    (p1 == bef_init_p and init_state == bef_init_state):
                    self.passed = True
                    break
                p1 = self.link[p1]

            self.num_opers_t += nv
            if self.passed is False:
                #print('t pass false')
                return 

        self.num_loops_t += self.t_loop_len

    def d_loop_update(self):

        ml = 50*self.l

        for _ in range(self.d_loop_len):
            nv = 0
            vert_num0 = -1
            init_state = -1

            while vert_num0 == -1 or init_state == -1:
                init_p = random.randrange(0,4*self.num_op)
                vp0 = init_p//4
                vert_num0 = self.vert[vp0]
                in_leg0 = init_p%4
                init_state = self.vx_leg[in_leg0, vert_num0]

            bef_init_p = self.link[init_p]
            vx = self.vert[bef_init_p//4]
            in_leg = bef_init_p%4
            bef_init_state = self.vx_leg[in_leg, vx]

            p1 = init_p
            self.passed = False

            for i in range(ml):
                vp = p1//4
                vx = self.vert[vp]
                in_leg = p1%4
                in_state = self.vx_leg[in_leg, vx]

                if i==0:
                    if init_state == 0:
                        instate_aft_flip = random.choice([self.act_dz[init_state],self.act_dp[init_state],\
                            self.act_dm[init_state]])
                        # r = random.random()
                        # if r <= 1/3:
                        #     instate_aft_flip = self.act_dz[init_state]
                        # elif r <= 2/3:
                        #     instate_aft_flip = self.act_dp[init_state]
                        # else:
                        #     instate_aft_flip = self.act_dm[init_state]
                    elif init_state == 1:
                        instate_aft_flip = self.act_dp[init_state]
                    elif init_state == 2:
                        instate_aft_flip = self.act_dz[init_state]
                    else:
                        instate_aft_flip = self.act_dm[init_state]
                else:
                    instate_aft_flip = outstate_aft_flip

                if p1 == init_p:
                    init_state = instate_aft_flip
                elif p1 == bef_init_p:
                    bef_init_state = instate_aft_flip

                r = random.random()
                for out_leg in range(4):
                    if r <= self.d_worm_prob[in_leg, out_leg, instate_aft_flip, vx]:
                        new_vx = self.vx_new[in_leg, out_leg, instate_aft_flip, vx]
                        outstate_aft_flip = self.vx_leg[out_leg, new_vx]
                        self.vert[vp] = new_vx
                        break
                    if out_leg == 3:
                        with open('error.txt','a') as file:
                            file.write("didn't find a suitable outleg for d worm!\n") #should have reached break statement previously
                    
                p1 = 4*vp + out_leg
                nv += 1

                if p1 == init_p:
                    init_state = outstate_aft_flip
                elif p1 == bef_init_p:
                    bef_init_state = outstate_aft_flip
                
                if (p1 == init_p and init_state == bef_init_state) or \
                    (p1 == bef_init_p and init_state == bef_init_state):
                    self.passed = True
                    break
                p1 = self.link[p1]

            self.num_opers_d += nv
            if self.passed is False:
                #print('d pass false')
                return 

        self.num_loops_d += self.d_loop_len

    def update_opstring(self):

        j=0
        for i in range(self.l):
            if self.opstring[i] != 0:
                self.opstring[i] = 6*(self.opstring[i]//6) + self.oper_from_vx_num[self.vert[j]]
                j += 1
        
        for i in range(self.nn):
            if self.first[i] != -1:
                in_leg = self.first[i]%2
                vp = self.first[i]//4
                self.state[i] = self.vx_leg[in_leg,self.vert[vp]]
                # if self.state[i] != 0 and self.state[i] != 1 \
                #     and self.state[i] != 2 and self.state[i] != 3:
                if self.state[i] not in [0,1,2,3]:
                    with open('error.txt','a') as file:
                        file.write("wrong state!\n")
                
            else:
                self.state[i] = random.choice([0,1,2,3])

    def adjust_trun_cutoff(self):

        temp_opstring = np.copy(self.opstring)
        dl = int(self.l/10) + 2
        if self.num_op < self.l-dl/2:
            return
        old_l = self.l
        self.l = self.l + dl

        self.opstring = np.zeros(self.l,dtype=np.int64)
        self.opstring[:old_l] = temp_opstring[:]

        self.vert = np.zeros(self.l,dtype=np.int64)
        self.link = np.zeros(4*self.l,dtype=np.int64)

        #print("new trun cutoff {0}".format(self.l))

    def adjust_loop_len(self):

        try:
            avg_op_per_loop_t = self.num_loops_t/self.num_opers_t
            #nl = 1+int(2*self.l/avg_op_per_loop_t)
            nl = 1+int(self.l/avg_op_per_loop_t)
            self.t_loop_len = int((self.t_loop_len+nl)/2)
        except ZeroDivisionError:
            pass

        try:
            avg_op_per_loop_d = self.num_loops_d/self.num_opers_d
            #nl = 1+int(2*self.l/avg_op_per_loop_d)
            nl = 1+int(self.l/avg_op_per_loop_d)
            self.d_loop_len = int((self.d_loop_len+nl)/2)
        except ZeroDivisionError:
            pass

        with open('log.txt','a') as file:
            file.write("new t and d loop lens are {0} and {1}\n".format(self.t_loop_len,self.d_loop_len))

    def write_observables(self):

        energy = -self.num_op_for_energy/(self.mc_steps * self.beta)
        energy += self.max_wgt * self.nb #diagonal shift
        energy /= self.nn #energy per dimer

        file = open('energy.txt','a')
        file.write(str(energy)+'\n')
        file.close()

    def set_zero(self):

        self.num_opers_t = 0
        self.num_loops_t = 0
        self.num_opers_d = 0
        self.num_loops_d = 0
        self.num_op_for_energy = 0

    def one_mc_step(self):

        self.passed = False
        while self.passed is False:
            self.diagonal_update()
            self.linked_list()
            self.t_loop_update()
            self.d_loop_update()
            self.update_opstring()


    def equilibration(self):
        with open('log.txt','a') as file:
            file.write('Starting equilibration\n')

        for i in range(1,int(self.mc_steps)+1):
            #print("equilibration step ", i)
            self.one_mc_step()
            self.adjust_trun_cutoff()
            # if i%(self.mc_steps//20) == 0: #call it 20 times
            #     self.write_conf(0)
                #self.adjust_loop_len()
                #self.set_zero()

        with open('log.txt','a') as file:
            file.write('Completed equilibration. L = {0}, t loop len = \
                {1}, d loop len = {2}\n'.format(self.l,self.t_loop_len, self.d_loop_len))

    def write_conf(self,i):
    
        with open('conf.txt','a') as file:
            file.write("Config for run: "+ str(i+1)+'\n')
            for i in range(self.nn):
                file.write(str(self.state[i])+" ")
                if (i+1)%self.nx == 0:
                    file.write("\n")
            file.write("\n")

    def main_mc_runs(self):

        for i in range(self.num_runs):

            with open('log.txt','a') as file:
                file.write('Starting run '+ str(i)+'\n')
            self.set_zero()

            for j in range(self.mc_steps):
                self.one_mc_step()

            self.write_observables()

            with open('log.txt','a') as file:
                file.write('Finished run '+ str(i)+'\n')
            self.write_conf(i)
    
    def main(self):

        #set up
        self.init_state()
        self.lattice()
        self.init_matrix_ele()
        self.pvect0()
        self.vxweight()
        self.initvrtx()
    
        #equilibrate
        self.equilibration()

        #run
        self.main_mc_runs()

            
