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
    ns_to_iq = np.zeros((4,4),dtype=np.int64)
    iq_to_ns = np.zeros((2,MAX_BOND_NUM),dtype=np.int64)

    tp_wgt = np.zeros(4,dtype=np.float64)
    tm_wgt = np.zeros(4,dtype=np.float64)
    tz_wgt = np.zeros(4,dtype=np.float64)
    dp_wgt = np.zeros(4,dtype=np.float64)
    dm_wgt = np.zeros(4,dtype=np.float64)
    dz_wgt = np.zeros(4,dtype=np.float64)

    act_tp = np.zeros(4,dtype=np.int64); act_tdp[:] = -1
    act_tm = np.zeros(4,dtype=np.int64); act_tdm[:] = -1
    act_tz = np.zeros(4,dtype=np.int64); act_tdz[:] = -1
    act_dp = np.zeros(4,dtype=np.int64); act_ddp[:] = -1
    act_dm = np.zeros(4,dtype=np.int64); act_ddm[:] = -1
    act_dz = np.zeros(4,dtype=np.int64); act_ddz[:] = -1

    nvx = 0 #counter for num of vertices

    def __init__(self,j1,j2,beta,nx):

        # Hamiltonian parameters
        self.j1 = j1
        self.j2 = j2
        self.beta = beta

        # lattice parameters
        self.nx = nx
        self.ny = nx #square lattice
        self.nn = self.nx * self.ny
        self.nb = 2 * self.nn

        # arrays that depend on lattice size, hence put into constructor
        state = np.zeros(self.nn,dtype=np.int8) #create 1D array
        cords_of_site = np.zeros((2,self.nn),dtype=np.int64) #given site number, what is the coordinates 
        site_of_cords = np.zeros((self.nx,self.ny),dtype=np.int64) #given coordinates, what is the site number 
        bond = np.zeros((2,self.nb),dtype=np.int64)

        np.random.seed(int(time.time())) #set random seed when constructor called

    def init_state(self):
        for i in range(len(self.state)):
            self.state[i] = np.random.randint(4) #initializing initial state of lattice
            #putting 2 ensures between 0 and 3 are generated

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
                x1 = x; y1 = y 
                x2 = (x1+1) % self.nx ; y2 = y1 
                self.bond[0,q] = self.site_of_cords[x1,y1]
                self.bond[1,q] = self.site_of_cords[x2,y2]
                q = q+1

        for x in range(self.nx): 
            for y in range(self.ny):
                x1 = x; y1 = y 
                x2 = x1 ; y2 = (y1+1)% self.ny 
                self.bond[0,q] = self.site_of_cords[x1,y1]
                self.bond[1,q] = self.site_of_cords[x2,y2]
                q = q+1

        for iq in range(self.MAX_BOND_NUM):
            iiq = iq
            ns0 = iiq % 2; iiq = iiq // 2
            ns1 = iiq % 2; iiq = iiq // 2
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

       self.tp_wgt[2] = np.sqrt(2)
       self.tm_wgt[2] = np.sqrt(2)
       self.dz_wgt[2] = 1

       self.tz_wgt[3] = 1
       self.tm_wgt[3] = np.sqrt(2)
       self.dm_wgt[3] = np.sqrt(2)

       self.act_dz[0] = 2
       self.act_dp[0] = 3
       self.act_dm[0] = 1

       self.act_tz[1] = 1
       self.act_tp[1] = 2
       self.act_dp[1] = 0

       self.act_tp[2] = 3
       self.act_tm[2] = 1
       self.act_dz[2] = 0

       self.act_tz[3] = 3
       self.act_tm[3] = 2
       self.act_dm[3] = 0

    def pvect0(self):

        for iq in range(self.MAX_BOND_NUM):
            s1 = self.iq_to_ns[0,iq]
            s2 = self.iq_to_ns[1,iq]
            self.wgt[iq] = 0.5 * (self.j1 + self.j2) * self.tdz_wgt[s1] * self.tdz_wgt[s2]            
            # if self.wgt[iq] > max_wgt:
            #     max_wgt = self.wgt[iq]

        # max_wgt = max_wgt + 1
        wgt = np.add(wgt,0.5*(self.j1+self.j2))
        awgt[:] = wgt[:]

        for iq in range(self.MAX_BOND_NUM):
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
                self.vx_matrix_ele[iiv] = 0.5*(j1+j2)

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
                self.vx_matrix_ele[iiv] = 0.5*(j1+j2)

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
                self.vx_matrix_ele[iiv] = 0.5*abs(j1-j2)

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
                self.vx_matrix_ele[iiv] = 0.5*abs(j1-j2)

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
                self.vx_matrix_ele[iiv] = 0.5*abs(j1-j2)

    def initvrtx(self):
        ns = np.zeros(4,dtype=np.int64)
        self.t_worm_prob[:,:,:,:] = 0
        self.d_worm_prob[:,:,:,:] = 0
        self.vx_new[:,:,:,:] = 0


        for i in range(1,nvx+1):
            iq = self.int_from_vx_num[i]
            ns[0] = self.iq_to_ns[0,iq]
            ns[1] = self.iq_to_ns[1,iq]
            ns[2] = self.iq_to_ns[2,iq]
            ns[3] = self.iq_to_ns[3,iq]

        for ic in range(4):
            instate = ns[ic]
            #====================== START OF T+ T- WORM============================#
            ns1[:] = np.copy(ns)
            
