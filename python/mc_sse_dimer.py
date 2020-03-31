import numpy as np
import random
import time

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

    tdp_wgt = np.zeros(4,dtype=np.float64)
    tdm_wgt = np.zeros(4,dtype=np.float64)
    tdz_wgt = np.zeros(4,dtype=np.float64)
    ddp_wgt = np.zeros(4,dtype=np.float64)
    ddm_wgt = np.zeros(4,dtype=np.float64)
    ddz_wgt = np.zeros(4,dtype=np.float64)

    act_tdp = np.zeros(4,dtype=np.int64); act_tdp[:] = -1
    act_tdm = np.zeros(4,dtype=np.int64); act_tdm[:] = -1
    act_tdz = np.zeros(4,dtype=np.int64); act_tdz[:] = -1
    act_ddp = np.zeros(4,dtype=np.int64); act_ddp[:] = -1
    act_ddm = np.zeros(4,dtype=np.int64); act_ddm[:] = -1
    act_ddz = np.zeros(4,dtype=np.int64); act_ddz[:] = -1

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
       
       self.ddz_wgt[0] = 1
       self.ddp_wgt[0] = np.sqrt(2)
       self.ddm_wgt[0] = np.sqrt(2)

       self.tdz_wgt[1] = -1
       self.tdp_wgt[1] = np.sqrt(2)
       self.ddp_wgt[1] = np.sqrt(2)

       self.tdp_wgt[2] = np.sqrt(2)
       self.tdm_wgt[2] = np.sqrt(2)
       self.ddz_wgt[2] = 1

       self.tdz_wgt[3] = 1
       self.tdm_wgt[3] = np.sqrt(2)
       self.ddm_wgt[3] = np.sqrt(2)

       self.act_ddz[0] = 2
       self.act_ddp[0] = 3
       self.act_ddm[0] = 1

       self.act_tdz[1] = 1
       self.act_tdp[1] = 2
       self.act_ddp[1] = 0

       self.act_tdp[2] = 3
       self.act_tdm[2] = 1
       self.act_ddz[2] = 0

       self.act_tdz[3] = 3
       self.act_tdm[3] = 2
       self.act_ddm[3] = 0


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


