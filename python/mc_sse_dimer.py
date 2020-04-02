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

    act_tp = np.zeros(4,dtype=np.int64); act_tp[:] = -1
    act_tm = np.zeros(4,dtype=np.int64); act_tm[:] = -1
    act_tz = np.zeros(4,dtype=np.int64); act_tz[:] = -1
    act_dp = np.zeros(4,dtype=np.int64); act_dp[:] = -1
    act_dm = np.zeros(4,dtype=np.int64); act_dm[:] = -1
    act_dz = np.zeros(4,dtype=np.int64); act_dz[:] = -1

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
        self.state = np.zeros(self.nn,dtype=np.int8) #create 1D array
        self.cords_of_site = np.zeros((2,self.nn),dtype=np.int64) #given site number, what is the coordinates 
        self.site_of_cords = np.zeros((self.nx,self.ny),dtype=np.int64) #given coordinates, what is the site number 
        self.bond = np.zeros((2,self.nb),dtype=np.int64)

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
            self.wgt[iq] = 0.5 * (self.j1 + self.j2) * self.tz_wgt[s1] * self.tz_wgt[s2]         
            # if self.wgt[iq] > max_wgt:
            #     max_wgt = self.wgt[iq]

        # max_wgt = max_wgt + 1
        self.wgt = np.add(self.wgt,0.5*(self.j1+self.j2))
        self.awgt[:] = self.wgt[:]

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


        counter = 0
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