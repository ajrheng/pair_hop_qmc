import numpy as np

class mc_sse_class:

    def __init__(self,t,tp,vv,vv2,mu,nx):

        # Hamiltonian parameters
        self.t = t
        self.tp = tp
        self.vv = vv
        self.vv2 = vv2
        self.mu = mu

        # lattice parameters
        self.nx = nx
        self.ny = nx #square lattice
        self.nn = self.nx * self.ny
        self.nq = self.nn
        self.nb = 2 * self.nn
        self.z = 4 #coordination number is always 4 for a 2D square lattice

        #simulation variables
        self.num_vx = 62
        self.tot_vx = 2 ** 8 - 1
        self.ntau = 100
        self.wgt = np.zeros(16,dtype=np.float64)
        self.awgt = np.zeros(16,dtype=np.float64)
        self.dwgt =np.zeros(16,dtype=np.float64)
        self.vx_num_from_int = np.zeros(self.tot_vx+1,dtype=np.int64) #ivx, +1 so index tot_vex is not illegal
        self.vx_num_from_int[:] = -1 #-1 to indicate it is an invalid vertex
        self.int_from_vx_num = np.zeros(self.num_vx,dtype=np.int64) #vxi
        self.op = np.zeros((7,16),dtype=np.int64)
        self.oper_from_vx_num = np.zeros(self.num_vx,dtype=np.int64) #vxoper
        self.vx_num_aft_op = np.zeros((7,16),dtype=np.int64) #vxcode
        self.vx_leg = np.zeros((8,self.num_vx),dtype=np.int64) 
        self.vx_new = np.zeros((8,8,self.num_vx),dtype=np.int64)
        self.vx_prob = np.zeros((8,8,self.num_vx),dtype=np.float64)
        self.ctra = np.zeros((7,2,2,2,2),dtype=np.int64)
        self.l = 20 #initial length of opstring and related arrays
        self.state = np.zeros(self.nn,dtype=np.int8) #create 1D array
        self.opstring = np.zeros(self.l,dtype=np.int64) 
        self.vert = np.zeros(self.l-1,dtype=np.int64)
        self.link = np.zeros(8*self.l-1,dtype=np.int64)

        np.random.seed(int(time.time())) #set random seed when constructor called

    def initconf(self):
        for i in range(len(self.state)):
            self.state[i] = np.random.randint(2) #initializing initial state of lattice
            #putting 2 ensures only 0 and 1 are generated


