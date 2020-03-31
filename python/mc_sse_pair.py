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
        self.wgt = np.zeros(16,dtype= np.float64)
        self.awgt = np.zeros(16,dtype= np.float64)
        self.dwgt =np.zeros(16,dtype= np.float64)
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
        self.cords_of_site = np.zeros((2,self.nn),dtype=np.int64) #given site number, what is the coordinates 
        self.site_of_cords = np.zeros((self.nx,self.ny),dtype=np.int64) #given coordinates, what is the site number 
        self.plqt = np.zeros((4,self.nq),dtype=np.int64)
        self.ns2iq = np.zeros((2,2,2,2),dtype=np.int64)
        self.iq2ns = np.zeros((4,16),dtype=np.int64)


        np.random.seed(int(time.time())) #set random seed when constructor called

    def initconf(self):
        for i in range(len(self.state)):
            self.state[i] = np.random.randint(2) #initializing initial state of lattice
            #putting 2 ensures only 0 and 1 are generated

    def lattice(self):
        i = 0
        for y in range(self.ny):
            for x in range(self.nx):
                self.cords_of_site[0,i] = x
                self.cords_of_site[1,i] = y
                self.site_of_cords[x,y] = i
                i = i+1

        for q in range(0,self.nn): 
            x1 = self.cords_of_site[0,q]; y1 = self.cords_of_site[1,q] #get cords of first site in plaquette i
            x2 = (x1+1) % self.nx ; y2 = y1 #imagine a plaquette
            x3 = x2; y3 = (y1+1) % self.ny
            x4 = x1; y4 = y3
            self.plqt[0,q] = self.site_of_cords[x1,y1]
            self.plqt[1,q] = self.site_of_cords[x2,y2]
            self.plqt[2,q] = self.site_of_cords[x3,y3]
            self.plqt[3,q] = self.site_of_cords[x4,y4]

        for iq in range(0,16):
            iiq = iq
            ns0 = iiq % 2; iiq = iiq // 2
            ns1 = iiq % 2; iiq = iiq // 2
            ns2 = iiq % 2; iiq = iiq // 2
            ns3 = iiq % 2
            self.ns2iq[ns0,ns1,ns2,ns3] = iq
            self.iq2ns[0,iq] = ns0
            self.iq2ns[1,iq] = ns1
            self.iq2ns[2,iq] = ns2
            self.iq2ns[3,iq] = ns3

    def pvect0(self):
        max_wgt = 0

        for iq in range(16):
            s1 = self.iq2ns[0,iq]
            s2 = self.iq2ns[1,iq]
            s3 = self.iq2ns[2,iq]
            s4 = self.iq2ns[3,iq]
            self.wgt[iq] = self.vv * float(s1*s2 + s2*s3 + s3*s4 + s4*s1) #nearest neighbor repulsion
            self.wgt[iq]  = self.wgt[iq] + self.vv2 * float(s1*s3 + s2*s4) #next nearest neighbor repulsion
            self.wgt[iq] = self.wgt[iq] - self.mu * float(s1+s2+s3+s4)/self.z
            
            if self.wgt[iq] > max_wgt:
                max_wgt = self.wgt[iq]

        max_wgt = max_wgt + 1
        for iq in range(16):
            self.awgt[iq] = max_wgt - self.wgt[iq]
            if self.awgt[iq] > 1e-6:
                self.dwgt[iq] = 1/self.awgt[iq]
            else:
                self.dwgt[iq] = 1e6 


