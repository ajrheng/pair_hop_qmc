import numpy as np
import time

#======================init params========================================
nx = 5; ny = 5; nn = nx * ny
nq = nn #num of plaquettes
nb = 2 * nn #num of bonds
z = 4 #coordination number

num_vx = 52 #num possible vertex
ivmax = 2 ** 8 - 1 #???
ntau = 100


l = 20 #initial length of opstring and related arrays
state = np.zeros(nn,dtype=np.int8) #create 1D array

np.random.seed(int(time.time()))

#===============initconf=====================================
for i in range(len(state)):
    state[i] = np.random.randint(2) #initializing initial state of lattice
    #putting 2 ensures only 0 and 1 are generated

opstring = np.zeros(l,dtype=np.int64) 
vert = np.zeros(l-1,dtype=np.int64)
link = np.zeros(8*l-1,dtype=np.int64)
#=============================================================

#============lattice==========================================
cords_of_site = np.zeros((2,nn),dtype=np.int64) #given site number, what is the coordinates 
site_of_cords = np.zeros((nx,ny),dtype=np.int64) #given coordinates, what is the site number 
plqt = np.zeros((4,nq),dtype=np.int64)
ns2iq = np.zeros((2,2,2,2),dtype=np.int64)
iq2ns = np.zeros((4,16),dtype=np.int64)

i = 0
for y in range(ny):
    for x in range(nx):
        cords_of_site[0,i] = x
        cords_of_site[1,i] = y
        site_of_cords[x,y] = i
        i = i+1

for q in range(0,nn): 
    x1 = cords_of_site[0,q]; y1 = cords_of_site[1,q] #get cords of first site in plaquette i
    x2 = (x1+1) % nx ; y2 = y1 #imagine a plaquette
    x3 = x2; y3 = (y1+1) % ny
    x4 = x1; y4 = y3
    plqt[0,q] = site_of_cords[x1,y1]
    plqt[1,q] = site_of_cords[x2,y2]
    plqt[2,q] = site_of_cords[x3,y3]
    plqt[3,q] = site_of_cords[x4,y4]

for iq in range(0,16):
    iiq = iq
    ns0 = iiq % 2; iiq = iiq // 2
    ns1 = iiq % 2; iiq = iiq // 2
    ns2 = iiq % 2; iiq = iiq // 2
    ns3 = iiq % 2
    ns2iq[ns0,ns1,ns2,ns3] = iq
    iq2ns[0,iq] = ns0
    iq2ns[1,iq] = ns1
    iq2ns[2,iq] = ns2
    iq2ns[3,iq] = ns3

#==============================================================