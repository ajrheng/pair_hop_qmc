import numpy as np
import time

#======================init params========================================
nx = 5
ny = 5
nn = nx * ny
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
cords_of_site = {} #given site number (key), what is the coordinates (value)
site_of_cords = {} #given coordinates (value), what is the site number (key)

i = 0
for y in range(ny):
    for x in range(nx):
        i = i+1
        cords_of_site[i] = (x,y)
        site_of_cords[(x,y)] = i

for q in range(1,nn+1): #+1 because end of range func is exclusive