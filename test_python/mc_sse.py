import numpy as np
import time

#======================init params========================================
nx = 5; ny = 5; nn = nx * ny
nq = nn #num of plaquettes
nb = 2 * nn #num of bonds
z = 4 #coordination number

num_vx = 52 #num non-illegal vertices, aka nvx
tot_vx = 2 ** 8 - 1 #num of vertices, including illegal ones, aka ivmax
ntau = 100

wgt = np.zeros(16,dtype=np.float64)
awgt = np.zeros(16,dtype=np.float64)
dwgt =np.zeros(16,dtype=np.float64)
vx_num_from_int = np.zeros(tot_vx+1,dtype=np.int64) #ivx, +1 so index tot_vex is not illegal
vx_num_from_int[:] = -1 #-1 to indicate it is an invalid vertex
int_from_vx_num = np.zeros(num_vx,dtype=np.int64) #vxi
op = np.zeros((7,16),dtype=np.int64)
oper_from_vx_num = np.zeros(num_vx,dtype=np.int64) #vxoper
vx_num_aft_op = np.zeros((7,16),dtype=np.int64) #vxcode
vx_leg = np.zeros((8,num_vx),dtype=np.int64) 
vx_new = np.zeros((8,8,num_vx),dtype=np.int64)
vx_prob = np.zeros((8,8,num_vx),dtype=np.float64)
ctra = np.zeros((7,2,2,2,2),dtype=np.int64)


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

#==================pvect0======================================
max_wgt = 0

vv = 1; mu = 1; vv2=1; t = 1; tp =1

for iq in range(16):
    s1 = iq2ns[0,iq]
    s2 = iq2ns[1,iq]
    s3 = iq2ns[2,iq]
    s4 = iq2ns[3,iq]
    wgt[iq] = vv * float(s1*s2 + s2*s3 + s3*s4 + s4*s1) #nearest neighbor repulsion
    wgt[iq]  = wgt[iq] + vv2 * float(s1*s3 + s2*s4) #next nearest neighbor repulsion
    wgt[iq] = wgt[iq] - mu * float(s1+s2+s3+s4)/z
    
    if wgt[iq] > max_wgt:
        max_wgt = wgt[iq]

max_wgt = max_wgt + 1
for iq in range(16):
    awgt[iq] = max_wgt - wgt[iq]
    if awgt[iq] > 1e-6:
        dwgt[iq] = 1/awgt[iq]
    else:
        dwgt[iq] = 1e6

#================vxweight================================

ns = np.zeros(8,dtype=np.int64)
vx_count = -1 #you want the first vertex to be 0 indexed.
#in the fortran code it didnt matter because fortran does not enforce
#strict array indexing. so no array out of bounds error for the last index.

for iq in range(16):
    for i in range(4):
        ns[i] = iq2ns[i,iq]
        ns[i+4] = ns[i]
    
    vx_int = 0
    for k in range(8):
        vx_int = vx_int + ns[k] * (2**k)
    vx_count += 1
    vx_num_from_int[vx_int] = vx_count
    int_from_vx_num[vx_count] = vx_int
    op[0,iq] = iq
    oper_from_vx_num[vx_count] = 0
    vx_num_aft_op[0,iq] = vx_count
    for k in range(8):
        vx_leg[k,vx_count] = ns[k]

#single boson hopping vertices
for iq in range(16):
    ns[:] = 0
    for i in range(4):
        ns[i] = iq2ns[i,iq]
        ns[i+4] = ns[i]

    for i in range(4):
        ns_temp = np.copy(ns) #shallow copy
        j = (i+1) % 4
        if ns_temp[i] != ns_temp[j]:
            ns_temp[i+4] = 1-ns[i]
            ns_temp[j+4] = 1-ns[j]


            vx_int = 0
            for k in range(8):
                vx_int = vx_int + ns_temp[k] * (2**k)
            vx_count += 1
            vx_num_from_int[vx_int] = vx_count 
            int_from_vx_num[vx_count] = vx_int
            jq = 0
            for k in range(4):
                jq = jq + ns_temp[k+4] * (2**k)
            op[i+1,iq] = jq
            oper_from_vx_num[vx_count] = i+1
            vx_num_aft_op[i+1,iq] = vx_count
            for k in range(8):
                vx_leg[k,vx_count] = ns_temp[k]

#pair hopping vertices
for iq in range(16):
    ns[:] = 0
    for i in range(4):
        ns[i] = iq2ns[i,iq]
        ns[i+4] = ns[i]

    if ns[0] == 0 and ns[1] == 0 and ns[2] == 1 and ns[3] == 1:
        ns[4] = 1; ns[5] =1; ns[6] = 0; ns[7] = 0
        vx_int = 0
        for k in range(8):
            vx_int = vx_int + ns[k] * (2**k)
        vx_count += 1
        vx_num_from_int[vx_int] = vx_count 
        int_from_vx_num[vx_count] = vx_int
        jq = 0
        for k in range(4):
            jq = jq + ns[k+4] * (2**k)
        op[5,iq] = jq
        oper_from_vx_num[vx_count] = 5
        vx_num_aft_op[5,iq] = vx_count
        for k in range(8):
            vx_leg[k,vx_count] = ns[k]

    elif ns[0] == 1 and ns[1] == 1 and ns[2] == 0 and ns[3] == 0:
        ns[4] = 0; ns[5] =0; ns[6] = 1; ns[7] = 1
        vx_int = 0
        for k in range(8):
            vx_int = vx_int + ns[k] * (2**k)
        vx_count += 1
        vx_num_from_int[vx_int] = vx_count 
        int_from_vx_num[vx_count] = vx_int
        jq = 0
        for k in range(4):
            jq = jq + ns[k+4] * (2**k)
        op[5,iq] = jq
        oper_from_vx_num[vx_count] = 5
        vx_num_aft_op[5,iq] = vx_count
        for k in range(8):
            vx_leg[k,vx_count] = ns[k]

    elif ns[0] == 1 and ns[1] == 0 and ns[2] == 0 and ns[3] == 1:
        ns[4] = 0; ns[5] =1; ns[6] = 1; ns[7] = 0
        vx_int = 0
        for k in range(8):
            vx_int = vx_int + ns[k] * (2**k)
        vx_count += 1
        vx_num_from_int[vx_int] = vx_count 
        int_from_vx_num[vx_count] = vx_int
        jq = 0
        for k in range(4):
            jq = jq + ns[k+4] * (2**k)
        op[6,iq] = jq
        oper_from_vx_num[vx_count] = 6
        vx_num_aft_op[6,iq] = vx_count
        for k in range(8):
            vx_leg[k,vx_count] = ns[k]

    elif ns[0] == 0 and ns[1] == 1 and ns[2] == 1 and ns[3] == 0:
        ns[4] = 1; ns[5] =0; ns[6] = 0; ns[7] = 1
        vx_int = 0
        for k in range(8):
            vx_int = vx_int + ns[k] * (2**k)
        vx_count += 1
        vx_num_from_int[vx_int] = vx_count 
        int_from_vx_num[vx_count] = vx_int
        jq = 0
        for k in range(4):
            jq = jq + ns[k+4] * (2**k)
        op[6,iq] = jq
        oper_from_vx_num[vx_count] = 6
        vx_num_aft_op[6,iq] = vx_count
        for k in range(8):
            vx_leg[k,vx_count] = ns[k]

#====================================================

#======================initvrtx=====================

for i in range(num_vx): #in fortran it loops from 1-nvx inclusive, here its 0 - nvx-1
    ns[:] = 0
    vx_int1 = int_from_vx_num[i]
    for j in range(8):
        ns[j] = vx_int1 % 2; vx_int1 = vx_int1//2
    
    for in_leg in range(8):
        ns_temp1 = np.copy(ns)
        ns_temp1[in_leg] = 1-ns_temp1[in_leg]
        ns_temp2 = np.copy(ns_temp1)

        for out_leg in range(8):
            ns_temp1[out_leg] = 1-ns_temp1[out_leg]
            vx_int2 = 0
            for k in range(8):
                vx_int2 = vx_int2 + ns_temp1[k] * (2**k)
            if vx_num_from_int[vx_int2] != -1:
                vx_new[out_leg,in_leg,i] = vx_num_from_int[vx_int2]
                oper = oper_from_vx_num[vx_num_from_int[vx_int2]]

                if oper == 0:
                    vx_int3 = 0
                    for k in range(4):
                        vx_int3 = vx_int3 + ns_temp1[k] * (2**k)
                    vx_prob[out_leg, in_leg, i] = awgt[vx_int3]
                
                elif oper == 1 or oper == 2 or oper ==3 or oper == 4:
                    vx_prob[out_leg,in_leg,i] = t
                    
                    if oper ==1 and (ns_temp1[1] > ns_temp1[0]):
                        ctra[oper,ns_temp1[0],ns_temp1[1],ns_temp1[2],ns_temp1[3]] = -1
                    elif oper ==1 and (ns_temp1[1] < ns_temp1[0]):
                        ctra[oper,ns_temp1[0],ns_temp1[1],ns_temp1[2],ns_temp1[3]] = 1
                    elif oper==2 and (ns_temp1[2] > ns_temp1[1]):
                        ctra[oper,ns_temp1[0],ns_temp1[1],ns_temp1[2],ns_temp1[3]] = -1
                    elif oper==2 and (ns_temp1[2] < ns_temp1[1]):
                        ctra[oper,ns_temp1[0],ns_temp1[1],ns_temp1[2],ns_temp1[3]] = 1
                    elif oper==3 and (ns_temp1[3] > ns_temp1[2]):
                        ctra[oper,ns_temp1[0],ns_temp1[1],ns_temp1[2],ns_temp1[3]] = 1
                    elif oper==3 and (ns_temp1[3] < ns_temp1[2]):
                        ctra[oper,ns_temp1[0],ns_temp1[1],ns_temp1[2],ns_temp1[3]] = -1
                    elif oper==4 and (ns_temp1[3] > ns_temp1[0]):
                        ctra[oper,ns_temp1[0],ns_temp1[1],ns_temp1[2],ns_temp1[3]] = -1
                    else:
                        ctra[oper,ns_temp1[0],ns_temp1[1],ns_temp1[2],ns_temp1[3]] = 1
                else:
                    vx_prob[out_leg,in_leg,i] = tp
                    if oper == 5 and (ns_temp1[1] < ns_temp1[2]):
                        ctra[oper,ns_temp1[0],ns_temp1[1],ns_temp1[2],ns_temp1[3]] = -2
                    elif oper ==5 and (ns_temp1[1] > ns_temp1[2]) :
                        ctra[oper,ns_temp1[0],ns_temp1[1],ns_temp1[2],ns_temp1[3]] = 2
                    elif oper == 6 and ns_temp1[0] > ns_temp1[1]:
                        ctra[oper,ns_temp1[0],ns_temp1[1],ns_temp1[2],ns_temp1[3]] = 2
                    else:
                        ctra[oper,ns_temp1[0],ns_temp1[1],ns_temp1[2],ns_temp1[3]] = -2
            ns_temp1 = np.copy(ns_temp2)


for i in range(num_vx):
    for in_leg in range(8):
        for out_leg in range(1,8):
            vx_prob[out_leg,in_leg,i] = vx_prob[out_leg,in_leg,i]+ vx_prob[out_leg-1,in_leg,i]
            

for i in range(num_vx):
    for in_leg in range(8):
        for out_leg in range(8):
            vx_prob[out_leg,in_leg,i] = vx_prob[out_leg,in_leg,i]/vx_prob[7,in_leg,i]
            if vx_prob[out_leg,in_leg,i] < 1e-6:
                vx_prob[out_leg,in_leg,i] = -1

#======================================================================================