import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

X=np.linspace(0,2,101)
Y=np.linspace(0,2,101)
strFac=np.zeros((101,101))

i=0
j=0
input = open('stgfull.dat',"r")
for line in input:
    s=line.lstrip()
    strFac[i,j]=abs(float(s))
    j=j+1
    if (j==101):
        j=0
        i=i+1
        

plot=plt.contourf(X,Y,strFac,cmap=cm.YlGnBu)
plt.colorbar(plot)

plt.xlabel('kx/$\pi$')
plt.ylabel('ky/$\pi$')
plt.title('t=1, tp=1, v=0, v2=8, $\mu$=20. Plot of S(kx,ky).')
plt.savefig('funny.png', dpi=300)


    