import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib

matplotlib.rcParams['font.sans-serif']="Arial"
matplotlib.rcParams['font.family'] = "sans-serif"
plt.rcParams['mathtext.fontset'] = 'stix'

X=np.linspace(0,2,9)
Y=np.linspace(0,2,9)
strFac=np.zeros((9,9))

i=0
j=0
input = open('stgfull.dat',"r")
for line in input:
	s=line.lstrip()
	strFac[j,i]=abs(float(s))
	a=9*i+j
	if (a==0 or a==8 or a==72 or a==80):
		strFac[j,i]=0
	j=j+1
	if (j==9):
		j=0
		i=i+1

plot=plt.contourf(X,Y,strFac,cmap=cm.inferno)
plt.colorbar(plot)

plt.xlabel('$k_1/\pi$')
plt.ylabel('$k_2/\pi$')
plt.tight_layout()
#plt.title('t=1, tp=4, v=4, v2=2, $\mu$=20. Plot of S(kx,ky).')
plt.savefig('full_str_plot.pdf', dpi=300)
