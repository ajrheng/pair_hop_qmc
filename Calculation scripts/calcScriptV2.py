import math
import os
import matplotlib.pyplot as plt
plt.rcParams['font.family']='serif'
plt.rcParams['mathtext.fontset']='stix'

def calcStg(inputList):
    zeropi = []
    pipi = []
    for each in inputList:
        os.chdir(each)
        with open('stgfull.dat','r') as datfile:
            counter=0
            zeropiSum=0
            pipiSum=0
            for item in datfile:
                if counter==8 or counter == 136:
                    zeropiSum+=float(item)
                elif counter==144:
                    pipiSum += float(item)
                counter+=1
            zeropi.append(zeropiSum)
            pipi.append(pipiSum)
        os.chdir('..')
    return zeropi,pipi
            

def mainCalculation(inputList,dattype):
    mean = []
    stddev=[]
    for each in inputList:
        os.chdir(each)
        with open(dattype,'r') as datfile:
            tempList=[]
            tempList1=[]
            for line in datfile:
                p=line.split()
                tempList.append(float(p[1]))
                tempList1.append(float(p[0]))
            if dattype!= 'uni.dat':
                tempList=[sum(x) for x in zip(tempList,tempList1)]
            mean.append(calcMean(tempList))
            stddev.append(calcStdDev(tempList))
    os.chdir('..')
    return mean,stddev


def listOfDirect():
    tempList=os.listdir()
    outputList=[]
    cwd=os.getcwd()
    xaxis=[]
    for item in tempList:
        if item[:2]=='mu':
            outputList.append(cwd+'/'+item)
            xaxis.append(float(item[3:]))
    outputList.sort(key=lambda x: int(''.join(filter(str.isdigit,x))))
    return outputList,sorted(xaxis)
    

def calcMean(inputList):
    length = len(inputList)
    sumOfItems = float(0)
    for item in inputList:
        sumOfItems+=float(item)
    sumOfItems/=length
    return sumOfItems

def calcStdDev(inputList):
    mean=calcMean(inputList)
    length=len(inputList)
    stdDevSum=float(0)
    for item in inputList:
        stdDevSum+=(mean-float(item))**2
    stdDevSum=math.sqrt(stdDevSum/(length*(length-1)))
    return stdDevSum


directList, xaxis = listOfDirect()
#==========Uni/rho plot================#
mean, stddev = mainCalculation(directList,'uni.dat')
densitymean=mean
fig, ax1 = plt.subplots()
ax1.set_xlabel('$\mu$')
ax1.set_ylabel('<N>',color='b')
ax1.errorbar(xaxis,mean,stddev,marker='o',label='<N>')
mean, stddev = mainCalculation(directList,'rho.dat')
ax2 = ax1.twinx()
ax2.set_ylabel('$\\rho_s$',color='m')
ax2.errorbar(xaxis,mean,stddev,color='tab:orange',marker='^',label='$\\rho_s$')
ax1.set_ylim(0,1)
ax2.set_ylim(0,1.6)
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc=0)
fig.tight_layout()
plt.savefig('unirho.png', dpi=300)

#=========rhot/rhotp plot================#
mean,stddev=mainCalculation(directList,'rhot.dat')
fig, ax1 = plt.subplots()
ax1.set_xlabel('$\mu$')
ax1.set_ylabel('$\\rho_t$')
ax1.errorbar(xaxis,mean,stddev,marker='.',label='$\\rho_t$')
mean, stddev = mainCalculation(directList,'rhotp.dat')
ax2 = ax1.twinx()
ax2.set_ylabel('$\\rho_{tp}$')
ax2.errorbar(xaxis,mean,stddev,color='tab:orange',marker='^',label='$\\rho_{tp}$')
ax2.set_ylim(0,1.5)
ax1.set_ylim(0,1.5)
lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc=0)
fig.tight_layout()
plt.savefig('rhotrhotp', dpi=300)

#========structure factor plot===========#
zeropi, pipi = calcStg(directList)
fig, ax1=plt.subplots()
ax1.plot(densitymean,pipi,marker='.')
ax1.plot(densitymean,zeropi,marker='^')
ax1.set_ylim([0,50])
ax1.legend(['$S(\pi,\pi)$','$S(0,\pi),(\pi,0)$'])
ax1.set_xlabel('<N>')
fig.tight_layout
plt.savefig('stgfull',dpi=300)


