#! /usr/bin/env python3
import traceback
import contextlib
import math
import csv
import pandas as pd
import numpy as np
from scipy.stats import binned_statistic

def return_binned(x,y,no_of_bins,log=True,type='mean'):
    """
    function takes two lists and returns the binned x and y
    """
    if log is True:
        #x_max=math.ceil(np.log10(max(x)))
        x_max=np.log10(np.max(x))
        #x_min=math.floor(np.log10(min(x)))
        x_min = np.min(x)
        if x_min != 0:
            x_min=np.log10(np.min(x))
        else:
            x_min=0
        s, edges, _ = binned_statistic(x,y,statistic=type,bins=np.logspace(x_min,x_max,no_of_bins))
    else:
        s, edges, _ = binned_statistic(x,y,statistic=type,bins=no_of_bins)

    #x=[]
    #for left_edge,right_edge in zip(edges,edges[1:]):
    #    x.append(np.mean([left_edge,right_edge]))
    #return x,s
    return edges[:-1]+np.diff(edges)/2, s, edges

def main(type,steps):
	if type == "pipi":
		dataframe =  pd.read_csv("sf_equil.dat",sep='\s+',names=['step no','pipi','pi 0'])
		sf = list(dataframe['pipi'])
	elif type == "pi 0":
		dataframe =  pd.read_csv("sf_equil.dat",sep='\s+',names=['step no','pipi','pi 0'])
		sf = list(dataframe['pi 0'])
	elif type == "rhot":
		dataframe =  pd.read_csv("rho_equil.dat",sep='\s+',names=['step no','rhot','rhotp'])
		sf = list(dataframe['rhot'])
	elif type == "rhotp":
		dataframe =  pd.read_csv("rho_equil.dat",sep='\s+',names=['step no','rhot','rhotp'])
		sf = list(dataframe['rhotp'])
	step_no = list(dataframe['step no'])

	outputfile = "".join([type,'.csv'])
	csv_list = []

	#i=0
	#while i < len(step_no):
	#	if i%steps == 0:
	#		csv_list.append([step_no[i],sf[i]])
	#	i+=1

	step_no, sf, _ = return_binned(step_no,sf,steps)

	for i in range(len(step_no)):
		csv_list.append([step_no[i],sf[i]])

	with open(outputfile,'w') as csvFile:
		writer = csv.writer(csvFile)
		writer.writerows(csv_list)
		print ("done with ",outputfile)

if __name__ == "__main__":
	accepted_type = ['pipi', 'pi 0','rhot','rhotp']
	while True:
		steps = input('How many logged bins do you want: ')
		steps = int(steps)
		if steps > 0:
			break
	for items in accepted_type:
		main(items,steps)
