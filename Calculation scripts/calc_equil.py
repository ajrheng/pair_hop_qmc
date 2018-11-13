#! /usr/bin/env python3
import traceback
import contextlib
import math
import csv
import pandas as pd
import numpy as np

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
		print(np.mean(sf[-60000:]))
	step_no = list(dataframe['step no'])

	outputfile = "".join([type,'.csv'])
	csv_list = []

	i=0
	while i < len(step_no):
		if i%steps == 0:
			csv_list.append([step_no[i],sf[i]])
		i+=1


	with open(outputfile,'w') as csvFile:
		writer = csv.writer(csvFile)
		writer.writerows(csv_list)
		print ("done with ",outputfile)

if __name__ == "__main__":
	accepted_type = ['pipi', 'pi 0','rhot','rhotp']
	while True:
		steps = input('In what step no intervals do you wanna calculate: ')
		steps = int(steps)
		if steps > 0:
			break
	for items in accepted_type:
		main(items,steps)
