#! /usr/bin/env python3
import traceback
import contextlib
import math
import csv
import pandas as pd
import numpy as np

def main(type,steps):
	dataframe =  pd.read_csv("sf_equil.dat",sep='\t',names=['step no','pipi','pi 0'])

	step_no = list(dataframe['step no'])

	if type is 'pipi':
		sf = list(dataframe['pipi'])
	elif type is 'pi 0':
		sf = list(dataframe['pi 0'])

	csv_list = []

	i=0
	while i < len(step_no):
		if i%steps == 0:
			csv_list.append([step_no[i],sf[i]])
		i+=1

	with open('sf_equil.csv','w') as csvFile:
		writer = csv.writer(csvFile)
		writer.writerows(csv_list)
		print ("done")

if __name__ == "__main__":
	sf_type = ['pipi', 'pi 0']
	while True:
		type = input('What sf do you wanna calculate (pipi or pi 0): ')
		if type in sf_type:
			while True:
				steps = input('In what step no intervals do you wanna calculate: ')
				if isinstance(steps,int):
					break
			break
	main(text,steps)
