#! /usr/bin/env python3
import os, os.path
import traceback
import contextlib
import math
import csv
import numpy as np
from scipy import stats

class getItems():
	def __init__(self):
		self.possible_dattype_x = ['rho_x','rho_tx','rho_tpx']
		self.possible_dattype_y = ['rho_y','rho_ty','rho_tpy']

	def get_individual_type(self,pwd,each,dattype):
		if dattype == 'rho_x' or dattype == 'rho_y':
			dattype_ = 'rho'
		elif dattype == 'rho_tx' or dattype == 'rho_ty':
			dattype_ = 'rhot'
		else:
			dattype_ = 'rhotp'
		#print(dattype_)
		individual_dir = "".join([pwd,'/',each,'/',dattype_,'.dat'])

		return individual_dir

	def get_dat(self,dattype):
		items = os.listdir(".")
		items_removed = []
		for x in items:
			# because there are some hidden folders that will be grabbed
			# as well as this script,
			# we check the first 2 letters of all directories grabbed if
			# it starts with mu and construct a new list to put them into.
			if x[:2] == "mu":
				m=float(x[3:]) #get mu number
				if m%1==0:#if it is whole number
					items_removed.append(int(m))
				else:#if mu has decimal eg. 13.4
					items_removed.append(m)


		mu_list = sorted(items_removed)
		for index,truncated in enumerate(mu_list):
			mu_list[index] = "mu="+str(truncated)

		#gets a list of all directories
		pwd = os.getcwd() #your current working drive
		return_list = []
		for each in mu_list:
			# /Users/Alvin/Desktop/cy2001 results/varymuwhenv=5/mu=30/uni.dat
			#each = subfolder.
			#{ subfolder | mean |  sd }
			individual_dir = self.get_individual_type(pwd,each,dattype)
			print(individual_dir)
			dict_to_append = {}

			with open (individual_dir,"r") as datfile:
				data_list = []
				datfilelist = datfile.read().split()
				counter = 0
				for item in datfilelist:
					if dattype in self.possible_dattype_x:
						if counter %2 == 0 : #if its left list meaning %2 == 0
							data_list.append(float(item))
					else:
						if counter %2 == 1 : #if its right list meaning %2 == 1
							data_list.append(float(item))
					counter+=1

				mean = self.calculate_mean(data_list)
				sd = self.standard_deviation(data_list)
				dict_to_append["subfolder"] = float(each[3:])
				dict_to_append["mean"] = mean
				dict_to_append["sd"] = sd

				return_list.append(dict_to_append)

		self.construct_csv(return_list,pwd,dattype)


	def calculate_mean(self,array_of_items):
		return np.mean(array_of_items)

	def standard_deviation(self,array_of_items):
		return stats.sem(array_of_items)

	def get_csv_headers(self,dattype):
		csv_List = [["mu",dattype+" mean","STDEV"]]
		#if dattype == "uni":
		#	csv_List = [["mu","Mean","STDEV"]]
		#elif dattype == "rho":
		#	csv_List = [["mu","rho mean","STDEV"]]
		#elif dattype == "rhot":
		#	csv_List = [["mu","rhot mean","STDEV"]]
		#elif dattype =="rhotp":
	#		csv_List = [["mu","rhotp mean","STDEV"]]
#		elif dattype == "stgpipi":
#			csv_List = [["mu","ssa mean","STDEV"]]
#		elif dattype == "curr_single":
#			csv_List = [["mu","single current","STDEV"]]
#		elif dattype == "curr_pair":
#			csv_List = [["mu","pair current","STDEV"]]

		return csv_List

	def get_file_name(self,pwd,dattype):
		csvfile = "".join([pwd,'/',dattype,'.csv'])

		return csvfile

	def construct_csv(self,array_of_calculations,pwd,dattype):
		csv_List = self.get_csv_headers(dattype)

		for row in array_of_calculations:
			csv_List.append([row['subfolder'],row['mean'],row['sd']])

		csvfilename = self.get_file_name(pwd,dattype)

		with open(csvfilename,'w') as csvFile:
			writer = csv.writer(csvFile)
			writer.writerows(csv_List)
			print ("done")

if __name__ == "__main__":
	accepted_inputs = ['rho_x','rho_tx','rho_tpx','rho_y','rho_ty','rho_tpy']
	while True:
		text = input("Enter a type of .dat file (rho_x|rho_y|rho_tx|rho_ty|rho_tpx|rho_tpy): ")
		if text in accepted_inputs:
			break

	getItems().get_dat(text)
