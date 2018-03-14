#! /usr/bin/env python3
import os, os.path
import traceback
import contextlib
import math
import csv
class getItems():
	def __init__(self):
		self.possible_dattype = ['zeropi','pizero']

	def get_individual_type(self,pwd,each,dattype):
		individual_dir = "".join([pwd,'/',each,'/','stgfull','.dat'])

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
				items_removed.append(int(x[3:]))
 
		items_removed = sorted(items_removed)
		for index,truncated in enumerate(items_removed):
			items_removed[index] = "mu="+str(truncated)

		#gets a list of all directories
		pwd = os.getcwd() #your current working drive
		return_list = []
		for each in items_removed:
			# /Users/Alvin/Desktop/cy2001 results/varymuwhenv=5/mu=30/uni.dat
			#each = subfolder.
			#{ subfolder | mean |  sd }
			individual_dir = self.get_individual_type(pwd,each,dattype)
			print(individual_dir)
			dict_to_append = {}

			with open (individual_dir,"r") as datfile:
				counter = 0
				for item in datfile:
					if dattype== "zeropi":
						if counter==50:
							stgList=float(item)
					elif dattype == "pizero":
						if counter==101*50:
							stgList=float(item)
					counter+=1
								
				dict_to_append["subfolder"] = int(each[3:])
				#dict_to_append["mean"] = mean
				#dict_to_append["sd"] = sd
				dict_to_append["stg"] = stgList			

				return_list.append(dict_to_append)

		self.construct_csv(return_list,pwd,dattype)


	def get_csv_headers(self,dattype):
		if dattype == "zeropi":
			csv_List = [["mu","Stg(0,pi)"]]
		elif dattype == "pizero":
			csv_List = [["mu","Stg(pi,0)"]]
		return csv_List

	def get_file_name(self,pwd,dattype):
		csvfile = "".join([pwd,'/',dattype,'.csv'])

		return csvfile

	def construct_csv(self,array_of_calculations,pwd,dattype):
		csv_List = self.get_csv_headers(dattype)

		for row in array_of_calculations:
			csv_List.append([row['subfolder'],row['stg']])
		
		csvfilename = self.get_file_name(pwd,dattype)

		with open(csvfilename,'w') as csvFile:
			writer = csv.writer(csvFile)
			writer.writerows(csv_List)
			print ("done")

if __name__ == "__main__":
	accepted_inputs = ['zeropi','pizero']
	while True:
		text = input("Enter a type of .dat file (zeropi,pizero)): ")
		if text in accepted_inputs:
			break

	getItems().get_dat(text)
