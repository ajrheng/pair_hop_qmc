#! /usr/bin/env python3
import os, os.path
import traceback
import contextlib
import math
import csv
class getItems():
	def __init__(self):
		self.possible_dattype = ['rho','rhotp','rhot']

	def get_individual_type(self,pwd,each,dattype):
		individual_dir = "".join([pwd,'/',each,'/',dattype,'.dat'])

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
				m=float(x[3:])
				if m%1==0:#if it is whole number
					items_removed.append(int(m))
				else:#if mu has decimal eg. 13.4
					items_removed.append(m)


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
				rightlist = []
				leftlist = []
				datfilelist = datfile.read().split()
				counter = 0
				for item in datfilelist:
					if counter%2 :
						rightlist.append(float(item))
					else:
						leftlist.append(float(item))
					counter +=1
				#if (dattype == "rho" or dattype == "rhot" or dattype == "rhotp"):
				#	rightlist = [3/4*x for x in rightlist]
				#	leftlist = [3/4*x for x in leftlist]
				# we want to calculate mean here
				rightmean = self.calculate_mean(rightlist)
				rightsd = self.standard_deviation(rightlist)
				dict_to_append["subfolder"] = float(each[3:])
				#dict_to_append["mean"] = mean
				#dict_to_append["sd"] = sd
				if dattype in self.possible_dattype:
					leftmean = self.calculate_mean(leftlist)
					leftsd = self.standard_deviation(leftlist)
					#dict_to_append["leftmean"] = leftmean
					#dict_to_append["leftsd"] = leftsd
					mean =rightmean+leftmean
					sd = math.sqrt(rightsd*rightsd+leftsd*leftsd)
					dict_to_append["mean"] = mean
					dict_to_append["sd"] = sd
				elif dattype == "uni":
					dict_to_append["mean"] = rightmean
					dict_to_append["sd"] = rightsd
				elif dattype == "stgpipi":
					leftmean = self.calculate_mean(leftlist)
					leftsd = self.standard_deviation(leftlist)
					dict_to_append["mean"] = leftmean
					dict_to_append["sd"] = leftsd

				return_list.append(dict_to_append)

		self.construct_csv(return_list,pwd,dattype)


	def calculate_mean(self,array_of_items):
		length_array = len(array_of_items)
		sum_of_items = float(0)
		for each in array_of_items:
			sum_of_items += float(each)
		#here we have the combined sum,
		mean = sum_of_items/length_array
		return mean

	def standard_deviation(self,array_of_items):
		#calculate mean, go through all items, for each item.
		#minus mean, square sum it divide 90
		mean = self.calculate_mean(array_of_items)
		counter = 0
		for index,value in enumerate(array_of_items):
			value = value-mean
			value *= value
			array_of_items[index] = value
			counter += 1
		sum_of_squared = sum(array_of_items)
		divide_of_sum = float(sum_of_squared/(counter*(counter-1)))
		sqrt_of_divided = math.sqrt(divide_of_sum)
		return sqrt_of_divided

	def get_csv_headers(self,dattype):
		if dattype == "uni":
			csv_List = [["mu","Mean","STDEV"]]
		elif dattype == "rho":
			csv_List = [["mu","rho mean","STDEV"]]
		elif dattype == "rhot":
			csv_List = [["mu","rhot mean","STDEV"]]
		elif dattype =="rhotp":
			csv_List = [["mu","rhotp mean","STDEV"]]
		elif dattype == "stgpipi":
			csv_List = [["mu","ssa mean","STDEV"]]
		elif dattype == "curr_single":
			csv_List = [["mu","single current","STDEV"]]
		elif dattype == "curr_pair":
			csv_List = [["mu","pair current","STDEV"]]

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
	accepted_inputs = ['uni','rho','rhot','rhotp','stgpipi']
	while True:
		text = input("Enter a type of .dat file (uni|rho|rhot|rhotp|stgpipi): ")
		if text in accepted_inputs:
			break

	getItems().get_dat(text)
