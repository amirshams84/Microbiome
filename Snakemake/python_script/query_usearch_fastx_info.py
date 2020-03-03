 # ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Oct-22-2018
# Email: amir.shams84@gmail.com
# Aim: Python script to extract specifc info from fastx_info report
# ################################### IMPORT ##################################### #

import os
import sys


def convert_si_to_number(x):
	total_stars = 0
	if 'k' in x:
		if len(x) > 1:
			total_stars = float(x.replace('k', '')) * 1000 # convert k to a thousand
	elif 'M' in x:
		if len(x) > 1:
			total_stars = float(x.replace('M', '')) * 1000000 # convert M to a million
	elif 'B' in x:
		total_stars = float(x.replace('B', '')) * 1000000000 # convert B to a Billion
	else:
		total_stars = int(x) # Less than 1000
	
	return int(total_stars)


def query_usearch_fastx_info(fastx_info_file_Path, query):
	"""
	"""
	fastx_info_Dict = {}
	f = open(fastx_info_file_Path, "r")
	for index, i in enumerate(f):
		i = i.rstrip()
		if index == 0:
			line = i.split(",")
			fastx_info_Dict["file_size"] = line[0].split(" ")[2]
			fastx_info_Dict["read_count"] = convert_si_to_number(line[1].split(" ")[1])
			fastx_info_Dict["letter_count"] = line[2].split(" ")[1]
		elif index == 1:
			line = i.split(",")
			fastx_info_Dict["min_length"] = line[0].split(" ")[2]
			fastx_info_Dict["lo_quartile_length"] = line[1].split(" ")[2]
			fastx_info_Dict["median_length"] = line[2].split(" ")[2]
			fastx_info_Dict["hi_quartile_length"] = line[3].split(" ")[2]
			fastx_info_Dict["max_length"] = line[4].split(" ")[2]
		elif index == 4:
			fastx_info_Dict["ascii_base"] = i.split("=")[1]
		elif index == 5:
			line = i.split(",")
			fastx_info_Dict["mean_ee"] = line[0].split(";")[0].split(" ")[2]
			fastx_info_Dict["min_ee"] = line[0].split(";")[1].split(" ")[2]
			fastx_info_Dict["lo_quartile_ee"] = line[1].split(" ")[2]
			fastx_info_Dict["median_ee"] = line[2].split(" ")[2]
			fastx_info_Dict["hi_quartile_ee"] = line[3].split(" ")[2]
			fastx_info_Dict["max_ee"] = line[4].split(" ")[2]
		
		else:
			pass
	f.close()
	if query in fastx_info_Dict:
		print(fastx_info_Dict[query])
	else:
		print("N/A")
	return True
	

query_usearch_fastx_info(sys.argv[1], sys.argv[2])
