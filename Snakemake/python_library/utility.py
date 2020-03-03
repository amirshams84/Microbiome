# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Oct-10-2019
# Email: amir.shams84@gmail.com
# Aim: Utility
# ################################### IMPORT #################################### #


import os
import sys
import pandas
import re
import pickle
import functools
# ################################### FUNCTIONS ################################# #


def is_file_exist(file_Path):
	"""
	"""
	if os.path.isfile(file_Path) and os.path.exists(file_Path) and os.access(file_Path, os.R_OK):
		return True
	else:
		return False


def read_file(the_file_Path):
	"""
	"""
	f = open(the_file_Path, "r")
	read_binary_target = f.read()
	f.close()
	return read_binary_target


def write_file(the_String, the_file_Path):
	#
	f = open(the_file_Path, "w")
	f.write(the_String)
	f.close()
	return True


def write_string_down(the_String, the_file_Path):
	#
	f = open(the_file_Path, "w")
	f.write(the_String)
	f.close()
	return True


def fix_path(the_Path):
	"""
	"""
	if the_Path[-1] == "/":
		the_Path = the_Path[:-1]
	return the_Path


def slugify(target_String):
	#replacing non-AlphNumeric with underscore
	slugified_String = re.sub('[^0-9a-zA-Z]+', '_', target_String)
	#slugified_String = slugified_String.replace("___", "__")
	return slugified_String


def merge_dict(first_dict, second_dict):
	"""
	"""
	metadata_Dict = {
		**first_dict,
		**second_dict
	}

	return metadata_Dict


def dict_to_pickle_converter(the_Dict, pickle_file_Path):
	"""
	"""
	pickle_Handle = open(pickle_file_Path, "wb")
	pickle.dump(the_Dict, pickle_Handle)
	pickle_Handle.close()
	return True


def build_metadata_dict(metadata_file_Path, index_column):
	"""
	"""
	if is_file_exist(metadata_file_Path) is True:
		pass
	metadata_Dict = {}

	metadata_DF = pandas.read_csv(
		metadata_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		delimiter="\t",
		index_col=index_column
	)
	transposed_metadata_DF = metadata_DF.transpose()
	metadata_Dict = transposed_metadata_DF.to_dict()
	return metadata_Dict


def get_cluster_info(parameters_List):
	"""
	"""
	for each_element in parameters_List:
		#
		if "--cluster=" not in each_element:
			#
			processors = 10
			memory = 10000
			return(processors, memory)
		if "--cluster=" in each_element:
			cluster_String = each_element
			cluster_List = cluster_String.split(" ")
			for each_cluster_element in cluster_List:
				if "--cpus-per-task=" in each_cluster_element:
					processors = int(each_cluster_element.split("=")[1])
				elif "--mem=" in each_cluster_element:
					memory = int(each_cluster_element.split("=")[1])
				else:
					pass
			else:
				pass
	else:
		pass
	return(processors, memory)


def merge_dataframe_file_list(dataframe_file_List, target_file_Path):
	"""
	"""
	dataframe_List = []
	for each_file in dataframe_file_List:
		##
		dataframe_DF = pandas.read_csv(
			each_file,
			encoding=None,
			skip_blank_lines=True,
			error_bad_lines=False,
			delimiter="\t",
			index_col="Distance"
		)
		dataframe_List.append(dataframe_DF)
	else:
		##for each_file in dataframe_file_List:
		pass
	
	merged_dataframe = functools.reduce(lambda x, y: pandas.merge(x, y, on='Distance'), dataframe_List)
	merged_dataframe.to_csv(
		target_file_Path,
		sep='\t',
		index=True,
		header=True
	)
	return True


def relabel_index(source_file_Path, target_file_Path, compression_type=None):
	"""
	"""
	dataframe_DF = pandas.read_csv(
		source_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		delimiter="\t",
		compression=compression_type,
		index_col=0
	)
	dataframe_DF.index.name = "#OTU ID"

	dataframe_DF = dataframe_DF.fillna(0).astype(int)
	dataframe_DF.to_csv(
		target_file_Path,
		sep='\t',
		index=True,
		header=True
	)
	return True


def concat_dataframe_file_list(dataframe_file_List, target_file_Path):
	"""
	"""
	dataframe_List = []
	for each_file in dataframe_file_List:
		##
		dataframe_DF = pandas.read_csv(
			each_file,
			encoding=None,
			skip_blank_lines=True,
			error_bad_lines=False,
			delimiter="\t",
			index_col=None
		)
		dataframe_List.append(dataframe_DF)
	else:
		##for each_file in dataframe_file_List:
		pass
	
	concat_dataframe = pandas.concat(dataframe_List, axis=0)
	concat_dataframe.to_csv(
		target_file_Path,
		sep='\t',
		index=False,
		header=True
	)
	return True


def print_dictionary(the_Dict):
	"""
	"""
	for x in the_Dict:
		print (x)
		for y in the_Dict[x]:
			print ("\t", y,':',the_Dict[x][y], "\n")
	return True
