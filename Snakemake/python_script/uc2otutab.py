 # ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Oct-22-2018
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for 16S microbiome analysis
# ################################### IMPORT ##################################### #

import os
import sys
import pandas
import numpy


def uclust_to_otutab_converter(uclust_file_Path, sample_delimiter, otutab_file_Path):
	#
	#new_df = pandas.DataFrame()
	#dfList = []
	#for chunk_DF in pandas.read_csv(uclust_file_Path, sep="\t", header=None, chunksize=1000000):
	uclust_DF = pandas.read_csv(uclust_file_Path, sep="\t", header=None)
	filtered_uclust_DF = uclust_DF.loc[uclust_DF[0] == 'H']
	selected_filtered_uclust_DF = filtered_uclust_DF[[8, 9]]
	sample_label_DF = selected_filtered_uclust_DF[8].str.split(sample_delimiter, 0, expand=True)
	sample_label_DF.drop(sample_label_DF.columns[1:], axis=1, inplace=True)

	sample_selected_filtered_uclust_DF = selected_filtered_uclust_DF.merge(sample_label_DF, left_index=True, right_index=True, how="inner")
	selected_sample_selected_filtered_uclust_DF = sample_selected_filtered_uclust_DF[[9, 0]]
	otutab_DF = pandas.crosstab(selected_sample_selected_filtered_uclust_DF[9], selected_sample_selected_filtered_uclust_DF[0], colnames=['0'])
	otutab_DF.index.names = ['#OTU ID']
	#dfList.append(otutab_DF)
	#new_df = pandas.concat(dfList, sort=False)

	otutab_DF.to_csv(otutab_file_Path, sep='\t', index=True, header=True)
	print("Conversion is successful")
	return True



uclust_to_otutab_converter(sys.argv[1], sys.argv[2], sys.argv[3])
