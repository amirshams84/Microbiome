
import os
import sys
import pandas


def usearch_to_mothur_converter(usearch_abundance_file_Path, title, mothur_abundance_file_Path):
	"""
	"""
	print("Reading abundance file...")
	abundance_DF = pandas.read_table(
		usearch_abundance_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		delimiter="\t",
		index_col="#OTU ID"
	)
	print("Done! :)")

	print("Generating MOTHUR abundance format...")
	OTU_name_List = abundance_DF.index.values.tolist()
	sample_name_List = abundance_DF.columns.values.tolist()

	transposed_abundance_DF = abundance_DF.transpose()
	transposed_abundance_DF = transposed_abundance_DF[OTU_name_List]
	NumOtus_value = str(len(OTU_name_List))
	NumOtus_List = [NumOtus_value] * len(sample_name_List)
	label_List = [title] * len(sample_name_List)
	shared_column_List = ['label', 'Groups', 'numOtus']
	shared_column_List.extend(OTU_name_List)
	transposed_abundance_DF = transposed_abundance_DF.assign(label=label_List, Groups=sample_name_List, numOtus=NumOtus_List)
	mothur_abundance_DF = pandas.DataFrame()
	mothur_abundance_DF = transposed_abundance_DF[shared_column_List]
	mothur_abundance_DF.to_csv(mothur_abundance_file_Path, sep='\t', index=False, header=True)
	print("Done! :)")
	print("Conversion is successful")
	return True


usearch_to_mothur_converter(sys.argv[1], sys.argv[2], sys.argv[3])