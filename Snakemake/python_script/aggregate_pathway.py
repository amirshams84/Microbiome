import os,sys
import pandas


def aggregate_pathway_level(top_level_file_Path, sec_level_file_Path, third_level_file_Path, aggregate_file_Path):
	"""
	"""
	top_DF = pandas.read_csv(
		top_level_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		delimiter="\t",
		index_col="pathway"
	)
	top_DF = top_DF[["description"]]

	top_DF.rename(columns = {'description':'L1'}, inplace = True)
	sec_DF = pandas.read_csv(
		sec_level_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		delimiter="\t",
		index_col="pathway"
	)
	sec_DF = sec_DF[["description"]]
	sec_DF.rename(columns = {'description':'L2'}, inplace = True)
	third_DF = pandas.read_csv(
		third_level_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		delimiter="\t",
		index_col="pathway"
	)
	third_DF = third_DF[["description"]]

	third_DF.rename(columns = {'description':'L3'}, inplace = True)
	aggregate_DF = pandas.concat([top_DF, sec_DF, third_DF], join='outer', axis=1)

	merge_DF = aggregate_DF.assign(Annotation = aggregate_DF.L1.astype(str) + ';' + aggregate_DF.L2.astype(str) + ';' + aggregate_DF.L3.astype(str) + ';')

	merge_DF = merge_DF[["Annotation"]]

	merge_DF.to_csv(aggregate_file_Path, sep='\t', index=True, header=True)

	return True


aggregate_pathway_level(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
