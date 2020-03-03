# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Oct-10-2019
# Email: amir.shams84@gmail.com
# Aim: Factory
# ################################### IMPORT #################################### #


import os
import sys
import pandas
import utility
# ################################### FUNCTIONS ################################# #


def process_abundance_metadata(abundance_file_Path, metadata_file_Path, processed_abundance_file_Path, processed_metadata_file_Path):
	"""
	"""
	print("Reading metadata file...")
	metadata_DF = pandas.read_csv(
		metadata_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		sep="\t",
		index_col="Sample_ID"
	)

	#
	print("Reading abundance file...")
	abundance_DF = pandas.read_csv(
		abundance_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		sep="\t",
		index_col="#OTU ID"
	)
	#
	metadata_sample_name_List = metadata_DF.index.values.tolist()
	abundance_sample_name_List = abundance_DF.columns.values.tolist()
	if len(metadata_sample_name_List) <= len(abundance_sample_name_List):
		sample_name_List = metadata_sample_name_List
	else:
		sample_name_List = abundance_sample_name_List
	#
	processed_abundance_DF = abundance_DF[sample_name_List]
	metadata_DF = metadata_DF.loc[sample_name_List, ]
	print(str(os.path.basename(abundance_file_Path)) + " :('rows','columns'): " + str(processed_abundance_DF.shape))
	print("Generating USEARCH abundance format...")
	processed_abundance_DF.to_csv(processed_abundance_file_Path, sep='\t', index=True, header=True, encoding='utf-8')
	print("Done! :)")
	metadata_DF.to_csv(processed_metadata_file_Path, sep='\t', index=True, header=True, encoding='utf-8')
	return True


def process_taxonomy(taxonomy_file_Path, processed_taxonomy_file_Path):
	"""
	"""
	print("Reading taxonomy file...")
	sintax_Handle = open(taxonomy_file_Path, "rU")
	taxonomy_String = "#OTU ID" + "\t" + "Taxonomy" + "\n"
	for each_line in sintax_Handle:
		each_line = each_line.rstrip()
		each_line_List = each_line.split(";")
		OTU_ID = each_line_List[0]
		taxon_String = each_line_List[1]
		processed_Taxon = process_sintax_taxonomy(taxon_String)

		taxonomy_String += OTU_ID + "\t" + processed_Taxon + "\n"
	else:
		pass
	utility.write_string_down(taxonomy_String, processed_taxonomy_file_Path)
	return True


def process_taxonomy_old(taxonomy_file_Path, processed_taxonomy_file_Path):
	"""
	"""
	print("Reading taxonomy file...")
	sintax_Handle = open(taxonomy_file_Path, "rU")
	taxonomy_String = "#OTU ID" + "\t" + "Taxonomy" + "\n"
	for each_line in sintax_Handle:
		#
		each_line = each_line.rstrip()
		each_line_List = each_line.split("\t")
		OTU_ID = each_line_List[0]
		taxon_String = each_line_List[1]
		processed_Taxon = process_sintax_taxonomy(taxon_String)

		taxonomy_String += OTU_ID + "\t" + processed_Taxon + "\n"
	else:
		pass
	utility.write_string_down(taxonomy_String, processed_taxonomy_file_Path)
	return True


def calculate_taxonomy_lcr(taxonomy_file_Path, rectfied_taxonomy_file_Path, taxonomy_format):
	"""
	"""
	print("Reading taxonomy file...")
	taxonomy_HD = open(taxonomy_file_Path, "rU")
	taxonomy_String = ""
	for each_line in taxonomy_HD:
		##
		each_line = each_line.rstrip()
		each_line_List = each_line.split("\t")
		OTU_ID = each_line_List[0]
		taxon_String = each_line_List[1]
		if taxonomy_format == "sintax":
			#
			processed_Taxon = process_sintax_taxonomy(taxon_String)
		else:
			#
			pass
		taxonomy_String += OTU_ID + "\t" + processed_Taxon + "\n"
	else:
		##
		pass
	utility.write_string_down(taxonomy_String, rectfied_taxonomy_file_Path)
	return True


def process_sintax_taxonomy(sintax_String):
	"""
	d:Bacteria(1.0000),p:"Bacteroidetes"(1.0000),c:"Bacteroidia"(0.9800),o:"Bacteroidales"(0.9604),f:"Porphyromonadaceae"(0.8932),g:Coprobacter(0.1697),s:Coprobacter_fastidiosus(0.0322)
	"""
	
	processed_sintax_String = ""
	processed_sintax_List = []
	lineage_column_List = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
	sintax_List = sintax_String.split(",")
	for each_lineage in lineage_column_List:
		#
		lineage_index = lineage_column_List.index(each_lineage)

		try:
			lineage_Value = sintax_List[lineage_index]
			lineage_Value = lineage_Value.lower().split(":")[1].split("(")[0].capitalize()
			lineage_Value = lineage_Value.replace('"', '')
			lineage_Value = utility.slugify(lineage_Value)
			processed_sintax_List.append(lineage_Value)
		except IndexError:
			processed_sintax_List.append(lineage_Value + "_" + each_lineage)
	processed_sintax_String = ";".join(processed_sintax_List) + ";"

	return processed_sintax_String


def process_lowest_common_rank(sintax_file_Path, processed_sintax_file_Path):
	"""
	"""
	sintax_List = []
	f = open(sintax_file_Path, "r")
	sintax_String = ""
	for i in f:
		i = i.rstrip()
		line = i.split("\t")
		processed_Taxon = process_sintax_taxonomy(line[1])
		if processed_Taxon not in sintax_List:
			#
			sintax_List.append(processed_Taxon)
			sintax_String += i + "\n"
		else:
			pass
	else:
		pass
	utility.write_string_down(sintax_String, processed_sintax_file_Path)
	return True


def process_otu_fasta(otu_file_Path, taxonomy_file_Path, processed_otu_fasta_file_Path):
	"""
	"""
	taxonomy_Dict = {}
	f = open(taxonomy_file_Path, "r")
	for i in f:
		i = i.rstrip()
		line = i.split("\t")
		taxonomy_Dict[line[0]] = line[1]
	else:
		pass
	f.close()
	fasta_String = ""
	f = open(otu_file_Path, "r")
	for i in f:
		if i.startswith(">"):
			i = i.rstrip()
			fasta_String += ">" + taxonomy_Dict[i[1:]] + "\n"
		else:
			fasta_String += i
	else:
		pass
	utility.write_string_down(fasta_String, processed_otu_fasta_file_Path)
	return True


def convert_qiime_to_mothur(qiime_abundance_file_Path, mothur_abundance_file_Path, title):
	"""
	"""
	print("Reading abundance file...")
	qiime_DF = pandas.read_csv(
		qiime_abundance_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		sep="\t",
		index_col="#OTU ID"
	)
	OTU_name_List = qiime_DF.index.values.tolist()
	sample_name_List = qiime_DF.columns.values.tolist()
	transposed_qiime_DF = qiime_DF.transpose()
	transposed_qiime_DF = transposed_qiime_DF[OTU_name_List].fillna(0).astype(int)

	NumOtus_value = str(len(OTU_name_List))
	NumOtus_List = [NumOtus_value] * len(sample_name_List)
	label_List = [title] * len(sample_name_List)

	shared_column_List = ['label', 'Groups', 'numOtus']
	shared_column_List.extend(OTU_name_List)
	transposed_qiime_DF = transposed_qiime_DF.assign(label=label_List, Groups=sample_name_List, numOtus=NumOtus_List)
	print("Generating MOTHUR abundance format...")

	mothur_abundance_DF = pandas.DataFrame()
	mothur_abundance_DF = transposed_qiime_DF[shared_column_List]
	mothur_abundance_DF.to_csv(mothur_abundance_file_Path, sep='\t', index=False, header=True)
	print("Done! :)")
	return True


def process_metadata(metadata_file_Path, sample_column, treatment_column, design_String, mothur_metadata_file_Path, usearch_metadata_file_Path):
	"""
	"""
	sample_column = utility.slugify(sample_column)
	treatment_column = utility.slugify(treatment_column)
	design_String = utility.slugify(design_String)
	print("Reading metadata file...")
	metadata_DF = pandas.read_csv(
		metadata_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		sep="\t",
		index_col=None
	)
	print(str(os.path.basename(metadata_file_Path)) + " :('rows','columns'): " + str(metadata_DF.shape))
	print("processeing metadata format...")
	metadata_header_List = metadata_DF.columns.values.tolist()
	if sample_column not in metadata_header_List:
		print("sample_column: " + sample_column + "is not available in metadata file.")
		raise Exception("ABORTING!!!")
		sys.exit(2)
		return False
	else:
		pass
	if treatment_column not in metadata_header_List:
		print("treatment_column: " + treatment_column + "is not available in metadata file.")
		raise Exception("ABORTING!!!")
		sys.exit(2)
		return False
	else:
		pass
	processed_metadata_DF = metadata_DF[[sample_column, treatment_column]]
	print(processed_metadata_DF.head())
	metadata_Dict = {}
	metadata_Dict[sample_column] = "Sample_ID"
	metadata_Dict[treatment_column] = "Treatment"
	processed_metadata_DF = processed_metadata_DF.rename(columns=metadata_Dict)

	processed_metadata_DF = processed_metadata_DF[processed_metadata_DF['Treatment'] == design_String]
	processed_metadata_DF.to_csv(mothur_metadata_file_Path, sep='\t', index=False, header=True, encoding='utf-8')
	processed_metadata_DF.to_csv(usearch_metadata_file_Path, sep='\t', index=False, header=False, encoding='utf-8')
	print("Done! :)")
	return True


def relabel_fasta(fasta_file_path, label_file_Path, relabel_fasta_file_Path):
	"""
	"""
	label_Dict = {}
	f = open(label_file_Path)
	for i in f:
		i = i = i.rstrip()
		line = i.split("\t")
		label_Dict[line[0]] = line[1]
	else:
		pass
	f.close()
	#
	relabel_fasta_String = ""
	f = open(fasta_file_path)
	for i in f:
		if i.startswith(">"):
			header = i.rstrip()
			relabel_fasta_String += ">" + label_Dict[header[1:]] + "\n"
		else:
			relabel_fasta_String += i
	else:
		pass
	o = open(relabel_fasta_file_Path, "w")
	o.write(relabel_fasta_String)
	o.close()
	
	return True


def relabel_abundance(abundance_file_path, label_file_Path, relabel_abundance_file_Path):
	"""
	"""
	label_Dict = {}
	f = open(label_file_Path)
	for i in f:
		i = i = i.rstrip()
		line = i.split("\t")
		label_Dict[line[0]] = line[1]
	else:
		pass
	f.close()
	relabel_abundance_String = ""
	f = open(abundance_file_path)
	for i in f:
		if i.startswith("#OTU ID"):
			relabel_abundance_String += i
		else:
			i = i.rstrip()
			line = i.split("\t")
			relabel_abundance_String += label_Dict[line[0]] + "\t" + "\t".join(line[1:]) + "\n"
	else:
		pass
	o = open(relabel_abundance_file_Path, "w")
	o.write(relabel_abundance_String)
	o.close()
	return True


def relabel_taxonomy(taxonomy_file_Path, relabel_taxonomy_file_Path):
	"""
	"""
	relabel_taxonomy_String = ""
	f = open(taxonomy_file_Path, "r")
	for i in f:
		#
		i = i.rstrip()
		line = i.split("\t")
		otu_label = line[0]
		taxonomy_annotation = line[1]
		sintax_List = taxonomy_annotation.split(",")
		annotation_List = []
		for each_lineage in sintax_List:
			#
			annotation_List.append(each_lineage.split("(")[0])
		else:
			pass
		annotation_String = "tax=" + ",".join(annotation_List) + ";"
		relabel_taxonomy_String += otu_label + "\t" + annotation_String + "\n"
	else:
		pass
	o = open(relabel_taxonomy_file_Path, "w")
	o.write(relabel_taxonomy_String)
	o.close()
	return True


def get_uniq_taxa(centroids_taxonomy_file_Path, taxonomy_file_Path, unique_taxonomy_file_Path):
	"""
	"""
	unique_taxonomy_String = ""
	f = open(taxonomy_file_Path, "r")
	rdp_species_Dict = {}
	rdp_genus_Dict = {}
	rdp_family_Dict = {}
	for i in f:
		i = i.rstrip()
		rdp_ID = i.split(";")[0]
		taxonomy = i.split(";")[1]
		if "s:" in taxonomy:
			try:
				rdp_species_Dict[taxonomy] = rdp_ID
			except KeyError:
				pass
		elif "g:" in taxonomy:
			try:
				rdp_genus_Dict[taxonomy] = rdp_ID
			except KeyError:
				pass
		else:
			try:
				rdp_family_Dict[taxonomy] = rdp_ID
			except KeyError:
				pass
	else:
		pass
	f.close()
	f = open(centroids_taxonomy_file_Path, "r")
	for i in f:
		#
		i = i.rstrip()
		line = i.split("\t")
		OTU_ID = line[0]
		OTU_taxonomy = line[1].split(";")[0]
		try:
			rdp_ID = rdp_species_Dict[OTU_taxonomy]

		except KeyError:
			try:
				rdp_ID = rdp_genus_Dict[OTU_taxonomy]
			except KeyError:
				rdp_ID = rdp_family_Dict[OTU_taxonomy]

		unique_taxonomy_String += OTU_ID + "\t" + rdp_ID + "\t" + process_sintax_taxonomy(OTU_taxonomy) + "\n"
	f.close()
	o = open(unique_taxonomy_file_Path, "w")
	o.write(unique_taxonomy_String)
	o.close()
	return True


def aggregate_metadata(metadata_file_List, aggregate_metadata_file_Path):
	"""
	"""
	metadata_dataframe_List = []
	for each_metadata in metadata_file_List:
		#
		each_metadata_DF = pandas.read_csv(each_metadata,
			encoding=None,
			skip_blank_lines=True,
			error_bad_lines=False,
			sep="\t",
			index_col=None
		)
		metadata_dataframe_List.append(each_metadata_DF)
	aggregate_metadata_DF = pandas.concat(metadata_dataframe_List)
	aggregate_metadata_DF.to_csv(aggregate_metadata_file_Path, sep='\t', index=False, header=True, encoding='utf-8')
	return True


def aggregate_lowest_common_rank_abundance(taxonomy_file_Path, abundance_file_Path, lowest_common_rank_abundance_file_Path):
	"""
	"""
	taxonomy_Dict = {}
	f = open(taxonomy_file_Path, "r")
	for i in f:
		##
		i = i.rstrip()
		line = i.split("\t")
		if line[1] not in taxonomy_Dict:
			#
			taxonomy_Dict[line[1]] = [line[0]]
		else:
			#
			taxonomy_Dict[line[1]].append(line[0])
	else:
		##
		pass
	f.close()
	
	abundance_DF = pandas.read_csv(abundance_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		sep="\t",
		index_col="#OTU ID"
	)

	OTU_name_List = abundance_DF.index.values.tolist()
	#sample_name_List = abundance_DF.columns.values.tolist()

	transposed_abundance_DF = abundance_DF.transpose()
	transposed_abundance_DF = transposed_abundance_DF[OTU_name_List].astype(int)

	lowest_common_rank_abundance_DF = pandas.DataFrame()
	lowest_common_rank_abundance_DF.index.names = ['#OTU ID']
	for taxonomy, otu_List in taxonomy_Dict.items():
		#
		
		if len(otu_List) > 1:
			#
			lowest_common_rank_abundance_DF[otu_List[0]] = transposed_abundance_DF[otu_List].sum(axis=1)
		else:
			#
			lowest_common_rank_abundance_DF[otu_List[0]] = transposed_abundance_DF[otu_List[0]]
	else:
		pass

	transposed_lowest_common_rank_abundance_DF = lowest_common_rank_abundance_DF.transpose()
	transposed_lowest_common_rank_abundance_DF.index.names = ['#OTU ID']
	transposed_lowest_common_rank_abundance_DF.to_csv(lowest_common_rank_abundance_file_Path, sep='\t', index=True, header=True, encoding='utf-8')
	
	return True


def filter_lowest_common_rank_taxonomy(taxonomy_file_Path, LCR_abundance_file_Path, LCR_taxonomy_file_Path):
	"""
	"""
	abundance_DF = pandas.read_csv(
		LCR_abundance_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		sep="\t",
		index_col="#OTU ID"
	)

	taxonomy_DF = pandas.read_csv(taxonomy_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		sep="\t",
		header=None,
		index_col=None
	)
	taxonomy_DF.columns = ["OTU_ID", "sintax", "strand", "tax"]
	taxonomy_DF = taxonomy_DF.set_index("OTU_ID")
	
	filter_taxonomy_DF = taxonomy_DF[taxonomy_DF.index.isin(abundance_DF.index)]
	
	filter_taxonomy_DF.to_csv(LCR_taxonomy_file_Path, sep='\t', index=True, header=False, quoting=3)
	return True


def convert_usearch_to_mothur_distance(usearch_distance, index_column_string, mothur_distance):
	"""
	"""
	distance_DF = pandas.read_table(usearch_distance, encoding=None, skip_blank_lines=True, error_bad_lines=False, delimiter="\t", index_col=index_column_string)
	file_Handle = open(mothur_distance, "w")
	file_Handle.write(str(distance_DF.shape[0]) + '\n')
	file_Handle.close()
	distance_DF.to_csv(mothur_distance, sep='\t', index=True, header=False, encoding=None, mode='a')
	return True


def concat_taxonomy_file_list(taxonomy_file_list, target_file_Path):
	"""
	"""
	taxonomy_Dict = {}
	for taxonomy_file in taxonomy_file_list:
		##
		f = open(taxonomy_file, "r")
		for i in f:
			##
			if i.startswith("#OTU ID"):
				continue
			i = i.rstrip()
			line = i.split("\t")
			if line[0] not in taxonomy_Dict:
				#
				taxonomy_Dict[line[0]] = line[1]
			else:
				pass
		else:
			##for i in f:
			pass
		f.close()
	else:
		##for taxonomy_file in taxonomy_file_list:
		pass
	taxonomy_String = "#OTU ID" + "\t" + "Annotation" + "\n"
	for each_tax_ID in taxonomy_Dict:
		##
		taxonomy_String += each_tax_ID + "\t" + taxonomy_Dict[each_tax_ID] + "\n"
	else:
		##for each_tax_ID in taxonomy_Dict:
		pass
	utility.write_string_down(taxonomy_String, target_file_Path)
	return True


def aggregate_by_taxonomy(usearch_abundance_file, taxonomy_file, lineage, lineage_aggregate_file):
	"""
	"""
	lineage_column_List = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
	if lineage not in lineage_column_List:
		#
		print("Not recognized lineage")
		print("Please choose one of these: " + str(lineage_column_List))
		return False
	else:
		pass
	taxonomy_df = pandas.read_csv(taxonomy_file, sep="\t", index_col="#OTU ID")
	#df.head()
	processed_taxonomy_df = build_lineage_dataframe(taxonomy_df)
	#processed_df.head()
	otu_df = pandas.read_csv(usearch_abundance_file, sep="\t", index_col="#OTU ID")
	taxonomic_abundance_DF = otu_df.merge(processed_taxonomy_df, left_index=True, right_index=True, how="inner")
	lineage_df = taxonomic_abundance_DF.groupby(lineage, as_index=True).sum()
	lineage_column = lineage_df.columns.values.tolist()[0:5]
	lineage_df = lineage_df.sort_values(lineage_column, ascending=False)
	lineage_df.index.name = "#OTU ID"
	lineage_df.to_csv(lineage_aggregate_file, sep='\t', index=True, header=True)
	#print("Aggregate was successful")
	return True


def build_lineage_dataframe(taxonomy_DF):
	"""
	"""
	max_group_length = taxonomy_DF["Annotation"].str.count(";").max()
	
	lineage_column_List = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
	#lineage_column_List.append("LCA")
	lineage_DF = taxonomy_DF['Annotation'].str.split(';', max_group_length, expand=True)
	lineage_DF.drop(lineage_DF.columns[len(lineage_DF.columns) - 1], axis=1, inplace=True)
	lineage_DF.columns = lineage_column_List[:max_group_length]
	return lineage_DF


def relabel_annotation_column(annotation_file_Path, relabeled_annotation_file_Path, compression_type=None):
	"""
	"""
	dataframe_DF = pandas.read_csv(
		annotation_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		delimiter="\t",
		compression=compression_type,
		index_col=0
	)
	label = dataframe_DF.columns.values.tolist()
	subset_dataframe_DF = dataframe_DF[[label[0]]]
	subset_dataframe_DF.index.name = "#OTU ID"
	
	#label_Dict = {}
	#label_Dict[label[0]] = "Annotation"
	#print(label_Dict)
	subset_dataframe_DF.columns = ["Annotation"]
	#dataframe_DF.rename(columns=label_Dict, inplace=True)
	#print(subset_dataframe_DF.head())
	subset_dataframe_DF.to_csv(
		relabeled_annotation_file_Path,
		sep='\t',
		index=True,
		header=True
	)
	return True


def process_annotation(annotation_file_Path, processed_annotation_file_Path, compression_type=None):
	"""
	"""
	annotation_DF = pandas.read_csv(
		annotation_file_Path,
		encoding=None,
		skip_blank_lines=True,
		error_bad_lines=False,
		sep="\t",
		compression=compression_type,
		index_col=None
	)
	label = annotation_DF.columns.values.tolist()
	annotation_DF["#OTU ID"] = annotation_DF[label[0]].astype(str) + "_" + annotation_DF[label[1]].astype(str)
	annotation_DF = annotation_DF.drop(columns=[label[0], label[1]])
	annotation_DF.set_index("#OTU ID", inplace=True)

	annotation_DF = annotation_DF.fillna(0).astype(int)

	annotation_DF.to_csv(
		processed_annotation_file_Path,
		sep='\t',
		index=True,
		header=True
	)
	return True


########################################################################

