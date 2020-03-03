# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Oct-10-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for create cluster
# ################################### IMPORT ##################################### #


import os
import sys
import glob
import re
import pandas
library_path = os.path.abspath(workflow.basedir + "/../python_library/")
sys.path.append(library_path)
import utility
import factory
# ################################### FUNCTIONS ################################## #


def get_abundance_list(wildcards):
	"""
	"""
	abundance_List = []
	for sample, sample_Dict in METADATA_DICT.items():
		##
		abundance_List.append(THE_PATH + "/uparse/cluster_taxonomy/{sample}.uparse.LCR_abundance.txt".format(sample=sample))
	else:
		##
		pass
	return abundance_List


def get_taxonomy_list(wildcards):
	"""
	"""
	taxonomy_List = []
	for sample, sample_Dict in METADATA_DICT.items():
		##
		taxonomy_List.append(THE_PATH + "/uparse/cluster_taxonomy/{sample}.uparse.LCR_taxonomy.txt".format(sample=sample))
	else:
		##
		pass
	return taxonomy_List


def get_EC_abundance(wildcards):
	"""
	"""
	EC_abundance_List = []
	for sample, sample_Dict in METADATA_DICT.items():
		##
		EC_abundance_List.append(THE_PATH + "/uparse/picrust/{sample}.uparse.EC_abundance.txt".format(sample=sample))
	else:
		##
		pass
	return EC_abundance_List


def get_picrust_abundance_Dict():
	"""
	"""
	picrust_List = ["EC", "KO", "PWY"]
	picrust_abundance_Dict = {}
	for each_method in picrust_List:
		##
		picrust_abundance_Dict[each_method] = []
		for sample, sample_Dict in METADATA_DICT.items():
			##
			method_file = THE_PATH + "/uparse/picrust/{sample}.uparse.{method}_abundance.txt".format(sample=sample, method=each_method)

			if utility.is_file_exist(method_file) is True:
				#
				picrust_abundance_Dict[each_method].append(method_file)
		else:
			##for sample, sample_Dict in METADATA_DICT.items():
			pass
	else:
		##for each_method in picrust_List:
		pass
	#print(picrust_abundance_Dict)
	return picrust_abundance_Dict


def get_picrust_annotation_Dict():
	"""
	"""
	picrust_List = ["EC", "KO", "PWY"]
	picrust_annotation_Dict = {}
	for each_method in picrust_List:
		##
		picrust_annotation_Dict[each_method] = []
		for sample, sample_Dict in METADATA_DICT.items():
			##
			method_file = THE_PATH + "/uparse/picrust/{sample}.uparse.{method}_annotation.txt".format(sample=sample, method=each_method)

			if utility.is_file_exist(method_file) is True:
				#
				picrust_annotation_Dict[each_method].append(method_file)
		else:
			##for sample, sample_Dict in METADATA_DICT.items():
			pass
	else:
		##for each_method in picrust_List:
		pass
	#print(picrust_annotation_Dict)
	return picrust_annotation_Dict

# ################################### CONFIGURATION ############################## #


# +++++++++++++++++++++++++++++++++++++
#PATH
R_SCRIPT_PATH = os.path.abspath(workflow.basedir + "/../R_script")
PYTHON_SCRIPT_PATH = os.path.abspath(workflow.basedir + "/../python_script")
# -------------------------------------
# +++++++++++++++++++++++++++++++++++++
#GENERAL
config_general_Dict = config["GENERAL"]
PROJECT = config_general_Dict["PROJECT"]
EXPERIMENT = config_general_Dict["EXPERIMENT"]
TITLE = config_general_Dict["TITLE"]
WORKDIR = utility.fix_path(config_general_Dict["WORKDIR"])
DATADIR = utility.fix_path(config_general_Dict["DATADIR"])
REFDIR = utility.fix_path(config_general_Dict["REFDIR"])
THE_PATH = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#DATA
config_data_Dict = config["DATA"]
PLATFORM = config_data_Dict["PLATFORM"].lower()
FORMAT = config_data_Dict["FORMAT"].lower()
LAYOUT = config_data_Dict["LAYOUT"].lower()
SPECIES = config_data_Dict["SPECIES"].lower()
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#CLUSTER
PROCESSORS, MEMORY = utility.get_cluster_info(sys.argv)
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#METADATA
config_metadata_Dict = config["METADATA"]
METADATA_FILE = config_metadata_Dict["METADATA_FILE"]
SAMPLE_COLUMN = config_metadata_Dict["SAMPLE_COLUMN"]
METADATA_DICT = utility.build_metadata_dict(METADATA_FILE, SAMPLE_COLUMN)
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#USEARCH
config_usearch_Dict = config["USEARCH"]
USEARCH_FASTQ_FILTER = config_usearch_Dict["fastq_filter"][PLATFORM]
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
lineage_List = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
# ################################### WILDCARDS ################################## #


pre_process_List = []
pre_cluster_List = []
uparse_List = []
for sample, sample_Dict in METADATA_DICT.items():
	#
	pre_process_List.append(THE_PATH + "/pre_process/{sample}.pre_process.fastq".format(sample=sample))
	pre_process_List.append(THE_PATH + "/pre_process/{sample}.pre_process.fasta".format(sample=sample))
	pre_process_List.append(THE_PATH + "/pre_process/report/{sample}.pre_process.info".format(sample=sample))
	#
	pre_cluster_List.append(THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.fasta".format(sample=sample))
	pre_cluster_List.append(THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.filter_low_complexity.fasta".format(sample=sample))
	pre_cluster_List.append(THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.filter_low_complexity.filter_length.fasta".format(sample=sample))
	pre_cluster_List.append(THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.filter_low_complexity.filter_length.filter_duplicate.fasta".format(sample=sample))
	pre_cluster_List.append(THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.filter_low_complexity.filter_length.filter_duplicate.filter_taxonomy.fasta".format(sample=sample))
	#
	uparse_List.append(THE_PATH + "/uparse/cluster/{sample}.uparse.centroids.fasta".format(sample=sample))
	uparse_List.append(THE_PATH + "/uparse/cluster_alignment/{sample}.uparse.multi_hit.txt".format(sample=sample))
	uparse_List.append(THE_PATH + "/uparse/cluster_abundance/{sample}.uparse.abundance.txt".format(sample=sample))
	uparse_List.append(THE_PATH + "/uparse/cluster_taxonomy/{sample}.uparse.LCR_taxonomy.txt".format(sample=sample))
	uparse_List.append(THE_PATH + "/uparse/picrust/{sample}.uparse.picrust_execution.flag".format(sample=sample))
	uparse_List.append(THE_PATH + "/uparse/picrust/{sample}.uparse.EC_abundance.txt".format(sample=sample))
	#
	uparse_List.append(THE_PATH + "/uparse/aggregate/{title}.uparse.OTU_abundance.txt".format(title=TITLE))
	uparse_List.append(THE_PATH + "/uparse/aggregate/{title}.uparse.OTU_taxonomy.txt".format(title=TITLE))
	uparse_List.append(THE_PATH + "/uparse/aggregate/{title}.uparse.picrust.flag".format(title=TITLE))
	#
	#uparse_List.append(THE_PATH + "/uparse/{sample}.uparse.OTU_abundance.txt".format(sample=sample))
	
# ################################### PIPELINE FLOW ############################## #


rule Endpoint:
	"""
	"""
	input:
		pre_cluster_List +
		uparse_List
# ################################### PIPELINE RULES ############################# #

if FORMAT == "fastq":
	#
	rule filter_quality_fastq:
		"""
		"""
		input:
			pre_process_fastq = THE_PATH + "/pre_process/{sample}.pre_process.fastq",
			pre_process_fastq_info = THE_PATH + "/pre_process/report/{sample}.pre_process.info",
		output:
			filter_quality_fastq = THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.fastq",
			filter_quality_fasta = THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.fasta",
			filter_quality_info = THE_PATH + "/uparse/pre_cluster/report/{sample}.pre_process.filter_quality.info"
		message: "filter_quality_fastq: {PROJECT}|{EXPERIMENT}|{TITLE}|{wildcards.sample}"
		threads: PROCESSORS
		resources:
			mem_mb = MEMORY
		run:
			shell("""
				#
				##
				RESULT_PATH={THE_PATH}/uparse/pre_cluster/
				mkdir -p $RESULT_PATH
				
				REPORT_PATH={THE_PATH}/uparse/pre_cluster/report
				mkdir -p $REPORT_PATH
				##
				#
				
				module load usearch/11.0.667 || exit 1
				#######################################
				ASCII_BASE=$(python {PYTHON_SCRIPT_PATH}/query_usearch_fastx_info.py {input.pre_process_fastq_info} 'ascii_base')
				
				EE_HIGH_QUARTILE=$(python {PYTHON_SCRIPT_PATH}/query_usearch_fastx_info.py {input.pre_process_fastq_info} 'hi_quartile_ee')
				
				usearch -fastq_filter {input.pre_process_fastq} -fastq_ascii $ASCII_BASE {USEARCH_FASTQ_FILTER} -fastq_maxee $EE_HIGH_QUARTILE -fastqout {output.filter_quality_fastq}  -fastaout {output.filter_quality_fasta}
				usearch -fastx_info {output.filter_quality_fastq} -output {output.filter_quality_info}
			""")

elif FORMAT == "fasta":
	#
	pass


rule filter_low_complexity:
	"""
	"""
	input:
		filter_quality_fasta = THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.fasta"
	output:
		filter_low_complexity_fasta = THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.filter_low_complexity.fasta",
		filter_low_complexity_info = THE_PATH + "/uparse/pre_cluster/report/{sample}.pre_process.filter_quality.filter_low_complexity.info",
		discard_low_complexity_fasta = THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.low_complexity.fasta",
	
	message: "filter_low_complexity: {PROJECT}|{EXPERIMENT}|{TITLE}|{wildcards.sample}"
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			##
			RESULT_PATH={THE_PATH}/uparse/pre_cluster/
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={THE_PATH}/uparse/pre_cluster/report
			mkdir -p $REPORT_PATH

			SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
			mkdir -p $SCRATCH_PATH
			##
			#
			module load mothur/1.42.3 || exit 1
			module load usearch/11.0.667 || exit 1
			#######################################

			cd $SCRATCH_PATH
			
			mothur "#set.dir(tempdefault=$RESULT_PATH, \
			output=$RESULT_PATH); \
			set.logfile(name=$REPORT_PATH/{wildcards.sample}.filter_low_complexity.logfile); \
			set.current(processors={threads}); \
			screen.seqs(fasta={input.filter_quality_fasta}, maxhomop=8);"

			SCREEN_GOOD_OUTPUT=$RESULT_PATH/{wildcards.sample}.pre_process.filter_quality.good.fasta

			if [ -f $SCREEN_GOOD_OUTPUT ];
			then
				sed 's/[ \\t]*$//' $SCREEN_GOOD_OUTPUT > {output.filter_low_complexity_fasta}
			else
				echo "filter_low_complexity failed!!"
				cp {input.filter_quality_fasta} {output.filter_low_complexity_fasta}
			fi

			
			SCREEN_BAD_OUTPUT=$RESULT_PATH/{wildcards.sample}.pre_process.filter_quality.bad.accnos
			cut -f1 $SCREEN_BAD_OUTPUT > $REPORT_PATH/{wildcards.sample}.bad_read.txt
			usearch -fastx_getseqs {input.filter_quality_fasta} -labels $REPORT_PATH/{wildcards.sample}.bad_read.txt \\
			-fastaout {output.discard_low_complexity_fasta}
			usearch -fastx_info {output.filter_low_complexity_fasta} -output {output.filter_low_complexity_info}

			rm $SCREEN_GOOD_OUTPUT $SCREEN_BAD_OUTPUT $REPORT_PATH/{wildcards.sample}.bad_read.txt

		""")


rule filter_length:
	input:
		filter_low_complexity_fasta = THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.filter_low_complexity.fasta",
		filter_low_complexity_info = THE_PATH + "/uparse/pre_cluster/report/{sample}.pre_process.filter_quality.filter_low_complexity.info",
	output:
		filter_length_fasta = THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.filter_low_complexity.filter_length.fasta",
		filter_length_info = THE_PATH + "/uparse/pre_cluster/report/{sample}.pre_process.filter_quality.filter_low_complexity.filter_length.info",
	
	message: "filter_length: {PROJECT}|{EXPERIMENT}|{TITLE}|{wildcards.sample}"
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			##
			RESULT_PATH={THE_PATH}/uparse/pre_cluster
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={THE_PATH}/uparse/pre_cluster/report/
			mkdir -p $REPORT_PATH
			##
			#

			module load usearch/11.0.667 || exit 1
			#######################################


			LENGTH_MINUMUM=$(python {PYTHON_SCRIPT_PATH}/query_usearch_fastx_info.py {input.filter_low_complexity_info} 'min_length')
			LENGTH_LOW_QUARTILE=$(python {PYTHON_SCRIPT_PATH}/query_usearch_fastx_info.py {input.filter_low_complexity_info} 'lo_quartile_length')

			if (( $LENGTH_LOW_QUARTILE - $LENGTH_MINUMUM > 10 )); then
				LENGTH_THRESHOLD=$LENGTH_LOW_QUARTILE
			else
				LENGTH_THRESHOLD=$LENGTH_MINUMUM
			fi

			if [ {PLATFORM} == illumina ]; then
				LENGTH_THRESHOLD=$(( LENGTH_THRESHOLD - 1))
			elif [ {PLATFORM} == ion_torrent ]; then
				LENGTH_THRESHOLD=$(( LENGTH_THRESHOLD - 10 ))
			else
				LENGTH_THRESHOLD=$LENGTH_THRESHOLD
			fi
			echo "##############################"
			echo "LENGTH THRESHOLD: $LENGTH_THRESHOLD"
			echo "##############################"

			usearch -fastx_truncate {input.filter_low_complexity_fasta} -trunclen $LENGTH_THRESHOLD -fastaout {output.filter_length_fasta}
			usearch -fastx_info {output.filter_length_fasta} -output {output.filter_length_info}
		""")


rule filter_duplicate:
	input:
		filter_length_fasta = THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.filter_low_complexity.filter_length.fasta",
	output:
		filter_duplicate_fasta = THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.filter_low_complexity.filter_length.filter_duplicate.fasta",
		filter_duplicate_info = THE_PATH + "/uparse/pre_cluster/report/{sample}.pre_process.filter_quality.filter_low_complexity.filter_length.filter_duplicate.info",
	
	message: "filter_duplicate: {PROJECT}|{EXPERIMENT}|{TITLE}|{wildcards.sample}"
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			##
			RESULT_PATH={THE_PATH}/uparse/pre_cluster/
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={THE_PATH}/uparse/pre_cluster/report
			mkdir -p $REPORT_PATH
			##
			#
			module load usearch/11.0.667 || exit 1
			#######################################

			usearch -fastx_uniques {input.filter_length_fasta} -minuniquesize 2 -sizeout -strand both -fastaout {output.filter_duplicate_fasta}
			usearch -fastx_info {output.filter_duplicate_fasta} -output {output.filter_duplicate_info}
		""")


rule filter_taxonomy:
	input:
		filter_duplicate_fasta = THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.filter_low_complexity.filter_length.filter_duplicate.fasta",
	output:
		filter_taxonomy_fasta = THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.filter_low_complexity.filter_length.filter_duplicate.filter_taxonomy.fasta",
		filter_taxonomy_info = THE_PATH + "/uparse/pre_cluster/report/{sample}.pre_process.filter_quality.filter_low_complexity.filter_length.filter_duplicate.filter_taxonomy.info",
	
	message: "filter_taxonomy: {PROJECT}|{EXPERIMENT}|{TITLE}|{wildcards.sample}"
	priority: 996
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			##
			RESULT_PATH={THE_PATH}/uparse/pre_cluster
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={THE_PATH}/uparse/pre_cluster/report
			mkdir -p $REPORT_PATH
			##
			#

			module load usearch/11.0.667 || exit 1
			#######################################

			usearch -usearch_global {input.filter_duplicate_fasta} -db {REFDIR}/usearch/silva/silva.udb -id 0.80 -uc {output.filter_taxonomy_fasta}.uc -top_hit_only -maxhits 1 -strand both
			grep '^H' {output.filter_taxonomy_fasta}.uc | cut -f9 > {output.filter_taxonomy_fasta}.hit.label
			usearch -fastx_getseqs {input.filter_duplicate_fasta} -labels {output.filter_taxonomy_fasta}.hit.label -fastaout {output.filter_taxonomy_fasta}
			usearch -fastx_info {output.filter_taxonomy_fasta} -output {output.filter_taxonomy_info}

			rm {output.filter_taxonomy_fasta}.hit.label {output.filter_taxonomy_fasta}.uc
		""")


rule uparse_cluster:
	"""
	"""
	input:
		filter_taxonomy_fasta = THE_PATH + "/uparse/pre_cluster/{sample}.pre_process.filter_quality.filter_low_complexity.filter_length.filter_duplicate.filter_taxonomy.fasta",
	output:
		cluster_centroids = THE_PATH + "/uparse/cluster/{sample}.uparse.centroids.fasta",
		cluster_centroids_info = THE_PATH + "/uparse/cluster/report/{sample}.uparse.centroids.info",
	message: "uparse_cluster: {PROJECT}|{EXPERIMENT}|{TITLE}|{wildcards.sample}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"

			RESULT_PATH={THE_PATH}/uparse/cluster
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={THE_PATH}/uparse/cluster/report
			mkdir -p $REPORT_PATH
			##
			#
			module load usearch/11.0.667 || exit 1
			#######################################
			usearch -cluster_otus {input.filter_taxonomy_fasta} -otus {output.cluster_centroids} -uparseout {output.cluster_centroids_info}
		""")


rule uparse_cluster_alignment:
	"""
	"""
	input:
		pre_process_fasta = THE_PATH + "/pre_process/{sample}.pre_process.fasta",
		cluster_centroids = THE_PATH + "/uparse/cluster/{sample}.uparse.centroids.fasta",
	output:
		cluster_multi_hit_abundance = THE_PATH + "/uparse/cluster_alignment/{sample}.uparse.multi_hit.txt",
		cluster_multi_hit_abundance_uc = THE_PATH + "/uparse/cluster_alignment/report/{sample}.uparse.multi_hit.uc"
	
	message: "uparse_cluster_alignment: {PROJECT}|{EXPERIMENT}|{TITLE}|{wildcards.sample}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"

			RESULT_PATH={THE_PATH}/uparse/cluster_alignment
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={THE_PATH}/uparse/cluster_alignment/report
			mkdir -p $REPORT_PATH
			##
			#
			module load usearch/11.0.667 || exit 1
			#######################################
			Centroids_Count=$(grep -c '>' {input.cluster_centroids})

			usearch -usearch_global {input.pre_process_fasta} -db {input.cluster_centroids} -id 0.80 -maxhits 10 -top_hits_only \\
			-maxaccepts $Centroids_Count -maxrejects $Centroids_Count -strand both -uc {output.cluster_multi_hit_abundance_uc}

			python {PYTHON_SCRIPT_PATH}/uc2otutab.py {output.cluster_multi_hit_abundance_uc} ';' {output.cluster_multi_hit_abundance}
		""")


rule uparse_cluster_abundance:
	"""
	"""
	input:
		cluster_multi_hit_abundance = THE_PATH + "/uparse/cluster_alignment/{sample}.uparse.multi_hit.txt",
		cluster_centroids = THE_PATH + "/uparse/cluster/{sample}.uparse.centroids.fasta",
	output:
		cluster_abundance = THE_PATH + "/uparse/cluster_abundance/{sample}.uparse.abundance.txt",
		cluster_abundance_info = THE_PATH + "/uparse/cluster_abundance/report/{sample}.uparse.abundance.info",
	
	message: "uparse_cluster_abundance: {PROJECT}|{EXPERIMENT}|{TITLE}|{wildcards.sample}"
	resources:
		mem_mb = MEMORY
	run:
		
		shell("""
			module load usearch/11.0.667 || exit 1
			#######################################
			usearch -otutab_sortotus {input.cluster_multi_hit_abundance} -output {output.cluster_abundance}
			usearch -otutab_stats {output.cluster_abundance} -output {output.cluster_abundance_info}
		""")


rule uparse_cluster_taxonomy:
	"""
	"""
	input:
		cluster_centroids = THE_PATH + "/uparse/cluster/{sample}.uparse.centroids.fasta",
		cluster_abundance = THE_PATH + "/uparse/cluster_abundance/{sample}.uparse.abundance.txt",
	output:
		cluster_taxonomy_info = THE_PATH + "/uparse/cluster_taxonomy/report/{sample}.uparse.taxonomy.info",
		LCR_cluster_abundance = THE_PATH + "/uparse/cluster_taxonomy/{sample}.uparse.LCR_abundance.txt",
		LCR_cluster_taxonomy = THE_PATH + "/uparse/cluster_taxonomy/{sample}.uparse.LCR_taxonomy.txt",
		LCR_cluster_centroids = THE_PATH + "/uparse/cluster_taxonomy/{sample}.uparse.LCR_centroids.fasta",
	
	message: "uparse_cluster_taxonomy: {PROJECT}|{EXPERIMENT}|{TITLE}|{wildcards.sample}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"

			RESULT_PATH={THE_PATH}/uparse/cluster_taxonomy
			mkdir -p $RESULT_PATH
			
			REPORT_PATH={THE_PATH}/uparse/cluster_taxonomy/report
			mkdir -p $REPORT_PATH
			##
			#
			module load usearch/11.0.667 || exit 1
			#######################################
			usearch -sintax {input.cluster_centroids} -db {REFDIR}/usearch/rdp/rdp.udb -tabbedout {output.cluster_taxonomy_info}.tmp -strand both -sintax_cutoff 0.1

			sort -Vk1,1 {output.cluster_taxonomy_info}.tmp > {output.cluster_taxonomy_info}
			rm -rf {output.cluster_taxonomy_info}.tmp
		""")
		factory.calculate_taxonomy_lcr(output.cluster_taxonomy_info, output.cluster_taxonomy_info + ".lcr", "sintax")
		factory.aggregate_lowest_common_rank_abundance(output.cluster_taxonomy_info + ".lcr", input.cluster_abundance, output.LCR_cluster_abundance)
		factory.filter_lowest_common_rank_taxonomy(output.cluster_taxonomy_info, output.LCR_cluster_abundance, output.cluster_taxonomy_info)
		factory.relabel_taxonomy(output.cluster_taxonomy_info, output.cluster_taxonomy_info + ".sintax")
		shell("""
			module load usearch/11.0.667 || exit 1
			#######################################
			cut -f2 {output.cluster_taxonomy_info}.sintax > {output.cluster_taxonomy_info}.sintax.tax
			usearch -fastx_getseqs {REFDIR}/usearch/rdp/rdp.fasta -labels {output.cluster_taxonomy_info}.sintax.tax -fastaout {output.cluster_taxonomy_info}.sintax.tax.rdp.fasta -label_substr_match
			usearch -fastx_getlabels {output.cluster_taxonomy_info}.sintax.tax.rdp.fasta -output {output.cluster_taxonomy_info}.sintax.tax.rdp.id
			sed -i 's/;.*//' {output.cluster_taxonomy_info}.sintax.tax.rdp.fasta
		""")
		factory.get_uniq_taxa(output.cluster_taxonomy_info + ".sintax", output.cluster_taxonomy_info + ".sintax.tax.rdp.id", output.cluster_taxonomy_info + ".sintax.processed")
		factory.relabel_abundance(output.LCR_cluster_abundance, output.cluster_taxonomy_info + ".sintax.processed", output.LCR_cluster_abundance + ".final")
		shell("""
			cut -f2 {output.cluster_taxonomy_info}.sintax.processed > {output.cluster_taxonomy_info}.sintax.processed.id

			module load usearch/11.0.667 || exit 1
			#######################################
			usearch -fastx_getseqs {output.cluster_taxonomy_info}.sintax.tax.rdp.fasta -labels {output.cluster_taxonomy_info}.sintax.processed.id -fastaout {output.LCR_cluster_centroids}
			cut -f2,3 {output.cluster_taxonomy_info}.sintax.processed > {output.LCR_cluster_taxonomy}
			mv {output.LCR_cluster_abundance}.final {output.LCR_cluster_abundance}
		""")


rule uparse_aggregate_OTU:
	"""
	"""
	input:
		abundance_List = get_abundance_list,
		taxonomy_List = get_taxonomy_list,
	output:
		aggregate_abundance = THE_PATH + "/uparse/aggregate/{title}.uparse.OTU_abundance.txt",
		aggregate_abundance_info = THE_PATH + "/uparse/aggregate/{title}.uparse.OTU_abundance.info",
		aggregate_taxonomy = THE_PATH + "/uparse/aggregate/{title}.uparse.OTU_taxonomy.txt",
	
	message: "uparse_aggregate_OTU: {PROJECT}|{EXPERIMENT}|{TITLE}"
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	run:
		if os.path.exists(output.aggregate_abundance) is True:
			os.remove(output.aggregate_abundance)
		else:
			pass
		#
		for sample_abundance_sublist in [input.abundance_List[x:x + 100] for x in range(0, len(input.abundance_List), 100)]:
			#
			if os.path.exists(output.aggregate_abundance) is True:
				#
				sample_abundance_sublist.append(output.aggregate_abundance)
			else:
				pass
			alignment_file_String = ",".join(sample_abundance_sublist)
			shell("""
				module load usearch/11.0.667 || exit 1
			#######################################
				usearch -otutab_merge {alignment_file_String} -output {output.aggregate_abundance}.tmp
			""")
		else:
			pass
		#
		shell("""
			module load usearch/11.0.667 || exit 1
			#######################################
			usearch -otutab_sortotus {output.aggregate_abundance}.tmp -output {output.aggregate_abundance}
			usearch -otutab_stats {output.aggregate_abundance} -output {output.aggregate_abundance_info}

			rm -rf {output.aggregate_abundance}.tmp
		""")

		if os.path.exists(output.aggregate_taxonomy) is True:
			os.remove(output.aggregate_taxonomy)
		else:
			pass
		for sample_taxonomy_sublist in [input.taxonomy_List[x:x + 100] for x in range(0, len(input.taxonomy_List), 100)]:
			#
			if os.path.exists(output.aggregate_taxonomy) is True:
				#
				sample_taxonomy_sublist.append(output.aggregate_taxonomy)
			else:
				pass
			factory.concat_taxonomy_file_list(sample_taxonomy_sublist, output.aggregate_taxonomy)
		else:
			pass
		
		for each_lineage in lineage_List:
			##
			lineage_file = THE_PATH + "/uparse/aggregate/{title}.uparse.{method}_abundance.txt".format(title=TITLE, method=each_lineage)
			factory.aggregate_by_taxonomy(output.aggregate_abundance, output.aggregate_taxonomy, each_lineage, lineage_file)


rule uparse_picrust:
	"""
	"""
	input:
		LCR_cluster_centroids = THE_PATH + "/uparse/cluster_taxonomy/{sample}.uparse.LCR_centroids.fasta",
		LCR_cluster_abundance = THE_PATH + "/uparse/cluster_taxonomy/{sample}.uparse.LCR_abundance.txt",
	output:
		picrust_execution_flag = THE_PATH + "/uparse/picrust/{sample}.uparse.picrust_execution.flag"
	
	message: "uparse_picrust: {PROJECT}|{EXPERIMENT}|{TITLE}|{wildcards.sample}"
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	run:
		shell("""

			#
			##
			total_start_time="$(date -u +%s)"

			RESULT_PATH={THE_PATH}/uparse/picrust
			mkdir -p $RESULT_PATH
			
			##
			#

			module load picrust/2.3.0-b || exit 1
			#######################################
			picrust2_pipeline.py -s {input.LCR_cluster_centroids} -i {input.LCR_cluster_abundance} \\
			-o $RESULT_PATH/{wildcards.sample} -p {threads} --in_traits EC,KO --stratified --verbose

			add_descriptions.py -i $RESULT_PATH/{wildcards.sample}/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz -m EC \\
			-o $RESULT_PATH/{wildcards.sample}/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

			add_descriptions.py -i $RESULT_PATH/{wildcards.sample}/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz -m KO \\
			-o $RESULT_PATH/{wildcards.sample}/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz

			add_descriptions.py -i $RESULT_PATH/{wildcards.sample}/pathways_out/path_abun_unstrat.tsv.gz -m METACYC \\
			-o $RESULT_PATH/{wildcards.sample}/pathways_out/path_abun_unstrat_descrip.tsv.gz

			touch {output.picrust_execution_flag}
		""")


rule uparse_picrust_process:
	"""
	"""
	input:
		picrust_execution_flag = THE_PATH + "/uparse/picrust/{sample}.uparse.picrust_execution.flag",
	output:
		EC_abundance = THE_PATH + "/uparse/picrust/{sample}.uparse.EC_abundance.txt"

	message: "uparse_picrust_process: {PROJECT}|{EXPERIMENT}|{TITLE}|{wildcards.sample}"
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	run:
		#EC_abundance_zipped = THE_PATH + "/uparse/picrust/{sample}/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz".format(sample=wildcards.sample)
		EC_annotation_zipped = THE_PATH + "/uparse/picrust/{sample}/EC_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz".format(sample=wildcards.sample)
		EC_abundance = THE_PATH + "/uparse/picrust/{sample}.uparse.EC_abundance.txt".format(sample=wildcards.sample)
		#EC_annotation = THE_PATH + "/uparse/picrust/{sample}.uparse.EC_annotation.txt".format(sample=wildcards.sample)
		if utility.is_file_exist(EC_annotation_zipped) is True:
			#
			factory.process_annotation(EC_annotation_zipped, EC_abundance, compression_type="gzip")
			#utility.relabel_index(EC_abundance_zipped, EC_abundance, compression_type="gzip")
			#factory.relabel_annotation_column(EC_annotation_zipped, EC_annotation, compression_type="gzip")
			
		else:
			#
			pass
		#
		#KO_abundance_zipped = THE_PATH + "/uparse/picrust/{sample}/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz".format(sample=wildcards.sample)
		KO_annotation_zipped = THE_PATH + "/uparse/picrust/{sample}/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz".format(sample=wildcards.sample)
		KO_abundance = THE_PATH + "/uparse/picrust/{sample}.uparse.KO_abundance.txt".format(sample=wildcards.sample)
		#KO_annotation = THE_PATH + "/uparse/picrust/{sample}.uparse.KO_annotation.txt".format(sample=wildcards.sample)
		if utility.is_file_exist(KO_annotation_zipped) is True:
			#
			factory.process_annotation(KO_annotation_zipped, KO_abundance, compression_type="gzip")
			
		else:
			#
			pass
		#
		#PWY_abundance_zipped = THE_PATH + "/uparse/picrust/{sample}/pathways_out/path_abun_unstrat.tsv.gz".format(sample=wildcards.sample)
		PWY_annotation_zipped = THE_PATH + "/uparse/picrust/{sample}/pathways_out/path_abun_unstrat_descrip.tsv.gz".format(sample=wildcards.sample)
		PWY_abundance = THE_PATH + "/uparse/picrust/{sample}.uparse.PWY_abundance.txt".format(sample=wildcards.sample)
		#PWY_annotation = THE_PATH + "/uparse/picrust/{sample}.uparse.PWY_annotation.txt".format(sample=wildcards.sample)
		if utility.is_file_exist(PWY_annotation_zipped) is True:
			#
			factory.process_annotation(PWY_annotation_zipped, PWY_abundance, compression_type="gzip")

		else:
			#
			pass


rule uparse_aggregate_picrust_abundance:
	"""
	"""
	input:
		EC_abundance = get_EC_abundance
	output:
		picrust_aggregate_abundance_flag = THE_PATH + "/uparse/aggregate/{title}.uparse.picrust.flag",
	
	message: "uparse_aggregate_picrust_abundance: {PROJECT}|{EXPERIMENT}|{TITLE}"
	threads: PROCESSORS
	resources:
		mem_mb = MEMORY
	run:
		picrust_abundance_Dict = get_picrust_abundance_Dict()
		#picrust_annotation_Dict = get_picrust_annotation_Dict()
		for each_method in picrust_abundance_Dict:
			##
			aggregate_abundance = THE_PATH + "/uparse/aggregate/{title}.uparse.{method}_abundance.txt".format(title=TITLE, method=each_method)
			#aggregate_annotation = THE_PATH + "/uparse/aggregate/{title}.uparse.{method}_annotation.txt".format(title=TITLE, method=each_method)
			if os.path.exists(aggregate_abundance) is True:
				os.remove(aggregate_abundance)
			else:
				pass
			for abundance_sublist in [picrust_abundance_Dict[each_method][x:x + 100] for x in range(0, len(picrust_abundance_Dict[each_method]), 100)]:
				#
				if os.path.exists(aggregate_abundance) is True:
					#
					abundance_sublist.append(aggregate_abundance)
				else:
					pass
				abundance_sublist_String = ",".join(abundance_sublist)
				shell("""
					module load usearch/11.0.667 || exit 1
				#######################################
					usearch -otutab_merge {abundance_sublist_String} -output {aggregate_abundance}
				""")
			else:
				pass
			#
			'''
			shell("""
				module load usearch/11.0.667 || exit 1
				#######################################
				usearch -otutab_sortotus {aggregate_abundance} -output {aggregate_abundance}.tmp
				
				mv {aggregate_abundance}.tmp {aggregate_abundance}
				rm -rf {aggregate_abundance}.tmp
			""")
			#
			if os.path.exists(aggregate_annotation) is True:
				os.remove(aggregate_annotation)
			else:
				pass
			for annotation_sublist in [picrust_annotation_Dict[each_method][x:x + 100] for x in range(0, len(picrust_annotation_Dict[each_method]), 100)]:
				#
				if os.path.exists(aggregate_annotation) is True:
					#
					annotation_sublist.append(aggregate_annotation)
				else:
					pass
				factory.concat_taxonomy_file_list(annotation_sublist, aggregate_annotation)
			else:
				pass
			'''

		else:
			##for each_method in picrust_abundance_Dict:
			pass
		shell("""
			touch {output.picrust_aggregate_abundance_flag}
			""")
#######################################################################