# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Oct-10-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for preprocessing input data
# ################################### IMPORT ##################################### #


import os
import sys
import glob
import re
library_path = os.path.abspath(workflow.basedir + "/../python_library/")
sys.path.append(library_path)
import utility
# ################################### FUNCTIONS ################################## #


def get_fastq(wildcards):
	"""
	sample delimiter will prevent from mixing sample name
	"""
	return glob.glob(DATADIR + "/" + wildcards.sample + SAMPLE_DELIMITER + "*" + "." + SAMPLE_SUFFIX, recursive=True)
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
THE_PATH = WORKDIR + "/" + PROJECT + "/" + EXPERIMENT + "/" + TITLE
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#CLUSTER
PROCESSORS, MEMORY = utility.get_cluster_info(sys.argv)
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#DATA
config_data_Dict = config["DATA"]
PLATFORM = config_data_Dict["PLATFORM"].lower()
FORMAT = config_data_Dict["FORMAT"].lower()
LAYOUT = config_data_Dict["LAYOUT"].lower()
SPECIES = config_data_Dict["SPECIES"].lower()
SAMPLE_DELIMITER = config_data_Dict["SAMPLE_DELIMITER"]
SAMPLE_SUFFIX = config_data_Dict["SAMPLE_SUFFIX"].lower()
####
if SAMPLE_SUFFIX[0] == ".":
	SAMPLE_SUFFIX = SAMPLE_SUFFIX[1:]
if SAMPLE_SUFFIX not in ["fastq", "fastq.gz"]:
	print("This pipeline can only be used with fastq or fastq.gz data")
	print("Aborting!!!!")
	sys.exit(2)
else:
	pass
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#METADATA
config_metadata_Dict = config["METADATA"]
METADATA_FILE = config_metadata_Dict["METADATA_FILE"]
SAMPLE_COLUMN = config_metadata_Dict["SAMPLE_COLUMN"]
METADATA_DICT = utility.build_metadata_dict(METADATA_FILE, SAMPLE_COLUMN)
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#PRE_PROCESS
config_pre_process_Dict = config["PRE_PROCESS"]
PRE_PROCESS_CUTADAPT = config_pre_process_Dict["CUTADAPT"]
# ------------------------------------
# ################################### WILDCARDS ################################## #


pre_process_List = []
for sample, sample_Dict in METADATA_DICT.items():
	#
	pre_process_List.append(THE_PATH + "/pre_process/{sample}.pre_process.fastq".format(sample=sample))
# ################################### PIPELINE FLOW ############################## #


rule Endpoint:
	"""
	"""
	input:
		pre_process_List
# ################################### PIPELINE RULES ############################# #


if FORMAT == "fastq":
	#
	rule pre_process_fastq:
		"""
		"""
		input:
			fastq_List = get_fastq
		output:
			pre_process_fastq = THE_PATH + "/pre_process/{sample}.pre_process.fastq",
			pre_process_fasta = THE_PATH + "/pre_process/{sample}.pre_process.fasta",
			pre_process_info = THE_PATH + "/pre_process/report/{sample}.pre_process.info"
		message: "pre_process_fastq: {PROJECT}|{EXPERIMENT}|{TITLE}|{wildcards.sample}"
		threads: PROCESSORS
		resources:
			mem_mb = MEMORY
		run:
			for each_fastq in input.fastq_List:
				#
				each_fastq_basename = os.path.basename(each_fastq)
				each_fastq_begining = re.sub("." + SAMPLE_SUFFIX, "", each_fastq_basename)
				shell("""
					#
					##
					total_start_time="$(date -u +%s)"

					RESULT_PATH={THE_PATH}/pre_process
					mkdir -p $RESULT_PATH
					
					REPORT_PATH={THE_PATH}/pre_process/report
					mkdir -p $REPORT_PATH
					##
					#
					
					module load cutadapt/2.5 || exit 1
					#######################################
					cutadapt {PRE_PROCESS_CUTADAPT} --cores={threads} {each_fastq} >> {output.pre_process_fastq}.tmp 2> $REPORT_PATH/{each_fastq_begining}.cutadapt.txt

				""")
			else:
				pass
			shell("""
				#
				##
				total_start_time="$(date -u +%s)"

				RESULT_PATH={THE_PATH}/pre_process
				mkdir -p $RESULT_PATH
				
				REPORT_PATH={THE_PATH}/pre_process/report
				mkdir -p $REPORT_PATH
				##
				#
				
				module load usearch/11.0.667 || exit 1
				module load fastqc/0.11.8 || exit 1
				#######################################
				usearch -fastq_filter {output.pre_process_fastq}.tmp -relabel '{wildcards.sample};Read_' -fastqout {output.pre_process_fastq} -fastaout {output.pre_process_fasta}

				fastqc -o $REPORT_PATH -f fastq --threads {threads} {output.pre_process_fastq}
				usearch -fastx_info {output.pre_process_fastq} -output {output.pre_process_info}
				usearch -fastx_uniques {output.pre_process_fastq} -uc $REPORT_PATH/{wildcards.sample}.fastx_uniques.report -strand both
				Unique_Count=$(grep -c "^S" $REPORT_PATH/{wildcards.sample}.fastx_uniques.report)
				READ_COUNT=$(python {PYTHON_SCRIPT_PATH}/query_usearch_fastx_info.py {output.pre_process_info} 'read_count')

				echo "READ_COUNT: $READ_COUNT" >> {output.pre_process_info}
				echo "UNIQUE_COUNT: $Unique_Count">> {output.pre_process_info}

				rm -rf {output.pre_process_fastq}.tmp
			""")

elif FORMAT == "fasta":
	#
	pass

