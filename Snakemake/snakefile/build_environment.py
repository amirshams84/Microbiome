# ################################### INFO ###################################### #
# Author: Amir Shams
# Date: Oct-10-2019
# Email: amir.shams84@gmail.com
# Aim: Snakemake workflow for buidling environment
# ################################### IMPORT ##################################### #

import os
import sys
library_path = os.path.abspath(workflow.basedir + "/../python_library/")
sys.path.append(library_path)
import utility
# ################################### FUNCTIONS ################################## #
# ################################### CONFIGURATION ############################## #

# ++++++++++++++++++++++++++++++++++++
#GENERAL
config_general_Dict = config["GENERAL"]
REFDIR = utility.fix_path(config_general_Dict["REFDIR"])
# -----------------------------------
# ++++++++++++++++++++++++++++++++++++
#CLUSTER
PROCESSORS, MEMORY = utility.get_cluster_info(sys.argv)
# ------------------------------------
# ++++++++++++++++++++++++++++++++++++
#REFERENCE
config_reference_Dict = config["REFERENCE"]
RDP_DB = config_reference_Dict["USEARCH"]["RDP16"]["SPECIES"]
SILVA_DB = config_reference_Dict["USEARCH"]["SILVA123"]["GENERAL"]
# ------------------------------------
# ################################### WILDCARDS ################################## #

build_environment_List = [REFDIR + "/build_database.flag"]
# ################################### PIPELINE FLOW ############################## #


rule endpoint:
	input:
		build_environment_List
# ################################### PIPELINE RULES ############################# #


rule build_database:
	"""
	"""
	output:
		#
		usearch_rdp_fasta = REFDIR + "/usearch/rdp/rdp.fasta",
		usearch_rdp_udb = REFDIR + "/usearch/rdp/rdp.udb",
		#
		usearch_silva_fasta = REFDIR + "/usearch/silva/silva.fasta",
		usearch_silva_udb = REFDIR + "/usearch/silva/silva.udb",
		#
		build_database_flag = REFDIR + "/build_database.flag",
	priority: 999
	threads: PROCESSORS
	message: "build_database: {PROJECT}|{EXPERIMENT}|{TITLE}"
	resources:
		mem_mb = MEMORY
	run:
		shell("""
			#
			##
			total_start_time="$(date -u +%s)"
			mkdir -p {REFDIR}/usearch/rdp
			mkdir -p {REFDIR}/usearch/silva

			SCRATCH_PATH=/lscratch/${{SLURM_JOB_ID}}
			mkdir -p $SCRATCH_PATH

			module load usearch/11.0.667 || exit 1

			cd $SCRATCH_PATH
			wget {RDP_DB} -O ${{PWD}}/usearch_rdp.gz
			gunzip < ${{PWD}}/usearch_rdp.gz > {output.usearch_rdp_fasta}
			usearch -makeudb_usearch {output.usearch_rdp_fasta} -output {output.usearch_rdp_udb}
			rm -rf ${{PWD}}/usearch_rdp.gz

			wget {SILVA_DB} -O ${{PWD}}/usearch_silva.gz
			gunzip < ${{PWD}}/usearch_silva.gz > {output.usearch_silva_fasta}
			usearch -makeudb_usearch {output.usearch_silva_fasta} -output {output.usearch_silva_udb}
			rm -rf ${{PWD}}/usearch_silva.gz

			touch {output.build_database_flag}

			
		""")
