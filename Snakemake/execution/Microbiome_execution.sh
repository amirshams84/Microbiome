#! /bin/bash
set -o pipefail
set -e
args=("$@")
if [ $# -eq 0 ]; then
    echo "No Input file provided!!!"
    echo "Aborting!!"
    exit 1

fi
####################################################
#

declare -a json_List
function parse_json {
	
	index=0
	IFS="&&"
	export PYTHONIOENCODING=utf8
	for line in $(cat $1 | python -c 'import os,sys,json; \
		data = json.load(sys.stdin); \
		print(\
		os.path.expanduser(data["GENERAL"]["WORKDIR"])+"&&"+\
		str(data["GENERAL"]["PROJECT"])+"&&"+\
		str(data["GENERAL"]["EXPERIMENT"])+"&&"+\
		str(data["GENERAL"]["TITLE"])
		)')
	do
		#echo $index
		#echo $line
		json_List[$index]=$line
		index=$(($index+1))
		
	done
	IFS=" "
	
};

function read_link() {
    local path=$1
    if [ -d $path ] ; then
        local abspath=$(cd $path; pwd -P)
    else
        local prefix=$(cd $(dirname -- $path) ; pwd -P)
        local suffix=$(basename $path)
        local abspath="$prefix/$suffix"
    fi
    if [ -e $abspath ] ; then
        echo $abspath
    else
        echo "$1 is not accessible"
		echo "make sure the path is correct(Please use absoulte path)!"
		echo "Aborting!!"
    	exit 1
    fi
};
####################################################
#
#Parse json file
JSON_CONFIG_FILE=${args[0]}
#
parse_json $JSON_CONFIG_FILE
WORK_DIR=${json_List[0]}
mkdir -p $WORK_DIR
WORK_DIR=$(read_link $WORK_DIR)
#
PROJECT=${json_List[2]}
EXPERIMENT=${json_List[4]}
TITLE=${json_List[6]}
#
LOG_DIR=$WORK_DIR/${PROJECT}/${EXPERIMENT}/${TITLE}/logs
mkdir -p $LOG_DIR
LOG_DIR=$(read_link $LOG_DIR)
#
MAIN_DIR=/data/shamsaddinisha/Development/Microbiome
SNAKEMAKE_DIR=${MAIN_DIR}/Snakemake
SNAKE_FILE_DIR=$SNAKEMAKE_DIR/snakefile
CLUSTER_CONFIG=${MAIN_DIR}/Snakemake/cluster_config
####################################################
#

EXECUTION_MODE=${args[1]}

if [ "$EXECUTION_MODE" == "DEVELOPMENT" ]
then
	PROCESSORS=10
	MEMORY=10000
elif [ "$EXECUTION_MODE" == "TEST" ]
then
	PROCESSORS=10
	MEMORY=50000
	CLUSTER_CONFIG_FILE=${CLUSTER_CONFIG}/cluster_test.yaml
elif [ "$EXECUTION_MODE" == "DEPLOYMENT" ]
then
	PROCESSORS=30
	MEMORY=100000
	CLUSTER_CONFIG_FILE=${CLUSTER_CONFIG}/cluster_deploy.yaml
fi
####################################################
#

module load snakemake || exit 1
####################################################
#


if [ "$EXECUTION_MODE" == "DEVELOPMENT" ]
then
	
	snakemake --snakefile $SNAKE_FILE_DIR/build_environment.py --configfile $JSON_CONFIG_FILE --unlock

	snakemake --snakefile $SNAKE_FILE_DIR/build_environment.py --configfile $JSON_CONFIG_FILE --keep-going --rerun-incomplete --local-cores=$PROCESSORS --cores=$PROCESSORS

	snakemake --snakefile $SNAKE_FILE_DIR/pre_process.py --configfile $JSON_CONFIG_FILE --keep-going --rerun-incomplete --local-cores=$PROCESSORS --cores=$PROCESSORS

	snakemake --snakefile $SNAKE_FILE_DIR/uparse.py --configfile $JSON_CONFIG_FILE --keep-going --rerun-incomplete --local-cores=$PROCESSORS --cores=$PROCESSORS

else
	#

	mkdir -p ${LOG_DIR}/build_environment
	snakemake --snakefile $SNAKE_FILE_DIR/build_environment.py --configfile $JSON_CONFIG_FILE --unlock
	snakemake --snakefile $SNAKE_FILE_DIR/build_environment.py --configfile $JSON_CONFIG_FILE --cluster-config ${CLUSTER_CONFIG_FILE} --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete --cluster="sbatch --cpus-per-task=$PROCESSORS --mem=$MEMORY --mincpus=$PROCESSORS \
	--partition={cluster.partition} --time={cluster.time} --mail-type=FAIL --job-name={cluster.jobname} \
	--output=${LOG_DIR}/build_environment/{cluster.output} --error=${LOG_DIR}/build_environment/{cluster.error} {cluster.extra}"
	
	#

	mkdir -p ${LOG_DIR}/pre_process
	snakemake --snakefile $SNAKE_FILE_DIR/pre_process.py --configfile $JSON_CONFIG_FILE --cluster-config ${CLUSTER_CONFIG_FILE} --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete --cluster="sbatch --cpus-per-task=$PROCESSORS --mem=$MEMORY --mincpus=$PROCESSORS \
	--partition={cluster.partition} --time={cluster.time} --mail-type=FAIL --job-name={cluster.jobname} \
	--output=${LOG_DIR}/pre_process/{cluster.output} --error=${LOG_DIR}/pre_process/{cluster.error} {cluster.extra}"

	#

	mkdir -p ${LOG_DIR}/uparse
	snakemake --snakefile $SNAKE_FILE_DIR/uparse.py --configfile $JSON_CONFIG_FILE --cluster-config ${CLUSTER_CONFIG_FILE} --local-cores=$PROCESSORS --cores=$PROCESSORS \
	--max-jobs-per-second=100 --latency-wait=120 --keep-going --rerun-incomplete --cluster="sbatch --cpus-per-task=$PROCESSORS --mem=$MEMORY --mincpus=$PROCESSORS \
	--partition={cluster.partition} --time={cluster.time} --mail-type=FAIL --job-name={cluster.jobname} \
	--output=${LOG_DIR}/uparse/{cluster.output} --error=${LOG_DIR}/uparse/{cluster.error} {cluster.extra}"
	
	

fi
#
####################################################