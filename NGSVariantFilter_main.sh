#!/bin/bash

############################################################
#
# Assembly Combine and filter gVCF Files
# Cathal Ormond 2017. 
#
# 
# 
# 
# usage: ./main.sh -f BASE [-i] [-g]
#
# -f BASE
#     BASE is the base file name for the GWAS files
# -i
#    Run interactively, i.e. confirm plink values at each 
#    step of the analysis
# -g
#	Create graphics
#
############################################################



### Pipefail
set -o pipefail



### Global directories and files
# The raw data lies here. 
# BED or PED/MAP files must live here. 
SOURCE_DIR=${PWD}

LOG_DIR=${SOURCE_DIR}/logs
RESULTS_DIR=${SOURCE_DIR}/results
TEMP_DIR=${SOURCE_DIR}/temp/
GRAPHICS_DIR=${SOURCE_DIR}/graphics

SCRIPT_DIR=/home/shared/scripts
CURR_SCRIPT_DIR=${SCRIPT_DIR}/ngs-variant-filter

TOOL_DIR=/home/shared/tools
REF_DIR=/home/shared/reference/ReferenceGenome
REFERENCE=GRCh38

MAIN_LOG_FILE=${LOG_DIR}/main.log



# Get functions defined elsewhere
. ${SCRIPT_DIR}/func.sh




# Log directory for all log files 
if [ ! -d ${LOG_DIR} ]
then
	mkdir -p ${LOG_DIR}
fi


if [ -f ${MAIN_LOG_FILE} ]
then
	mv ${MAIN_LOG_FILE} ${MAIN_LOG_FILE}.$(date +%Y-%m-%d_%H.%M.%S)
	touch ${MAIN_LOG_FILE}
	log "Renaming previous log file" 4 
else
	touch ${MAIN_LOG_FILE}
	log "Creating the main log file: ${MAIN_LOG_FILE}" 4  
fi


### Argument processing
BUILD=GRCh38
CLEAN=false
DNSNP=""
E=100
F=0
GATK_FILE=${TOOL_DIR}/gatk/GenomeAnalysisTK_3.8.jar
JAVA_MAX="6G"
JAVA_MIN="6G"
OUTPUT_FILE="Variant_Filter"
S=0
NPROCS=4  # change to appropriate default
VERBOSE=3
SNP_TRANCHE_LEVEL=99.9
INDEL_TRANCHE_LEVEL=99.0
MAF_FILTER=0.02
HWE_FILTER=1e-6
LMISS_FILTER=0.02
FAM_FILE=""


cmd(){
	echo `basename $0`
}


usage(){
echo "\
`cmd` [OPTION ...]"

echo -e "\
-b, --build ; hg19 or GRCh38 ; [GRCh38]
-c, --clean; ; Remove all input files once finished; [${CLEAN}] 
-d, --dbsnp; <FILE>; dbSNP VCF file; [v 150]
-e, --end; <INT>; Tool to finish on; [${E}]
-f, --file; <FILE>; Name of source file
-g, --gatk; <FILE>; GATK file; [GATK 3.8]
-h, --help; ; Output this message
-m, --mem-min; <INT>; Java minimum memory; [${JAVA_MIN}]
-o, --out; <FILE>; Output file prefix; [${OUTPUT_FILE}] 
-r, --reference; <FILE>; FASTA Reference file; [GRCh38]
-s, --start; <INT>; Tool to start on; [${S}]
-t, --threads; <INT>; Number of threads for multithreaded processes; [${NPROCS}]
-v, --verbose; <INT>; Set verbosity level; [${VERBOSE}]
-x, --mem-max; <INT>; Java maximum memory; [${JAVA_MAX}]
    --snpTS; <FLOAT>; Tranche Sensitivity for SNPs; [${SNP_TRANCHE_LEVEL}]
    --indelTS; <FLOAT>; Tranche Sensitivity for indels; [${INDEL_TRANCHE_LEVEL}]
    --fam ; <FILE> ; FAM file ; []
" | column -t -s ";"
}


OPTS=`getopt -o b:cd:e:f:g:hm:o:r:s:t:v:x: \
	--long build:,file:,start:,end:,threads:,help,verbose:,dbsnp:,reference:,mem-min:,mem-max:,clean,out:,snpTS:,indelTS:,gatk:,fam: \
	-n '$(cmd)' -- "$@"`

if [ $? != 0 ]
then 
	echo "Error with arguments. Terminating ..." >&2
	exit 1
fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$OPTS"



while true; do
	case "$1" in
		-a | --apply)
			APPLY=true
			shift
			;;

		-b | --build)
			BUILD="${2}"
			shift 2
			;;

		-c | --clean)
			CLEAN=true
			shift
			;;

		-d | --dbsnp)
			DBSNP="$2"
			shift 2
			;;

		-e | --end)
			E="$2"
			shift 2
			;;

		-f | --file)
			SOURCE_FILE="$2"
			shift 2
			;;

		-g | --gatk)
			GATK_FILE="$2"
			shift 2
			;;

		-h | --help)
			usage
			exit 0
			;;

		-m | --mem-min)
			JAVA_MIN="$2"
			shift 2
			;;

		-o | --out)
			OUTPUT_FILE="$2"
			shift 2
			;;

		-p | --ped)
			PED_FILE="$2"
			shift 2
			;;

		-r | --reference)
			REF_FASTA="$2"
			shift 2
			;;

		-s | --start)
			S="$2"
			shift 2
			;;

		-t | --threads)
			NPROCS="$2"
			shift 2
			;;

		-v | --verbose)
			VERBOSE="$2"
			shift 2
			;;

		-x | --mem-max)
			JAVA_MAX="$2"
			shift 2
			;;

		--snpTS)
			SNP_TRANCHE_LEVEL="$2"
			shift 2
			;;

		--indelTS)
			INDEL_TRANCHE_LEVEL="$2"
			shift 2
			;;

		--calculate)
			CALCULATE=true
			shift 
			;;

		--fam)
			FAM_FILE="$2"
			shift 2
			;;

		--)
			shift
			break
			;;

		\?)
			err "Invalid flags. Exiting ... "
			exit 12
			break
			;;

		:)
			err "Flag requires an argument. Exiting ... "
			exit 13
			break
			;;

		*)
			err "Error with flags. Exiting ... "
			exit 14
			break
			;;
	esac
done





### Call the initialisation script

log "Calling initialisation script" 3 
log " " 3
. ${CURR_SCRIPT_DIR}/initialise.sh ${SOURCE_FILE} 

if [ $? -ne 0 ]
then
	log "Error - initialisation script returned an error: $?" 1
	exit 11
fi




### Main section of pipeline
log $(printf '#%.0s' $(seq 1 $(($(tput cols)-35)) ) ) 3 2> /dev/null
log " " 3 2>/dev/null
log "Main pipeline section" 3
log " " 3


if [[ ${APPLY} = true ]]
then
	log "Applying VQSR to all variants. " 0
	log "SNP Tranche Sensitivity: ${SNP_TRANCHE_LEVEL} " 0
	log "Indel Tranche Sensitivity: ${INDEL_TRANCHE_LEVEL} " 0
else
	log "Calculating VQSLOD scores and tranches only" 0
fi

log " " 3
log " " 3


CURR=0


# Index the reference FASTA File 
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	NAME="01_joint"
	log $NAME 3
	(. ${CURR_SCRIPT_DIR}/$NAME.sh $NAME) 
	testResult $? $NAME
fi



# Create the reference dictionary for the reference FASTA file
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	NAME="02_dp"
	log $NAME 3
	(. ${CURR_SCRIPT_DIR}/$NAME.sh $NAME) 
	testResult $? $NAME
fi



# Create the reference dictionary for the reference FASTA file
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	NAME="03_split"
	log $NAME 3
	(. ${CURR_SCRIPT_DIR}/$NAME.sh $NAME) 
	testResult $? $NAME
fi



# Map the reads and convert to a BAM file
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	NAME="04_recalSNP"
	log $NAME 3
	(. ${CURR_SCRIPT_DIR}/$NAME.sh $NAME) 
	testResult $? $NAME
fi



# Edit the Read Group
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	NAME="05_recalIND"
	log $NAME 3
	(. ${CURR_SCRIPT_DIR}/$NAME.sh $NAME) 
	testResult $? $NAME
fi



# Reorder the BAM file
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	NAME="06_applySNP"
	log $NAME 3
	(. ${CURR_SCRIPT_DIR}/$NAME.sh $NAME) 
	testResult $? $NAME
fi



# Sort the BAM file
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	NAME="07_applyIND"
	log $NAME 3
	(. ${CURR_SCRIPT_DIR}/$NAME.sh $NAME) 
	testResult $? $NAME
fi



# Validate the BAM file
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	NAME="08_hardOTH"
	log $NAME 3
	(. ${CURR_SCRIPT_DIR}/$NAME.sh $NAME) 
	testResult $? $NAME
fi



# Create recalibration table for BQSR
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	NAME="09_merge"
	log $NAME 3
	(. ${CURR_SCRIPT_DIR}/$NAME.sh $NAME) 
	testResult $? $NAME
fi



# Create recalibration table for BQSR
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	NAME="10_validate"
	log $NAME 3
	(. ${CURR_SCRIPT_DIR}/$NAME.sh $NAME) 
	testResult $? $NAME
fi



# Create recalibration table for BQSR
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	NAME="11_table"
	log $NAME 3
	(. ${CURR_SCRIPT_DIR}/$NAME.sh $NAME) 
	testResult $? $NAME
fi



# Create recalibration table for BQSR
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	NAME="12_rename"
	log $NAME 3
	(. ${CURR_SCRIPT_DIR}/$NAME.sh $NAME) 
	testResult $? $NAME
fi



# Create recalibration table for BQSR
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	NAME="13_files"
	log $NAME 3
	(. ${CURR_SCRIPT_DIR}/$NAME.sh $NAME) 
	testResult $? $NAME
fi



# Create recalibration table for BQSR
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	NAME="14_filter"
	log $NAME 3
	(. ${CURR_SCRIPT_DIR}/$NAME.sh $NAME) 
	testResult $? $NAME
fi



# Create recalibration table for BQSR
CURR=$(( CURR + 1 ))
if [[ ${S} -le ${CURR} && ${E} -ge ${CURR} ]]
then
	NAME="15_report"
	log $NAME 3
	(. ${CURR_SCRIPT_DIR}/$NAME.sh $NAME) 
	testResult $? $NAME
fi




