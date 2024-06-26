#!/bin/sh

set -euo pipefail

# Global Variable
BASHRC=~/.bashrc

shopt -s expand_aliases
source $BASHRC
source scripts/template/banner.sh

eval "$(conda shell.bash hook)"

VERSION=0.1.1
THREAD=20
REF=~/bin/GenomonITDetector38/GRCh38.d1.vd1.fa
REF_NAME=GRCh38.d1.vd1.fa
PINDEL_ENV=master-pindel
FILE_NAME=unset
PREFIX=unset
TCGA_ID=unset
OUT_DIR=./
verbose=false

usage(){
>&2 cat << EOF
Usage: $0
   [ -V | --version ]
   [ -v | --verbose ]
   [ -h | --help ]
   [ -f | --file_name args ]
   [ -p | --prefix args]
   [ -t | --TCGA_ID  args ]
   [ -o | --out_dir args; default=./ ]
EOF
}

args=$(getopt -a -o Vvht:f:p:o: --long version,verbose,help,TCGA_ID:,file_name:,prefix:,out_dir: -- "$@")

if [[ $? -gt 0 ]]; then
  usage
fi

eval set -- ${args}
while :
do
  case $1 in
    -V | --version)        echo $VERSION ; exit 1 ;;
    -v | --verbose)        verbose=true ; shift   ;;
    -h | --help)           wanglab_banner2 ;usage ; exit 1 ;;
    -f | --file_name)      FILE_NAME=$2 ; shift 2;;
    -p | --prefix)         PREFIX=$2 ; shift 2;;
    -t | --TCGA_ID)        TCGA_ID=$2 ; shift 2;;
    -o | --out_dir)        OUT_DIR=$2   ; shift 2 ;;

    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break ;;
    *) >&2 echo Unsupported option: $1
       usage 
       exit 1;;
  esac
done

if [  "$PREFIX" = "" ]; then
    BAM_FILE=${FILE_NAME}.bam
else
    BAM_FILE=${PREFIX}/${FILE_NAME}.bam
fi

CWD=$(pwd)
WORK_DIR=`mktemp -d -p "$CWD"`

if [ $verbose = true ]; then   
   wanglab_banner2
   >&2 echo "version          : ${VERSION}"
   >&2 echo "TCGA ID          : ${TCGA_ID}" 
   >&2 echo "File Name        : ${FILE_NAME}"
   >&2 echo "BAM File         : ${BAM_FILE}"
   >&2 echo "Output Directory : ${OUT_DIR}"
fi

conda activate $PINDEL_ENV

CONFIG_FILE=${WORK_DIR}/${FILE_NAME}_bam_config.txt
if [ ! -f $CONFIG_FILE ]
        then
                touch $CONFIG_FILE
                echo -e "${BAM_FILE}\t250\t${TCGA_ID}" >> $CONFIG_FILE
        else

                echo "${TCGA_ID}: ${file_ID} done bam config"
        fi

pindel -f $REF -i $CONFIG_FILE -o ${WORK_DIR}/${FILE_NAME} -T ${THREAD}
pindel2vcf -p ${WORK_DIR}/${FILE_NAME}_TD -r ${REF} -R ${REF_NAME} -d  $(date +'%Y%m%d') -v ${OUT_DIR}/${FILE_NAME}_TD.vcf &
pindel2vcf -p ${WORK_DIR}/${FILE_NAME}_SI -r ${REF} -R ${REF_NAME} -d  $(date +'%Y%m%d') -v ${OUT_DIR}/${FILE_NAME}_SI.vcf &
wait

rm -r ${WORK_DIR}/${FILE_NAME}_*

rm -rf $WORK_DIR
