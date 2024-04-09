#!/bin/sh

set -euo pipefail

# Global Variable
BASHRC=unset

shopt -s expand_aliases
# source $BASHRC
source scripts/template/banner.sh
VERSION=0.1.0
GENOMON_ITD_ENV=master-genomonITD
GENOMONITDetector38_DIR=~/bin/GenomonITDetector38
FILE_NAME=unset
PREFIX=""
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


# conda activate GENOMON_ITD_ENV

# ITD detection with genomon-ITDetector
cd $GENOMONITDetector38_DIR


./detectITD.sh $bam_file $WORK_DIR ${TCGA_ID}
wait # done genomonITD detection

# cleanup
cd ${CWD}
mv $WORK_DIR/itd_list.tsv ${OUT_DIR}/${FILE_NAME}_itd_list.tsv
rm -rf $WORK_DIR

