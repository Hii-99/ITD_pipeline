#!/bin/sh

set -euo pipefail

# Global Variable
BASHRC=~/.bashrc

shopt -s expand_aliases
source $BASHRC
source scripts/template/banner.sh

eval "$(conda shell.bash hook)"

VEP_ENV=master-vep
VERSION=0.1.0
FILE_NAME=unset
VCF_DIR=unset
OUT_DIR=./

verbose=false

usage(){
>&2 cat << EOF
Usage: $0
   [ -V | --version ]
   [ -v | --verbose ]
   [ -h | --help ]
   [ -f | --file_name args]
   [ -p | --prefix args]
   [ -o | --out_dir args; default=./ ]
EOF
}

args=$(getopt -a -o Vvhf:g:p:o: --long version,verbose,help,file_name:,prefix:,out_dir: -- "$@")

if [[ $? -gt 0 ]]; then
  usage
fi

eval set -- ${args}
while :
do
  case $1 in
    -V | --version)         echo $VERSION ; exit 1 ;;
    -v | --verbose)         verbose=true ; shift   ;;
    -h | --help)            wanglab_banner2 ;usage ; exit 1 ;;
    -f | --file_name)       FILE_NAME=$2 ; shift 2;;
    -p | --prefix)          VCF_DIR=$2 ; shift 2;;
    -o | --out_dir)         OUT_DIR=$2   ; shift 2 ;;

    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break ;;
    *) >&2 echo Unsupported option: $1
       usage ;;
  esac
done

if [ $verbose = true ]; then   
   wanglab_banner2
   >&2 echo "version                : ${VERSION}"
   >&2 echo "File Name              : ${FILE_NAME}"
   >&2 echo "Pindel VCF Directory   : ${VCF_DIR}"
   >&2 echo "Output Directory       : ${OUT_DIR}"
fi

conda activate $VEP_ENV

input_vcf=${PREFIX}/${FILE_ID}.pindel.vcf
output_vcf=${OUT_DIR}/${FILE_ID}.pindel.vep103.vcf
DIR_CACHE=~/bin/vep #vep lib
CACHE_VERSION=103
ASSEMBLY=GRCh38
THREAD=10
ASSEMBLY_FASTA=~/bin/gatk_index/Homo_sapiens_assembly38.fasta

vep -i $input_vcf -o $output_vcf --vcf --cache --dir_cache $DIR_CACHE --cache_version $CACHE_VERSION \
--assembly $ASSEMBLY --force_overwrite --fork $THREAD --offline --no_stats --everything -pick \
--fasta $ASSEMBLY_FASTA