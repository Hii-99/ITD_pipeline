#!/bin/sh

# Global Variable
BASHRC= /home/hiiluann99/.bashrc
OUT_DIR=~/data_10T/



############### preprocessing ###############


shopt -s expand_aliases
source $BASHRC
source banner.sh



cwd=$(pwd)
cd OUT_DIR

[[ ! -d "download" ]] && mkdir -p "download"
[[ ! -d "tmp" ]] && mkdir -p "tmp"


gdc_download(){

        [[ ! -d "tmp/download" ]] && mkdir -p "tmp/download"
        [[ ! -f  tmp/gdc_download.out ]] && touch tmp/gdc_download.out

        FILE_NAME=${2}
        GDC_ID=${1}
        TCGA_ID=${3}

        echo "start Downloading WXS BAM files"
        gdc-client download -t $GDC_TOKEN -d tmp/download $GDC_ID

        wait

        # rename file

        file_name=($(echo ${FILE_NAME} | tr "."  ' '))
        mv tmp/download/${GDC_ID}/${file_name[0]}.bam download/${GDC_ID}.bam
        mv tmp/download/${GDC_ID}/${file_name[0]}.bai download/${GDC_ID}.bai

        # clean up
        rm -r tmp/download/${GDC_ID}

        echo -e "${GDC_ID}\tdownload/${file_name}.bam\t${TCGA_ID}" >> tmp/gdc_download.out

        echo "Done Downloading WXS BAM files"

}

############### main ##############

# [1]: GDC sample sheet

# Check the number of arguments
#if [ $# -eq 0 ]; then
#       echo "No arguments provided."
#       echo "gdc_download.sh {gdc_sample_sheet.tsv}"
#       exit
#fi



# checking file existence
[[ ! -f $1 ]] && {
        echo "$sample does not exist"
        exit
}

# checking bash_script running and lock gdc_download.sh
if [ ! -f tmp/gdc_download.lock ]
then echo -e "GDC download script is running on $(date +'%Y%m%d')" >> tmp/gdc_download.lock
else
        less tmp/gdc_download.lock | echo
        exit;
fi
# reading sample sheet
tail -n +2 $1 | while IFS=$'\t' read -r -a sample
do
        # [0]: file ID
        # [1]: file name
        # [6]: sample ID
        gdc_download ${sample[0]} ${sample[1]} ${sample[6]}
done

rm tmp/gdc_download.lock

