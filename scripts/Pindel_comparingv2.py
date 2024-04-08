import sys
import pandas as pd
import os
import argparse
import numpy as np
import io
from natsort import natsort_keygen

def Compare(Genomon_path, Pindel_path, fileID, sampleID, OutputDir):
    # Read the result file from GenomonITDetector
    GenomonITD = readITD(Genomon_path, fileID)
    Pindel_tidy = readPindel(Pindel_path, fileID, sampleID )
    PindelDict = splitPindelDict(Pindel_tidy)

    # check each row in GenomonITD results
    for i in range(GenomonITD.shape[0]): # for each ITD 
        sample_id = GenomonITD.iloc[i,6]
        chr = GenomonITD.iloc[i,1]
        s_pos = GenomonITD.iloc[i,2]
        e_pos = GenomonITD.iloc[i,4]
        seq = GenomonITD.iloc[i,12]
        
        if chr in PindelDict:
            thisChrPindel = PindelDict[chr]
        else:
            continue

        # Start to compare
        for j in range(thisChrPindel.shape[0]):
            if (chr == thisChrPindel.iloc[j,0]) and ((abs(s_pos-thisChrPindel.iloc[j,1]) <= 5) and (abs(e_pos-thisChrPindel.iloc[j,3]) <= 5)):
                #GenomonITD.iat[i,14] = 1
                GenomonITD.loc[i,"Pindel"] = 1
                GenomonITD.loc[i,"Pindel_SV_type"] = thisChrPindel.iloc[j,5]
                GenomonITD.loc[i,"Pindel_chr"] = thisChrPindel.iloc[j,0]
                GenomonITD.loc[i,"Pindel_s_pos"] = thisChrPindel.iloc[j,1]
                GenomonITD.loc[i,"Pindel_e_pos"] = thisChrPindel.iloc[j,3]
                GenomonITD.loc[i,"min_diff"] = min(abs(s_pos - thisChrPindel.iloc[j,1]),abs(e_pos - thisChrPindel.iloc[j,3]))
                GenomonITD.loc[i,"Pindel_read_counts"] = thisChrPindel.iloc[j,7]
                GenomonITD.loc[i,"Pindel_length"] = thisChrPindel.iloc[j,4]
                lengths = (-1) * int(thisChrPindel.iloc[j,4])
                GenomonITD.loc[i,"Sequence"] = thisChrPindel.iloc[j,2][lengths:]
                break
    
    file_name = OutputDir + fileID + ".pindel.csv"
    GenomonITD.to_csv(file_name, index=False)

def readITD(Genomon_path, fileID):
    GenomonITD = pd.read_csv(Genomon_path + fileID + '.tidy.csv', header=0, index_col=None)
    GenomonITD["Pindel"] = 0
    GenomonITD["Pindel_SV_type"] = None
    GenomonITD["Pindel_chr"] = None
    GenomonITD["Pindel_s_pos"] = np.nan
    GenomonITD["Pindel_e_pos"] = np.nan
    GenomonITD["min_diff"] = np.nan
    GenomonITD["Sequence"] = None
    GenomonITD["Pindel_read_counts"] = np.nan
    GenomonITD["Pindel_length"] = np.nan

    return GenomonITD

def splitPindelDict(Pindel_tidy):
    PindelDict = dict()
    for chromosome in set(Pindel_tidy.CHROM):
           this_data = Pindel_tidy[Pindel_tidy.CHROM == chromosome]
           PindelDict[chromosome]=this_data
    return PindelDict

def readPindel(Pindel_path, fileID, sampleID):
    Pindel_tidy = pd.DataFrame(columns=["CHROM", "Start Position", "ALT", "End Position", "Length", "SV Type", "HOMSeq", "Read_Counts"])
    pindel_SI_path = Pindel_path + fileID + "_SI.vcf" # open relative pindel SI vcf 
    pindel_TD_path = Pindel_path + fileID + "_TD.vcf" # open relative pindel TD vcf 
    if (os.path.isfile(pindel_SI_path)) and (os.path.getsize(pindel_SI_path) != 0):
        Pindel_tidy = pd.concat([Pindel_tidy, TidyPindel(pindel_SI_path, sampleID)])
    if (os.path.isfile(pindel_TD_path)) and (os.path.getsize(pindel_TD_path) != 0):
        Pindel_tidy = pd.concat([Pindel_tidy, TidyPindel(pindel_TD_path, sampleID)])
    Pindel_tidy = Pindel_tidy.sort_values(by="CHROM",key=natsort_keygen(), ascending=True)
    return Pindel_tidy

def TidyPindel(pindel_file,sampleID):
    # Tidy the information of Pindel result
    Pindel = read_vcf(pindel_file)
    Pindel = Pindel[["CHROM","POS", "ALT", "INFO", sampleID]]
    all_all_data = []
    for k in range(Pindel.shape[0]):
        keep_data = [Pindel.iloc[k,0], int(Pindel.iloc[k,1]), Pindel.iloc[k,2]]
        to_add = TransformINFO(Pindel.iloc[k,3])
        to_add_2 = TransformCounts(Pindel.iloc[k,4])
        all_info = keep_data + to_add + to_add_2
        if all_info[5] == "INS":
            all_info[3] = all_info[1] + all_info[4]
        all_all_data.append(all_info)
    Pindel_tidy = pd.DataFrame(all_all_data, columns=["CHROM", "Start Position", "ALT", "End Position", "Length", "SV Type", "HOMSeq", "Read_Counts"])
    return Pindel_tidy

def read_vcf(files):
    with open(files, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})

def TransformINFO(infomation):
    info_split = []
    first_split = infomation.split(";")
    for i in range(len(first_split)):
        info_split.append(first_split[i].split("="))
    
    to_add = []
    HOMSeq = "NaN"
    for j in range(len(info_split)):
        category = info_split[j][0]
        
        if category == "END":
            end_position = int(info_split[j][1])
            #print("End Position: ",end_position)
        elif category =="SVLEN":
            length = int(info_split[j][1])
            #print("length: ",length)
        elif category =="SVTYPE":
            SVTYPE = info_split[j][1]
            #print("SV type: ", SVTYPE)
        elif category == "HOMSEQ":
            HOMSeq = info_split[j][1]

    to_add = [end_position, length, SVTYPE, HOMSeq]

    return to_add

def TransformCounts(CountsInfo):
    first = CountsInfo.split(":")
    to_second = first[1].split(",")
    return [to_second[1]]


#getArguments from execution
def get_parser():
    parser = argparse.ArgumentParser(prog="Compare_between_Tools.py",description="Compare the results from Genomon-ITDetector and Pindel")
    parser.add_argument('-g', '--GenomonPath', required=True, help='The path of Genomon-ITDetector results')
    parser.add_argument('-p', '--PindelPath', required=True, help='The path of the Pindel files')
    parser.add_argument('-f', '--FileID', required=True, help='The File ID ')
    parser.add_argument('-s', '--SampleID', required=True, help='The Sample ID')
    parser.add_argument('-o', '--OutputDir', required=False, help="the Output Directory", default='./')
    return parser


# main
def main():
    parser = get_parser()
    args = parser.parse_args()

    GenomonITDFile = args.GenomonPath
    PindelFile = args.PindelPath
    FileID = args.FileID
    SampleID = args.SampleID
    OutputDir = args.OutputDir

    # 確定可以get到參數
    print("Check for the arguments!")
    print("----------")
    print("Genomon-ITDetector results:")
    print(GenomonITDFile)
    print("Pindel results:")
    print(PindelFile)
    Compare(GenomonITDFile, PindelFile, FileID, SampleID, OutputDir)


if __name__ == '__main__':
    main()
