import pandas as pd
import numpy as np
import os
import io
import re
import argparse
from natsort import natsort_keygen
from template.banner import wanglab_banner

VERSION: str="0.1.1"
THRESHOLD=5
VCF_COLNAME: list[str] = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']
TIDY_COLNAME: list[str] = ["CHROM", "Start Position", "ALT", "End Position", "Length", "SV_Type", "HOMSeq", "Read_Counts"]

def path(*args : str) -> str:
    return "/".join(args)



def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="pindel_filtering.py",description="using pindel SI and TD result to identify ITD candidate")
    parser.add_argument('-s', '--SampleSheet', required=True, help="Tumor Sample Sheet")
    parser.add_argument('-g', '--genomonITD_dir', required=True, help="path to the directory of genomon-ITD result")
    parser.add_argument('-p', '--pindel_dir', required=True, help="path to the directory of genomon-ITD result")
    parser.add_argument('-o', '--output_dir', required=False, default="./", help="the output directory, default is the running directory")
    parser.add_argument('-V', '--Version', action='version', version=VERSION, help="print out the Version of this program and exit")
    parser.add_argument('-v', '--verbose', action='store_true', required=False, default=False, help="turn on the verbosity of this program")
    return parser

def initiatePindel(data: pd.DataFrame) -> pd.DataFrame:
    data["Pindel"] = 0
    data["Pindel_SV_type"] = None 
    data["Pindel_chr"] = None
    data["Pindel_s_pos"] = np.nan
    data["Pindel_e_pos"] = np.nan
    data["min_diff"] = np.nan
    data["Sequence"] = None
    data["Pindel_read_counts"] = np.nan
    data["Pindel_length"] = np.nan

    return data

class PindelObject():
    def __init__(self, Pindel_dir, file_name) -> None:
        self.path = Pindel_dir
        self.file_name = file_name
        self.set_path()
        self.read()
        self.rawDict = self.splitPindelDict(self.raw)
        self.tidyDict = self.splitPindelDict(self.tidy)
        self.vcf_output = pd.DataFrame(columns= VCF_COLNAME)

    def read(self) -> None:
        self.set_vcf_header()
        self.raw = pd.concat([self.read_vcf(self.SI_path), self.read_vcf(self.TD_path)]).reset_index(drop=True)
        self.tidy = self.vcf_tidyup(self.raw)

    def splitPindelDict(self,data: pd.DataFrame) -> dict[pd.DataFrame]:
        thisDict = dict()
        for chromosome in set(data.CHROM):
            this_data = data[data.CHROM == chromosome]
            thisDict[chromosome]=this_data
        return thisDict

    def set_vcf_header(self) -> None:
        self.path_checkExistence(self.SI_path)
        with open(self.SI_path) as f:
            self.vcf_header = [l for l in f if  l.startswith('##')]

    def set_path(self) -> None:
        self.SI_path = path(self.path, self.file_name + "_SI.vcf") # open relative pindel SI vcf
        self.TD_path = path(self.path, self.file_name + "_TD.vcf") # open relative pindel TD vcf

    def path_checkExistence(self, path: str) -> None:
        assert (os.path.exists(path)), f"{path} is not exist"
        assert (os.path.getsize(path) != 0), f"{path} is empty"

    def read_vcf(self,file):
        self.path_checkExistence(file)

        data = pd.read_csv(file, sep='\t', comment='#', header=None)
        data.columns = VCF_COLNAME

        return data

    def vcf_tidyup(self, vcf: pd.DataFrame) -> None:
        def TransformINFO(information: str) -> list:
        
            variant_dict =dict(re.findall(r'(\w+)=(\w+)', information))
            keys_to_list = ["END", "SVLEN", "SVTYPE", "HOMSEQ"]

            return [variant_dict[key] for key in keys_to_list]

        def TransformCounts(CountsInfo: str) -> int:
            return int(re.search('(?<=:)-?\d+,-?\d+', CountsInfo).group().split(",")[1])

        tidy_vcf  = pd.DataFrame(vcf["INFO"].map(TransformINFO).to_list(), columns=TIDY_COLNAME[3:7]) # "End Position", "Length", "SV_Type", "HOMSeq"
        tidy_vcf = pd.DataFrame(vcf[["CHROM","POS","ALT"]]).merge(tidy_vcf, left_index=True, right_index=True) #CHROM", "Start Position", "ALT",
        tidy_vcf.loc[tidy_vcf['SV_Type'] == 'INS', 'End Position'] = tidy_vcf.loc[tidy_vcf['SV_Type'] == 'INS', 'POS'] + tidy_vcf.loc[tidy_vcf['SV_Type'] == 'INS', 'Length']
        tidy_vcf["Read_Counts"] = vcf["SAMPLE"].map(TransformCounts)
        tidy_vcf.columns = TIDY_COLNAME
        return tidy_vcf

    def vcf_append(self, row: pd.Series):
        self.vcf_output.loc[len(self.vcf_output.index)] = row

    def to_csv_mod(self, file_name) -> None:
        with open(file_name, "w") as f:
            for line in self.vcf_header:
                f.write(line)
            this_vcf = self.vcf_output.rename(columns={'CHROM': '#CHROM'})
            vcf_str = this_vcf.to_csv(sep='\t', index=False)
            f.write(vcf_str)
            f.close()
    

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    sample_sheet = args.SampleSheet
    tumor_dir = args.genomonITD_dir
    pindel_dir = args.pindel_dir
    output_dir = args.output_dir
    verbose = args.verbose
 

    if verbose:
        wanglab_banner()
        print(f"""
{'-'*20}
Sample Sheet            : {sample_sheet}
Tumor Sample Directory  : {tumor_dir}
Pindel Directory        : {pindel_dir}
Output Directory        : {output_dir}
{'-'*20}\n
""")
    
    samples = pd.read_table(sample_sheet, index_col = None)["File ID"].to_list()

    for sample in samples:
        g_data = initiatePindel(pd.read_csv(path(tumor_dir,sample+"_tidy.csv"), header=0, index_col=None))
        Pindel = PindelObject(pindel_dir, sample)

        for i, g_row in g_data.iterrows():
            chrom, s_pos, e_pos, seq = g_row[["chromosome1", "position1", "position2", "observed_inserted_nucleotide"]]

            if chrom in Pindel.tidyDict:
                thisChrPindel: pd.DataFrame = Pindel.tidyDict[chrom]
                thisChrRaw: pd.DataFrame = Pindel.rawDict[chrom]
            else:
                continue

            for j, p_row in thisChrPindel.iterrows():
                start_diff: int = abs(int(p_row["Start Position"])-int(s_pos))
                end_diff: int =  abs(int(p_row["End Position"])-int(e_pos))
                if start_diff <= THRESHOLD and end_diff <= THRESHOLD:
                    g_data.loc[i,"Pindel"] = 1
                    g_data.loc[i,"Pindel_SV_type"] = p_row["SV_Type"]
                    g_data.loc[i,"Pindel_chr"] = chrom
                    g_data.loc[i,"Pindel_s_pos"] = p_row["Start Position"]
                    g_data.loc[i,"Pindel_e_pos"] = p_row["End Position"]
                    g_data.loc[i,"min_diff"] = max(start_diff, end_diff)
                    g_data.loc[i,"Sequence"] = p_row["ALT"][-1::-1]
                    g_data.loc[i,"Pindel_read_counts"] = p_row["Read_Counts"]
                    g_data.loc[i,"Pindel_length"] = p_row["Length"]

                    Pindel.vcf_append(thisChrRaw.iloc[j])
                    break
        
        file_name = path(output_dir, sample+".pindel.csv")
        vcf_name = path(output_dir, sample+".pindel.vcf")

        g_data.to_csv(file_name, index=None)
        Pindel.to_csv_mod(vcf_name)
