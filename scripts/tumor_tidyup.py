import pandas as pd
import os
import argparse
import re
from natsort import natsort_keygen
from template.banner import wanglab_banner

VERSION: str="0.1.0"

POSITION_PATTERN = "([a-zA-Z0-9_]+):\\+([0-9]+)-([a-zA-Z0-9_]+):\\-([a-zA-Z0-9]+)"
COL_NAME: list[str] = ["tumor_type", "chromosome1", "position1", "chromosome2", "position2", 
                       "sample_ID", "max_read_counts", "grade", 
                       "RefSeq_gene_Exon", "RefSeq_gene_Intron", "Ensembl_gene",
                       "assembled_contig_sequence", "observed_inserted_nucleotide",
                       "OIN_length", "PDN_length", "in_normal", "in_pindel"]
DATA_ORDER: list[int] = [34,22,23,28,10,15,16,17]

def path(*args : str) -> str:
    return "/".join(args)

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="Normal_Pooling.py",description="Pooling Normal ITD data")
    parser.add_argument('-t', '--TumorType', required=True, help='Tumor Type')
    parser.add_argument('-s', '--SampleSheet', required=True, help="Normal Sample Sheet")
    parser.add_argument('-g', '--genomonITD_dir', required=True, help="path to the directory of genomon-ITD result")
    parser.add_argument('-o', '--output_dir', required=False, default="./", help="the output directory, default is the running directory")
    parser.add_argument('-V', '--Version', action='version', version=VERSION, help="print out the Version of this program and exit")
    parser.add_argument('-v', '--verbose', action='store_true', required=False, default=False, help="turn on the verbosity of this program")
    return parser

def extract_position(data: pd.DataFrame) -> pd.DataFrame:
    positions = data.iloc[:,0].map(apply_position_regex).to_list() 
    return pd.DataFrame(positions, columns=COL_NAME[1:5])

def apply_position_regex(string: str) -> list[str]:
    result = re.match(POSITION_PATTERN,string) # pattern: chr2, pos2, chr1 , pos1
    return result.group(3), result.group(4), result.group(1), result.group(2) 

def max_read_counts(row) -> int:
    return int(max(sum([row[1], row[2]]), sum([row[4], row[5]])))

def get_grade(data: pd.DataFrame) -> pd.DataFrame:
    pass

if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    tumor_type = args.TumorType
    sample_sheet = args.SampleSheet
    tumor_dir = args.genomonITD_dir
    output_dir = args.output_dir
    verbose = args.verbose
 

    if verbose:
        wanglab_banner()
        print(f"""
{'-'*20}
Tumor Type              : {tumor_type}
Sample Sheet            : {sample_sheet}
Tumor Sample Directory : {tumor_dir}
Output Directory        : {output_dir}
{'-'*20}\n
""")
        
    tumor_samples = list(pd.read_csv(sample_sheet, index_col = None)["File ID"])
    

    for sample in tumor_samples:
        this_file = path(tumor_dir,sample+"_itd_list.tsv")
        
        if (os.path.exists(this_file)) and (os.path.getsize(this_file) != 0):
            this_sample = pd.read_table(this_file, header=None, index_col=None)
        else:
            print(f"{this_file} is not exists or has size of zero !!!")
            continue
    
        this_data = extract_position(this_sample)   
        this_data["tumor_type"] = tumor_type
        this_data = this_data[COL_NAME[:5]]
        this_data["sample_ID"] = sample
        this_data["max_read_counts"]= this_sample.apply(lambda row: max_read_counts(row), axis=1)
        other_data = pd.DataFrame(this_data[:,DATA_ORDER], columns=COL_NAME[7:15])
        this_data = pd.concat([this_data, other_data], ignore_index=True)
        this_data["in_normal"] = 0
        this_data["in_pindel"] = 0

        this_data = this_data.sort_values("chromosome1", key=natsort_keygen(), ascending=True)
        this_data.to_csv(f"{output_dir}/{sample}_tidy.csv", index=None)


    

    
