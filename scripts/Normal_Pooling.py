import sys
import pandas as pd
import os
import argparse
import re
from natsort import natsort_keygen
from template.banner import wanglab_banner

VERSION: str="0.1.1"

POSITION_PATTERN = "([a-zA-Z0-9_]+):\\+([0-9]+)-([a-zA-Z0-9_]+):\\-([a-zA-Z0-9]+)"
COL_NAME: list[str] = ["chromosome1", "position1", "chromosome2", "position2", "max_read_counts", "sample_ID"]

def path(*args : str):
    return "/".join(args)

def get_parser():
    parser = argparse.ArgumentParser(prog="Normal_Pooling.py",description="Pooling Normal ITD data")
    parser.add_argument('-t', '--TumorType', required=True, help='Tumor Type')
    parser.add_argument('-s', '--SampleSheet', required=True, help="Normal Sample Sheet")
    parser.add_argument('-g', '--genomonITD_dir', required=True, help="path to the directory of genomon-ITD result")
    parser.add_argument('-o', '--output_dir', required=False, default="./", help="the output directory, default is the running directory")
    parser.add_argument('-V', '--Version', action='version', version=VERSION, help="print out the Version of this program and exit")
    parser.add_argument('-v', '--verbose', action='store_true', required=False, default=False, help="turn on the verbosity of this program")
    return parser

def extract_position(data: pd.DataFrame):
    positions = data.iloc[:,0].map(apply_position_regex).to_list() 
    return pd.DataFrame(positions, columns=COL_NAME[:4])

def apply_position_regex(string: str):
    result = re.match(POSITION_PATTERN,string) # pattern: chr2, pos2, chr1 , pos1
    return result.group(3), result.group(4), result.group(1), result.group(2) 

def max_read_counts(row):
    return int(max(sum([row[1], row[2]]), sum([row[4], row[5]])))

def remove_duplicate(data):
    print(COL_NAME[4])
    print(COL_NAME[5])
    data[COL_NAME[4]] = data.groupby(COL_NAME[:4])[COL_NAME[4]].transform(lambda x: ', '.join(map(str, x)))
    data[COL_NAME[5]] = data.groupby(COL_NAME[:4])[COL_NAME[5]].transform(lambda x: ', '.join(map(str, x)))

    return data.drop_duplicates(subset=COL_NAME[:4], keep='first')


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    tumor_type = args.TumorType
    sample_sheet = args.SampleSheet
    normal_dir = args.genomonITD_dir
    output_dir = args.output_dir
    verbose = args.verbose
 

    if verbose:
        wanglab_banner()
        print(f"""
{'-'*20}
Tumor Type              : {tumor_type}
Sample Sheet            : {sample_sheet}
Normal Sample Directory : {normal_dir}
Output Directory        : {output_dir}
{'-'*20}\n
""")
        
    normal_samples = list(pd.read_csv(sample_sheet, index_col = None)["File ID"])
    
    data = pd.DataFrame(columns=COL_NAME)

    for sample in normal_samples:
        this_file = path(normal_dir,sample+"_itd_list.tsv")
        
        if (os.path.exists(this_file)) and (os.path.getsize(this_file) != 0):
            this_sample = pd.read_table(this_file, header=None, index_col=None)
        else:
            print(f"{this_file} is not exists or has size of zero !!!")
            continue
        
        this_data = extract_position(this_sample)
        this_data["max_read_counts"]= this_sample.apply(lambda row: max_read_counts(row), axis=1)
        this_data["sample_ID"] = sample
        data = pd.concat([data,this_data], ignore_index=True)
    
    data = remove_duplicate(data.sort_values(by="chromosome1", key=natsort_keygen(), ascending=True)).reset_index(drop=True)

    for chromosome in set(data.chromosome1):
        data[data.chromosome1 == chromosome].to_csv(f"{output_dir}/{tumor_type}_Normal_Candidate_{chromosome}.csv")
    
    data.to_csv(f"{output_dir}/{tumor_type}_Normal_Candidate.csv")
