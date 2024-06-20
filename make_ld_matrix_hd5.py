import pandas as pd
import argparse
import sys
import warnings
warnings.filterwarnings('ignore')

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Creates an LD matrix file in hd5 format to speed things up a little ' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help='The name of the vcor ld file that plink2 generates',required='True')
    parser.add_argument('-o', '--output', help='output name default = Dhindsa_ld_matrix.hd5 ',default="Dhindsa_ld_matrix.hd5")
    results = parser.parse_args(args)
    return ( results.input, results.output) 

def main(input_name, output): 
    ld_matrix = pd.read_csv(input_file, sep="\t", low_memory=False)
    ld_matrix.rename(columns={"#CHROM_A" : "CHR"} , inplace=True)
    for i in range(1,23): 
        temp = ld_matrix[ld_matrix["CHR"] == str(i) ]
        temp.head(5)
        #temp["key"] = temp.ID_A + "_" + temp.ID_B 
        key_name = "chr" + str(i) 
        temp.to_hdf(output, key=key_name, mode='a')
    # temp = ld_matrix[ld_matrix["CHR"] == "Y" ]
    # #temp["key"] = te mp.ID_A + "_" + temp.ID_B 
    # key_name = "chrY" 
    # temp.to_hdf(output, key=key_name, mode='a')
    temp = ld_matrix[ld_matrix["CHR"] == "X" ]
    temp['CHR'] = temp['CHR'].replace({'X': 23 })

    #temp["key"] = temp.ID_A + "_" + temp.ID_B 
    key_name = "chr23" 
    temp.to_hdf(output, key=key_name, mode='a')
    


if __name__ == '__main__':
    input_file, output,  =check_arg(sys.argv[1:])