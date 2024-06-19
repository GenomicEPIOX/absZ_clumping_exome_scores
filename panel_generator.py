import pandas as pd
import argparse
import sys
import os
import glob
from subprocess import Popen, PIPE
import shlex
import re
import warnings
warnings.filterwarnings('ignore')

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Creates the exome-score panels, selecting variant with the largest ABS_Z in LD > 0.01 ' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-p', '--protein_name', help='Protein name',required='True')
    parser.add_argument('-o', '--output', help='the path where to create the results, will create a folder with protein_name, default= "./results/" ',default="./results/")
    parser.add_argument('-ld', '--ld_matrix', help='plink LD matrix in hd5 format with chrXX as the key',required='True')
    parser.add_argument('-r2', '--r2', help='R-square threashold',default=0.05)
    parser.add_argument('-ex', '--excel_file', help='Excel file from Dhindsa et al, Supplementary Table 2',required='True')
    parser.add_argument('-mhc', '--mhc_filter', help='To filter out MHC variants set to 1, will not filter cis MHC',default=0)
    parser.add_argument('-v', '--variant_type', help='select the type of variants to create the panel from - cis, trans, combined or trans_cis_pos',required='True')
    parser.add_argument('-m', '--missing', help='File with variants to be removed (missing from plink exome dataset). Single variant ID per a row, no header',required='True')
    results = parser.parse_args(args)
    return ( results.protein_name, results.output, results.ld_matrix, results.r2, results.excel_file, results.mhc_filter, results.variant_type, results.missing ) 


def make_temp_file(Dhindsa_excel_file, protein, variant_type, missing_variants) :
    missing = pd.read_csv(missing_variants, sep="\t", header=None)
    missing.columns = ["Genotype"]
    raw = pd.read_excel(Dhindsa_excel_file) 
    small = raw[raw["model"] == "genotypic"]
    small = small[["Protein", "Genotype" , "beta", "se", "Match(?Cis)" , "cis_trans_position_1mb_coding_region" ]].copy()
    small = small[small["Protein"] == protein ]
    small = small[~small.Genotype.isin(missing.Genotype)]
    if small.empty : 
                print (f"Unable to find {protein}" ) 
                return None 
    small["Genotype"] = small.Genotype.str.replace("-",":" )
    small["Genotype"] = small.Genotype.str.replace("X:","23:" )
    small["CHR"] = small.Genotype.str.split(":",expand=True)[0]
    small["BP"] = small.Genotype.str.split(":",expand=True)[1]
    small["BP"] = small.BP.astype(int) 
    small["REF"] = small.Genotype.str.split(":",expand=True)[2]
    small["ALT"] = small.Genotype.str.split(":",expand=True)[3]
    small["Z"] = small.beta / small.se
    small["abs_Z"] = abs(small.Z)
    cis_only = small[small["Match(?Cis)"] == "Yes" ] 
    trans_only = small[small["Match(?Cis)"] == "No" ] 
    trans_cis_pos = small[small["cis_trans_position_1mb_coding_region"] == "cis-position, trans-gene"] 
    if variant_type == "cis" : 
        cis_only = small[small["Match(?Cis)"] == "Yes" ] 
        if cis_only.empty :
            print (f"No cis variants for protein : {protein}" ) 
            return cis_only
        return cis_only
    elif variant_type == "trans" :
        trans_only = small[small["Match(?Cis)"] == "No" ] 
        if trans_only.empty :
                print (f"No trans variants for protein : {protein}" ) 
                return trans_only
        return trans_only
    elif variant_type == "combined" :
        if small.empty: 
                print (f"No genetic variants at all for protein : {protein}" ) 
                return small
        return small

    elif variant_type == "trans_cis_pos" : 
        trans_cis_pos = small[small["cis_trans_position_1mb_coding_region"] == "cis-position, trans-gene"]
        if trans_cis_pos.empty:
                print (f"No cis-position, trans-gene variants for protein : {protein}" ) 
                return None
        return trans_cis_pos

    else : 
        print ("You have not selected cis, trans, combined or trans_cis_pos - check your input for variant type") 
        return None


def process_pair(df, i, m, ld_matrix, r2):
    temp_df = pd.DataFrame(columns=["SNP1" , "SNP2" , "DISTANCE","R2", "NOTE", "FLAG" ]) 
    SNP_dict = { "SNP1" : i , "SNP2" : m }
    print (SNP_dict) 
    result = None
    if i == m : 
        SNP_dict["NOTE"] = "Same SNP" 
        SNP_dict["FLAG"] = "OK" 
        SNP_DF = pd.DataFrame([SNP_dict])
        temp_df = pd.concat([temp_df , SNP_DF],  ignore_index=True)
        return temp_df, None
#try:
    subset_i = df[df["Genotype"] == i]
    subset_m = df[df["Genotype"] == m]
    
    if subset_i.empty or subset_m.empty: 
        SNP_dict["NOTE"] = "Dataframe for one of the SNPs was empty" 
        SNP_dict["FLAG"] = "ERROR" 
        SNP_DF = pd.DataFrame([SNP_dict])
        temp_df = pd.concat([temp_df , SNP_DF],  ignore_index=True)
        return temp_df, None

    bp1 = subset_i.BP.values[0]
    bp2 = subset_m.BP.values[0]
    distance = abs(bp1 - bp2)
    SNP_dict["DISTANCE"] = distance 

    if distance > 10000000:
        SNP_dict["NOTE"] = "SNPs are not nearby" 
        SNP_dict["FLAG"] = "ok" 
        SNP_DF = pd.DataFrame([SNP_dict])
        temp_df = pd.concat([temp_df , SNP_DF],  ignore_index=True)
        return temp_df, None
        

    if 0 < distance <= 10000000:
        try: 
            ld_info = float(ld_matrix.loc[(ld_matrix["ID_A"] == i) & (ld_matrix["ID_B"] == m), "PHASED_R2"].values[0] ) 
            #ld_info = ld_matrix[(ld_matrix["ID_A"] == i) &(ld_matrix["ID_B"] == m) ].PHASED_R2.values[0]
        except:
            try:
                ld_info = float(ld_matrix.loc[(ld_matrix["ID_A"] == m) & (ld_matrix["ID_B"] == i), "PHASED_R2"].values[0] ) 
                #ld_info = ld_matrix[(ld_matrix["ID_A"] == m) &(ld_matrix["ID_B"] == i) ].PHASED_R2.values[0]
            except: 
                #### add code if empty from LD matrix drop that SNP
                if ld_matrix[(ld_matrix["ID_A"] == m)  ].empty: 
                    SNP_dict["NOTE"] =  (f"{m} is missing from the LD matrix - Keeping {i} " )
                    SNP_dict["FLAG"] = "failed" 
                    result = i
                    SNP_DF = pd.DataFrame([SNP_dict])
                    temp_df = pd.concat([temp_df , SNP_DF],  ignore_index=True)
                    return temp_df, result
                
                if ld_matrix[(ld_matrix["ID_A"] == i) ].empty: 
                    SNP_dict["NOTE"] =  (f"{i} is missing from the LD matrix - Keeping {m} " )
                    SNP_dict["FLAG"] = "failed" 
                    result = m
                    SNP_DF = pd.DataFrame([SNP_dict])
                    temp_df = pd.concat([temp_df , SNP_DF],  ignore_index=True)
                    return temp_df, result

        SNP_dict["R2"] = ld_info
    
        if ld_info >= r2:
            i_absz = subset_i.abs_Z.values[0]
            m_absz = subset_m.abs_Z.values[0]
            if m_absz > i_absz:
                result = i
                SNP_dict["NOTE"] =  (f"{m} wins with abs z {m_absz} compared to {i} with {i_absz} " ) 
            elif m_absz < i_absz:
                result = m
                SNP_dict["NOTE"] = (f"{i} wins with abs z {i_absz} compared to {m} with {m_absz} " ) 

        if ld_info < r2 : 
            SNP_dict["NOTE"] = (f"{i} and {m} are independent " )
            result = None
            
    
    SNP_DF = pd.DataFrame([SNP_dict])
    temp_df = pd.concat([temp_df , SNP_DF],  ignore_index=True)
    return temp_df, result

def process_chromosome(chromosome, ld_file, df, r2): 
    current_chr_df = df[df["CHR"] == chromosome] 
    chr_key = "chr" + str(chromosome) 
    ld_matrix = pd.read_hdf(ld_file, key=chr_key)
    list1 = current_chr_df.Genotype.unique()

    remove_ld_snps = []
    log = pd.DataFrame(columns=["SNP1" , "SNP2" , "DISTANCE","R2", "NOTE", "FLAG" ]) 
    for i in list1:
        for m in list1:
            log_from_test, result = process_pair(current_chr_df, i, m, ld_matrix, r2)
            log = pd.concat([log , log_from_test],  ignore_index=True)
            if result is not None:
                remove_ld_snps.append(result)

    result = list(set(remove_ld_snps)) ##this removes duplicated snps in remove_ld_snps list 
    keep = [i for i in list1 if i not in result]

    return log , result, keep

def process_protein(protein_name, r2,ld_file, variant_type, filter_MHC, file_path, PROTEIN_DF):
    ###filter MHC is 0 or 1 
    ###variant type matches how the directory is created. 
    ### file path is the name of the directory where the filters files live - file_path/ where the variant type folers sit eg 
    ## file_path/cis/ file_path/trans etc

    #PROTEIN = pd.read_csv(file_path + variant_type + "/" + protein_name + ".txt", sep="\t") 
    
    if PROTEIN_DF.empty:
        print (f"{protein_name} contains no {variant_type} hits" ) 
        return None
    PROTEIN_DF = PROTEIN_DF.sort_values(["CHR", "BP"]) 
    PROTEIN_DF['CHR'] = PROTEIN_DF['CHR'].replace({'chr23': 'chrX', 23 : "X" })
    current_chr = PROTEIN_DF.CHR.unique() 
    results = {}
    main_log = pd.DataFrame(columns=["SNP1" , "SNP2" , "DISTANCE","R2", "NOTE", "FLAG" ]) 
    for chromosome in current_chr:
        if chromosome == 6 and filter_MHC == 1 :
            PROTEIN_DF = filter_MHC(PROTEIN_DF)
        log, result, keep = process_chromosome(chromosome, ld_file, PROTEIN_DF, r2)
        main_log = pd.concat([main_log , log],  ignore_index=True)
        
        results[chromosome] = keep
    flattened_list = [item for sublist in results.values() for item in sublist]
    KEEP_DF = PROTEIN_DF[PROTEIN_DF.Genotype.isin(flattened_list)] 
    score_file = make_score(KEEP_DF)
    folder_name = file_path + protein_name + "/" 
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    out_name = f"{folder_name}{protein_name}_{variant_type}_{r2}.panel" 
    log_outname = f"{folder_name}{protein_name}_{variant_type}_{r2}.log" 
    main_log.to_csv(log_outname, sep="\t", index=None, header=None)
    if not score_file.empty:
        score_file.to_csv(out_name, sep="\t", index=None, header=None)

def make_score(df): 
    small = df[["Genotype", "REF", "ALT", "beta"]].copy()
    small_bad = small[small["beta"] < 0 ] 
    small_good = small[small["beta"] > 0 ] 
    small_bad["weight"] = small_bad["beta"] * -1 
    small_bad["A1"] = small_bad["REF"] 
    small_bad = small_bad[["Genotype", "A1", "weight"]].copy()
    small_good["weight"] = small_good["beta"]
    small_good["A1"] = small_good["ALT"]
    small_good = small_good[["Genotype", "A1", "weight"]].copy()
    SCORE = pd.concat([small_good, small_bad], axis=0 ) 
    SCORE["Genotype"] = SCORE["Genotype"].str.replace("-",":")
    return(SCORE)

def filter_MHC(df): 
    ##chr6:28,510,120-33,480,577
    df = df[(df["CHR"] == 6) & (df["BP"] < 28510120) & (df["BP"] > 33480577) ] 
    return df

def main(protein_name, output, ld_matrix, r2, excel_file, mhc_filter, variant_type, missing_variants):
     r2 = float(r2) ## its parsing it in as a string
     variant_df = make_temp_file(excel_file, protein_name, variant_type, missing_variants)
     process_protein(protein_name, r2, ld_matrix, variant_type, mhc_filter , output, variant_df)

if __name__ == '__main__':
    protein_name, output, ld_matrix, r2, excel_file, mhc_filter, variant_type, missing_variants  =check_arg(sys.argv[1:])
    main (protein_name, output, ld_matrix, r2, excel_file, mhc_filter, variant_type, missing_variants)