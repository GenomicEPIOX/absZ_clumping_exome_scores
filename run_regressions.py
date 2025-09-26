import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import argparse
import sys
import os
import glob
from subprocess import Popen, PIPE
import numpy as np
from sklearn.preprocessing import StandardScaler
import statsmodels.api as sm
import statsmodels.formula.api as smf
import shlex
import re
scaler = StandardScaler()
from subprocess import Popen, PIPE


def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Creates the exome-score panels, selecting variant with the largest ABS_Z in LD > 0.01 ' ,formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-p', '--protein_name', help='Protein name',required='True')
    parser.add_argument('-o', '--output', help='the path where the results will be stored, default= "./results/" ',default="./results/")
    parser.add_argument('-plink', '--plink_file', help='The plink files to run the exome-score on',required='True')
    parser.add_argument('-sex', '--sex_specific', help='If cancer site is sex-specific set this to 1, will remove sex from the model',default=0)
    parser.add_argument('-pheno', '--pheno_path', help='folder where cancer phenotype files are kept',required='True')
    parser.add_argument('-c', '--cancer', help='cancer site name to find pheno file plus add labels {pheno_path}{cancer}.pheno',required='True')
    parser.add_argument('-cov', '--covar_file', help='The full path including file name for the covariate file',required='True')
    parser.add_argument('-panel', '--panel_path', help='The path of where the panel files are stored',required='True')
    parser.add_argument('-v', '--variant_type', help='select the type of variants to create the panel from - cis, trans, combined or trans_cis_pos',required='True')
    parser.add_argument('-r2', '--r2', help='R-square threashold',required='True')
    results = parser.parse_args(args)
    return (results.protein_name, results.output, results.plink_file, results.sex_specific, results.pheno_path , results.cancer, results.covar_file, results.panel_path, results.variant_type, results.r2) 


def make_p(cancer, protein,covar_file, pheno_path, variant_type, r2,plink_file):
    print (f"Performing normal regression for {protein} on {cancer} " )
    cancer_file_name =f"{pheno_path}{cancer}.pheno"
    protein_scores = pd.read_csv(plink_file, sep="\t")
    cancer_pheno = pd.read_csv(cancer_file_name, sep="\t")
    cancer_pheno = cancer_pheno.dropna()
    merged = protein_scores.merge(cancer_pheno, on="IID")
    small = merged[["IID","SCORE1_SUM", cancer]].copy()
    covar = pd.read_csv(covar_file, sep="\t")
    dataset = small.merge(covar, on="IID")
    dataset['scaled_score'] = scaler.fit_transform(dataset[['SCORE1_SUM']])
    dataset["cc"] = dataset[cancer]
    formula = "cc ~ scaled_score +  SEX + age_recruitment + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 "
    results = do_stuff(dataset,formula, protein, cancer, variant_type, r2)
    return results
    
def make_p_sex(cancer, protein, covar_file, pheno_path, variant_type, r2, plink_file):
    print (f"Performing sex-specific regression (eg sex is left out of the model) for {protein} on {cancer} " )
    cancer_file_name =f"{pheno_path}{cancer}.pheno"
    protein_scores = pd.read_csv(plink_file, sep="\t")
    cancer_pheno = pd.read_csv(cancer_file_name, sep="\t")
    cancer_pheno = cancer_pheno.dropna()
    merged = protein_scores.merge(cancer_pheno, on="IID")
    small = merged[["IID","SCORE1_SUM", cancer]].copy()
    covar = pd.read_csv(covar_file, sep="\t")
    dataset = small.merge(covar, on="IID")
    dataset['scaled_score'] = scaler.fit_transform(dataset[['SCORE1_SUM']])
    dataset["cc"] = dataset[cancer]
    formula = "cc ~ scaled_score + age_recruitment + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 "
    results = do_stuff(dataset,formula, protein, cancer, variant_type, r2)
    return results    

def do_stuff(test_set, formula1, protein, cancer, variant_type, r2):
    n = test_set[test_set["cc"] == 1 ].shape[0]
    n_cont = test_set[test_set["cc"] == 0 ].shape[0]
    #formula1 = "cc ~ scaled_score +  SEX + age_recruitment + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 "
    logitfit = smf.logit(formula = formula1, data = test_set).fit(disp=0)
    params = logitfit.params
    conf = logitfit.conf_int()
    conf['Odds Ratio'] = params
    conf.columns = ['5%', '95%', 'Odds Ratio']
    moo = (round(np.exp(conf),3))
    OR = moo.iloc[1].values[2] 
    UPP = moo.iloc[1].values[1]
    LOW = moo.iloc[1].values[0]
    P = logitfit.pvalues[1]
    SE = logitfit.bse[1]
    Z = round(np.log(OR) / SE,2)
    log_results = pd.DataFrame( {'Protein': protein, 'Cancer': cancer, 'panel_type' : variant_type, 'clumping_R2': r2,
     'OR': OR,
     'lower': LOW, 
     'upper' : UPP ,
     'p' : P , 
     'SE' : SE ,
     'Z' : Z, 
    }, index=[0])
    ### converting to str for plots
    log_results["n"] = str(n)
    log_results["n_cont"] = str(n_cont)
    log_results["P_val"] = log_results["p"].apply(lambda x:  '%0.2e' % x) 
   
    return log_results

def execute(cmd):
    """ to run commands on the shell and pipe out the display """
    plink_cmd_run = shlex.split(cmd)
    p = Popen(plink_cmd_run, shell=False, stdout=PIPE , stderr=PIPE)
    out, err = p.communicate()
    return (out)

def run_plink(protein, plink_file, panel_path, variant_type, r2):
    plink_cmd = (f"./plink2 --pfile {plink_file} --score {panel_path}{protein}_{variant_type}_{r2}.panel 1 2 3 list-variants  ignore-dup-ids cols=+scoresums --out  {protein}" )
    p = execute(plink_cmd)

def main(protein_name, output, plink_file, sex, pheno_path , cancer, covar_file, panel_path, variant_type, r2 ):
    #run_plink(protein_name, plink_file, panel_path, variant_type, r2)

    if sex == str(1) : 
        temp_df = make_p_sex(cancer, protein_name, covar_file, pheno_path,variant_type, r2 ,plink_file )
    else: 
        temp_df = make_p(cancer, protein_name, covar_file, pheno_path, variant_type, r2 , plink_file)

    output_folder = f"{output}{protein_name}"
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    output_name = f"{output_folder}/{protein_name}.{cancer}.{variant_type}.{r2}.{sex}.txt"
    temp_df.to_csv(output_name, sep="\t", index=None, header=None)


if __name__ == '__main__':
    protein_name, output, plink_file, sex, pheno_path , cancer, covar_file, panel_path, variant_type, r2 =check_arg(sys.argv[1:])
    main (protein_name, output, plink_file, sex, pheno_path , cancer, covar_file, panel_path, variant_type, r2)