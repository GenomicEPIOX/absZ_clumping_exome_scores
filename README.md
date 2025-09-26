# absZ_clumping_exome_scores
Code used to create PRS panels for the exome-scores and regression scripts

Need to create a LD matrix on the RAP - Extract all pQTLs from the genotypic model from Dhindsa et al   eg..
~~~bash
./plink2 --bfile Dhindsa_plink1 --r2-phased--ld-window-kb 10000 --ld-window-r2 0.01 --out Dhindsa_LD_matrix
~~~

Turn this into a HD5 using make_ld_matrix_hd5.py eg..
~~~bash
python3 make_ld_matrix_hd5.py -i pQLT_LD_matrix_9_1_24.vcor
~~~

Create panel using panel_generator.py 


Need to update this 

~~~bash
or i in $(cat cis/cis_protiens.txt); do python run_regressions.py -p ${i} -o cis_results/ -c prostate_inc -pheno pheno/ -cov UKBB_covar_witharray_relative_kept.txt -panel cis/ -v cis -r2 0.01 --plink cis_results/${i}/${i}_cis.sscore --sex 1; done
~~~

~~~bash
for i in $(cat cis/cis_protiens.txt); do python run_regressions.py -p ${i} -o cis_results/ -c prostate_inc -pheno pheno/ -cov UKBB_covar_witharray_relative_kept.txt -panel cis/ -v cis -r2 0.01 --plink cis_results/${i}/${i}_cis.sscore --sex 1; done
~~~