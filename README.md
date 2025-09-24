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


