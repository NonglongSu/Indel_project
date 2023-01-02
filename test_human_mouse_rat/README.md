# Profliling of Indel phases in coding regions
* Find the **real indels ** that cannot be determined due to the limitation of current software.  
* Use indel phases  

## Requirements
* R 3.4.4 + 
* Bioconductor 3.6 +
* Ensembl Genes 100 (Human Genes GRCh38.p13)
* Mafft --version V.7.407 
* Prank --version v.170427  
 
## Workflow
## Part I
### 1. Find homologs in multiple species (Ex. human/mouse/rat), and extract cds  
make homoCall  
make cds  
* human/mouse/rat            -- 14547
* human/hamster/mouse/rat    -- 13890
* mouse/human/macaque        -- 14422

### 2. Quality control test
make test_mul_3  
make test_N  
make test_stop  
make test_preStop  
   
make filter_cds   
* human/mouse/rat            -- 13827  
* human/hamster/mouse/rat    -- 13881
* mouse/human/macaque        -- 14418


### 3. Translation of nucleotide to amino acid for each sequence.
make aa  

### 4. Multiple sequence alignment of cds.
make mafft  
make prank  

### 5. Abnormal mapping (From aa --> dna)
make mapped_mafft  
make mapped_prank  
make pair_mafft
make pair_prank


### 6. Quality control plot
make mafft_sim   
make prank_sim  


## Part II
### 7. Include gaps of length = {3, 6, 9, 12} only.
make mafft_upc   
make prank_upc  

### 8. Find the "true" gaps based on sliding-window method. 
make mafft_sw  
make prank_sw  

* h/m/r >>>>
* Data_3   : 2034
* Data_6   : 1712 
* Data_9   : 1652
* Data_12  : 1605

* h/c/m/r >>>>
* Datat_6.2 : 1101   


## Figure
### 1. Count the proportion of phase 0, 1, 2 indels of focal species together.
make mafft_phase  
make prank_phase  

### 2. Generate displacement bar plot between prior-post indels.  
make mafft_dis  
make prank_dis  

* test_reverse proves that it is due to the mafft alignment bias.  

### 3. Generate edge-bias HTML report
make mafft_edge  
make prank_edge   
make mafft_edge_report  
make prank_edge_report  

### 4. Count the proportion of phase 0, 1, 2 indels of focal species seperately.  
make mafft_phase_sep  
make prank_phase_sep    

### 5. Count the effective proportion of phase 1 & 2 indels.   
make mafft_phase_eff   
make prank_phase_eff   

### 6. Find the common ancestor based on maximum parsimony method.  
make mafft_anc  
make prank_anc  

### 7. Calculate the genome-wide dnds of focal species
make dnds 

### 8. Calculate the geomewide znzs of focal species
make znzs     

## In DataB_6  
###1 Find the best phase and update inside the folder
make mafft_mix

###2. Find the correlation between dnds and znzs (not working)
make omega_ladder  
make w_zed_corr     

###3. Simulate the znzs ratio within 6-mers
make sim_zedv







