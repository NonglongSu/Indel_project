#!/bin/bash
# housekeeping % vs dnds correlation
# indel position distribution

default: all
.PHONY:  default all clean help
.PHONY:  rm_ref coati_best id_hexplot filter_coati
.PHONY:  dnds_znzs_sum hk_logit edge_bias_ecm plot_indel_pos


RSCRIPT = Rscript --vanilla
RM      = rm -i
SCR     = ../Script

###################################
#remove the reference
#rm_ref:$(SCR)/remove_ref.sh Raw_data/cds_seq Raw_data/pair_cds
#	bash $< $(word 2,$^) $(word 3,$^)
#
##coati alignment
#coati_best:$(SCR)/coati_best.sh Raw_data/pair_cds Raw_data/coati_align
#	bash $< $(word 2,$^) $(word 3,$^)
#
##alignment filtering
#id_hexplot:../Script/QC/indel_hex_plot.R Raw_data/coati_align
#	$(RSCRIPT) $< $(word 2,$^) Raw_data/QC/id.hex.txt Raw_data/Figure/id.hex.pdf
#
#filter_coati: ../Script/QC/filter_coati.sh Raw_data/QC/id.hex.txt
#	bash $< $(word 2,$^)
#
##Generate a table of sum stas: Zn,Zs,ZN,ZS,Nd,N,Sd,S,lenA,lenB
#dnds_znzs_sum:$(SCR)/dnds_znzs_sum.R Data_6/Mafft/mapped_cds
#	$(RSCRIPT) $< $(word 2,$^) Results/dnds_znzs_sum.tsv
#

#housekeeping gene regression
hk_logit:$(SCR)/housekeeping_gene_track.R housekeeping_genes_human.csv Results/dnds_znzs_sum.tsv Raw_data/cds_seq Data_6/Mafft/mapped_cds 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^)  $(word 5,$^) Results/hk_logit.pdf

#check if ecm related to the edge bias issue.
#Print "Percentage of original alignment in ECM" on terminal
edge_bias_ecm:$(SCR)/edge_ecm_cat.R
	$(RSCRIPT) $< 


#plot the indel position distribution
plot_indel_pos:$(SCR)/indel_pos.R  
	$(RSCRIPT) $< Results/Phase.mafft.pos.pdf




  
  

