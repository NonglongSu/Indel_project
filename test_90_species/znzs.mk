#!/bin/bash
default: all

.PHONY:  default all clean
.PHONY:  rm_gap qMark test_cds filter_cds filter_more neutral_ZnZs GC_count mle_90 plot_phase_prop plot_phase_dist plot_dnds_ZnZs plot_ZnZs_prop plot_omega_dnds
.PHONY:  aa align_mafft mapBack jcdis_sw id_hex_sw filter_sw update_gap sw_gap  phase_sw_90 ZnZs_dnds_sw roll_ZnZs_sw
.PHONY:  coati_max jcdis_max id_hex_max filter_coati_max phase_coati_max_90 ZnZs_dnds_max roll_ZnZs_max roll_ZnZs_max2
.PHONY:  coati_sampling filter_coati_sam phase_coati_sample_90 phase_sample_ensembl plot_phase_sample ZnZs_dnds_sample roll_ZnZs_sam
.PHONY: supp

RSCRIPT  = Rscript --vanilla
MAFFT    = mafft --maxiterate 1000 --globalpair
RM       = rm -i
SCR      = ../Script/90
WINDOW   = 6
WALL     = 12
namelist:= $(shell test -f Species.txt && cat Species.txt)
pat1     = sw
pat2     = max
pat3     = sample


################################################Data quality control
#Remove all gaps: "copy raw_cds/ to cds/"
rm_gap: $(SCR)/gap_rm.sh Raw_data/cds
	bash $< $(word 2,$^)

#Record fasta files with pre-stop codons.
test_cds: $(SCR)/test_cds_90.R Raw_data/cds Raw_data/QC/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^)

#remove all disqualified files "include aa.fasta"
filter_cds: $(SCR)/filter_cds.sh Raw_data/QC Raw_data/cds
	bash $< $(word 2,$^) $(word 3,$^)

#remove .mfa file,
filter_more: $(SCR)/filter_more.sh Raw_data/cds
	bash $< $(word 2,$^)

#neutral typeN/(typeN+typeS)
neutral_ZnZs:$(SCR)/neutral_ZnZs.R Raw_data/cds
	$(RSCRIPT) $< $(word 2,$^) Results/ZnZs/neutral_ZnZs.txt

#GC% of data
GC_count:$(patsubst %, Results/GC/%.tsv, $(namelist))
Results/GC/%.tsv:$(SCR)/GC_perc.R Raw_data/cds/%
	$(RSCRIPT) $< $(word 2,$^) $@

#ensembl the mle est of chapter 3 as the default for coati-max/sampling
mle_90:$(SCR)/mle_90_ensemble.R
	$(RSCRIPT) $<

#plot the phase-proportions of 90 species
plot_phase_prop:$(SCR)/plot_phase_prop.R Results/phases
	$(RSCRIPT) $< $(word 2,$^) Figure/ZnZs/phase_prop.pdf

#plot the distribution of phases across 90 species
plot_phase_dist:$(SCR)/plot_phase_dist.R Results/phases
	$(RSCRIPT) $< $(word 2,$^) Figure/ZnZs/phase_dist.pdf

#plot dnds & typeN/typeNS dist
plot_dnds_ZnZs:$(SCR)/plot_dnds_ZnZs_dist.R Results/ZD_sum
	$(RSCRIPT) $< $(word 2,$^) Figure/ZnZs/dnds_ZnZs.pdf

#plot dnds & typeN/typeNS dist
plot_ZnZs_prop:$(SCR)/plot_ZnZs_prop.R Results/ZD_sum
	$(RSCRIPT) $< $(word 2,$^) Figure/ZnZs/ZnZs_prop.pdf

#plot omega vs dnds
plot_omega_dnds:$(SCR)/omega_vs_dnds_90.R ../chapter3/90/Results/PISE Results/dnds
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) Figure/ZnZs/omega_dnds.pdf





##################################################################################Mafft+sw
#Translate dna to aa
aa:$(patsubst %, Raw_data/aa/%, $(namelist))
Raw_data/aa/%:$(SCR)/dna2aa_90.R Raw_data/cds/%
	$(RSCRIPT) $< $(word 2,$^) $@

#Mafft alignment
align_mafft:$(patsubst %, Raw_data/align_mafft/%, $(namelist))
Raw_data/align_mafft/%:$(SCR)/align_90.sh Raw_data/aa/%
	bash $< $(word 2,$^) $@

#Map alignment back to dna
mapBack:$(patsubst %, Raw_data/mapped_cds/%, $(namelist))
Raw_data/mapped_cds/%:$(SCR)/map_cds2aa_90.R  Raw_data/align_mafft/% Raw_data/cds/%
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $@

filter_sw:$(SCR)/filter_align.sh Raw_data/QC2 Raw_data/mapped_cds
	bash $< $(word 2,$^) $(word 3,$^)


#keep independent gaps only
update_gap:$(patsubst %, Data/up_cds/%, $(namelist))
Data/up_cds/%:$(SCR)/update_gap_90.R Raw_data/mapped_cds/%
	$(RSCRIPT) $< $(word 2,$^) $@ $(WINDOW) $(WALL)

#Slide window method
sw_gap:$(patsubst %, Data/sw_cds/%, $(namelist))
Data/sw_cds/%:$(SCR)/sw_gap_90.R Data/up_cds/%
	$(RSCRIPT) $< $(word 2,$^) $@ $(WINDOW) $(WALL)

#phase analysis of sw
phase_sw_90:$(SCR)/phase_sw_90.R Data/sw_cds
	$(RSCRIPT) $< $(word 2,$^) Results/phases/phase_sw.txt

#Generate a table of sum stas: tyN,tyS,Nd,N,Sd,S,lenA,lenB
ZnZs_dnds_sw:$(patsubst %, Results/ZD_sum/%.sw.tsv, $(namelist))
Results/ZD_sum/%.sw.tsv:$(SCR)/ZnZs_dnds_sw.R Data/sw_cds/%
	$(RSCRIPT) $< $(word 2,$^) $@


#################################################################################coati-align-max
#max alignment
coati_max:$(patsubst %, Raw_data/align_max/%, $(namelist))
Raw_data/align_max/%:$(SCR)/coati_best.sh Raw_data/cds MLE/%.txt
	bash $< $(word 2,$^) $@ $(word 3,$^)

#generate table of:species_pair,locus,lenA,lenB,match_count,mismatch_count,gapA_count,gapB_count,gapA_len,gapB_len
jcdis_max:$(patsubst %, Raw_data/JCdis_sum/%.max.tsv, $(namelist))
Raw_data/JCdis_sum/%.max.tsv:$(SCR)/jcdis_sum.R Raw_data/align_max/%
	$(RSCRIPT) $< $(word 2,$^) $@ $(pat2)

#indel hex plot (fadj=0.5)
id_hex_max:$(SCR)/fadj_plot_90.R Raw_data/JCdis_sum
	$(RSCRIPT) $< $(word 2,$^) Raw_data/QC2/ Figure/QC2/id.hex.max.pdf $(pat2) 0.5

filter_coati_max:$(SCR)/filter_align.sh Raw_data/QC2 Raw_data/align_max
	bash $< $(word 2,$^) $(word 3,$^)

#phase analysis of coati-max
phase_coati_max_90:$(SCR)/phase_coati_max_90.R Raw_data/align_max
	$(RSCRIPT) $< $(word 2,$^) Results/phases/phase_coati_max.txt


#Generate a table of sum stas: Zn,Zs,Nd,N,Sd,S,lenA,lenB
ZnZs_dnds_max:$(patsubst %, Results/ZD_sum/%.max.tsv, $(namelist))
Results/ZD_sum/%.max.tsv:$(SCR)/ZnZs_dnds_max.R Raw_data/align_max/%
	$(RSCRIPT) $< $(word 2,$^) $@

#rolling ZnZs 
roll_ZnZs_max:$(SCR)/roll_ZnZs_nls.R Results/ZD_sum Results/GC Results/ZnZs/neutral_ZnZs.txt MLE90.tab.csv
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(word 5,$^) Figure/ZnZs/roll_ZnZs_max.pdf Results/ZnZs/coef_max.tsv $(pat2)

roll_ZnZs_max2:$(SCR)/roll_ZnZs_nlsLM.R Results/ZD_sum Results/GC Results/ZnZs/neutral_ZnZs.txt MLE90.tab.csv
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(word 5,$^) Figure/ZnZs/roll_ZnZs_max2.pdf Results/ZnZs/coef_max.tsv $(pat2)

#plot all-nls related coefficient
#plot_coef_max:$(SCR)/plot_nls_coef.R Results/ZnZs/coef_max.tsv
#	$(RSCRIPT) $< $(word 2,$^) Figure/ZnZs/coef_max.pdf



################################################################################coati-align-sampling
#sample alignment
coati_sampling:$(patsubst %, Raw_data/align_sampling/%, $(namelist))
Raw_data/align_sampling/%:$(SCR)/coati_sampling.sh Raw_data/cds MLE/%.txt
	bash $< $(word 2,$^) $@ $(word 3,$^)

filter_coati_sam:$(SCR)/filter_align.sh Raw_data/QC2 Raw_data/align_sampling
	bash $< $(word 2,$^) $(word 3,$^)


#phase analysis of coati-sample
phase_coati_sample_90:$(patsubst %, Results/phase_sample/%.txt, $(namelist))
Results/phase_sample/%.txt:$(SCR)/phase_coati_sample_90.R Raw_data/align_sampling/%
	$(RSCRIPT) $< $(word 2,$^) $@

#ensemble all phases into one file
phase_sample_ensembl:$(SCR)/phase_sample_ensembl.R
	$(RSCRIPT) $<

#Generate a table of sum stas: tyN,tyS,Nd,N,Sd,S,lenA,lenB >>>>
ZnZs_dnds_sample:$(patsubst %, Results/ZD_sum/%.sample.tsv, $(namelist))
Results/ZD_sum/%.sample.tsv:$(SCR)/ZnZs_dnds_sample.R Raw_data/align_sampling/%
	$(RSCRIPT) $< $(word 2,$^) $@

#rolling ZnZs
roll_ZnZs_sam:$(SCR)/roll_ZnZs_nlsLM.R Results/ZD_sum Results/GC Results/ZnZs/neutral_ZnZs.txt MLE90.tab.csv
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(word 5,$^) Figure/ZnZs/roll_ZnZs_sam2.pdf Results/ZnZs/coef_sam.tsv $(pat3)


###################################################################
supp:$(SCR)/supp_cha4_tab.R
	$(RSCRIPT) $<




######################
clean:
	@$(RM) Raw_data/raw_cds/*
	@$(RM) Raw_data/cds/*
	@$(RM) Raw_data/aa/*
	@$(RM) Raw_data/mafft_align/*
	@$(RM) Raw_data/mapped_cds/*
	@$(RM) Raw_data/align_max/*
	@$(RM) Raw_data/align_sampling/*
	@$(RM) Raw_data/QC/*
	@$(RM) Raw_data/QC2/*
