#!/bin/bash
#need to generate 100 folders of GS, GS_trim, JsonD beforehand.

default:all
all:

.PHONY: para_space gill100 gsam100 rm_gap gimp100 plot_sim plot_im plot_sim_err plot_im_err
.PHONY: get_max_Id gimp100_90 comp_r_prop_90 supp

#command
RSCRIPT= Rscript --vanilla
RM     = rm -i

#path
SCR  = ../Script/chapter3

len        := 999
omega_z    := 1
sim_size   := 500
sample_size:= 1000
ncore      := 4
target     := NULL

#################################################100 simulations
#Pick a parameter space (100)
para_space:$(SCR)/para_space.R
	$(RSCRIPT) $<

#simulations
nameL  := $(shell test -f index.txt && cat index.txt)
gill100:  $(patsubst %, Gs/%, $(nameL))
gsam100:  $(patsubst %, Results/Gse/%.est.json, $(nameL))
gimp100:  $(patsubst %, Results/Pise/%.est.json,$(nameL))
Gs:          gill100
Results/Gse: gsam100
Results/Pise:gimp100

#100 gillespie simulations
Gs/%: $(SCR)/haha_gillespie100.R trueP.100.txt
	$(RSCRIPT) $< $@ $(word 2,$^) $(len) $(omega_z) $(sim_size)

#100 gillespie simu estimates
Results/Gse/%.est.json: $(SCR)/gillespie_sim_est2.R Gs/%
	$(RSCRIPT) $< $(word 2,$^) $@

#generate original seqs
rm_gap:runbash100.sh
	bash $< 100 500

#100/10k switch
Results/Pise/%.est.json: $(SCR)/phase_im_sampling_1e4.R Gs_trim/% JsonD/%/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $@ $(sample_size) 64


#plot the true vs sim estimates
plot_sim: $(SCR)/plot_trueVsSim.R
	$(RSCRIPT) $< $(target) 1

plot_sim_err: $(SCR)/plot_err.R
	$(RSCRIPT) $< $(target) 1

plot_im: $(SCR)/plot_trueVsSim.R
	$(RSCRIPT) $< $(target) 0

plot_im_err: $(SCR)/plot_err.R
	$(RSCRIPT) $< $(target) 0

##############################################90 species
#hmr data
#sw data_3,_6,_9,_12
# get_sw_Id:$(SCR)/get_sw_Id.sh
# 	bash $<
#
# name_hmr:= $(shell test -f hmr/index.txt && cat hmr/index.txt)
# gimp100_hmr:$(patsubst %, hmr/Results/%.est.json, $(name_hmr))
# hmr/Results:gimp100_hmr
# hmr/Results/%.est.json: $(SCR)/phase_im_sampling_1e4.R hmr/Data/% hmr/JsonD/%/
# 	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $@ $(sample_size)

#grab the raw cds filterd from align_max alignments
get_max_Id:$(SCR)/get_max_Id.sh
	bash $<

#90 species MLE
name90:= $(shell cat ../test_90_species/Species.txt)
gimp100_90:$(patsubst %, 90/Results/PISE/%.est.json, $(name90))
90/Results/PISE/%.est.json:$(SCR)/phase_im_sampling_1e4.R ../test_90_species/Raw_data/cds/% 90/JsonD/%/ 90/species/%.txt
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $@ $(sample_size) 61


#compare indel-rate vs indel-proportions
comp_r_prop_90:$(SCR)/indel_rate_vs_phase_90.R 90/Results/PISE ../test_90_species/Results/phases
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) 90/Figure/r_vs_prop.pdf


####################################################supplementry table
supp:$(SCR)/supp_cha3_tab.R
	$(RSCRIPT) $<


################################
clean:
