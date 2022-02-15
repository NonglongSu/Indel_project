#!/bin/bash

default:all
all:

.PHONY: para_space gill100 gsam100 gimp100

################################################
RSCRIPT=Rscript --vanilla
RM=rm -i

#path
SCR=../Script/chapter3

len        := 999
sim_size   := 500
sample_size:= 100
omega_z    := 0
omega_l    := 12

#################################################100 simulations
#Pick a parameter space (100)
para_space:$(SCR)/para_space_wz.R
	$(RSCRIPT) $<

nameL  := $(shell test -f index.txt && cat index.txt)
gill100:  $(patsubst %, Gs/%, $(nameL))
gsam100:  $(patsubst %, Results/Gse/%.est.json, $(nameL))
gimp100:  $(patsubst %, Results/Pise/%.est.json,$(nameL))
Gs:          gill100
Results/Gse: gsam100
Results/Pise:gimp100


#100 gillespie simulations
Gs/%: $(SCR)/haha_gillespie100_wz.R trueP.100.txt
	$(RSCRIPT) $< $@ $(word 2,$^) $(len) $(omega_z) $(sim_size) $(omega_l)

#100 gillespie simu estimates
Results/Gse/%.est.json: $(SCR)/gillespie_sim_est_wz.R Gs/%
	$(RSCRIPT) $< $(word 2,$^) $@ $(omega_l)



#100 importance sampling
#Results/Pise/%.est.json: $(SCR)/phase_im_sampling_100.R Gs_trim/% JsonD/%/
#	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $@

#10k
Results/Pise/%.est.json: $(SCR)/phase_im_sampling_1e4.R Gs_trim/% JsonD/%/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $@ $(sample_size)



#################################
#plot the true vs sim estimates | sim vs im
