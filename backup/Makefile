#!/bin/bash

default:all
all:

.PHONY: sim_data gill_sampler impo_sampler impo_sampler_pal
.PHONY: para_space gill100 gsam100 gimp100
.PHONY: mafft mapCds ups sws  get_sw_Id

#command
RSCRIPT=Rscript --vanilla
RM=rm -i
MAFFT=mafft --maxiterate 1000 --globalpair

#path
SCR=../Script/chapter3

len        := 2001
omega_z    := 1
sim_size   := 100
sample_size:= 100
ncore      := 4
##############################################Model construction

sim_data: $(SCR)/haha_gillespie4.R
	$(RSCRIPT) $< $(word 2,$^) haha/ $(len)  $(omega_z) $(sim_size)
#alternative:haha_gillespie3.R

gill_sampler: $(SCR)/gillespie_sim_est2.R haha
	$(RSCRIPT) $< $(word 2,$^) gse.json
#alternative:gillespie_sim_est.R

impo_sampler: $(SCR)/phase_im_sampling_multi.R fasta Json/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) par.est.json weight_mat.txt

impo_sampler_pal: $(SCR)/phase_im_sampling_parallel.R fasta Json/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) par.est.json weight_mat.txt $(ncore)

#################################################100 simulations
#Pick a parameter space (100)
para_space:$(SCR)/para_space.R
	$(RSCRIPT) $<

nameL  := $(shell test -f index.txt && cat index.txt)
gill100:  $(patsubst %, Gs/%, $(nameL))
gsam100:  $(patsubst %, Results/Gse/%.est.json, $(nameL))
gimp100:  $(patsubst %, Results/Pise/%.est.json, $(nameL))
Gs:          gill100
Results/Gse: gsam100
Results/Pise:gimp100


#100 gillespie simulations
Gs/%: $(SCR)/haha_gillespie100.R trueP.100.txt
	$(RSCRIPT) $< $@ $(word 2,$^) $(len) $(omega_z) $(sim_size)

#100 gillespie simu estimates
Results/Gse/%.est.json: $(SCR)/gillespie_sim_est2.R Gs/%
	$(RSCRIPT) $< $(word 2,$^) $@

#plot the true vs sim estimates
plot_trueVsSim: $(SCR)/plot_trueVsSim.R
	$(RSCRIPT) $<

#100 importance sampling
#Results/Pise/%.est.json: $(SCR)/phase_im_sampling_100.R Gs_trim/% JsonD/%/
#	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $@

#10k
Results/Pise/%.est.json: $(SCR)/phase_im_sampling_1e4.R Gs_trim/% JsonD/%/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $@ $(sample_size)


##############################################part II
DATABASE:=hmr/Raw/nameList.txt
NAMELIST:= $(shell test -f $(DATABASE) && cat $(DATABASE))
window  := 6
wall    := 12

hmr/Raw/align:           mafft
hmr/Raw/align/%.fa:      |hmr/Raw/aa/%.fa
mafft:                   $(patsubst %, hmr/Raw/align/%.fa, $(NAMELIST))

hmr/Raw/mapped_cds:	     mapCds
hmr/Raw/mapped_cds/%.fa: |hmr/Raw/align/%.fa
mapCds:			             $(patsubst %, hmr/Raw/mapped_cds/%.fa, $(NAMELIST))

hmr/Data_6/ups:		     ups
hmr/Data_6/ups/%.fa:	 |hmr/Raw/mapped_cds/%.fa
ups:                   $(patsubst %, hmr/Data_6/ups/%.fa, $(NAMELIST))


hmr/Data_6/sws:		     sws
hmr/Data_6/sws/%.fa:	 |hmr/Data_6/ups/%.fa
sws:                   $(patsubst %, hmr/Data_6/sws/%.fa, $(NAMELIST))



#pairwise alignment
hmr/Raw/align/%.fa: hmr/Raw/aa/%.fa
	$(MAFFT) $< > $@
.SECONDARY:hmr/Raw/align/%.fa

#map aa->dna
hmr/Raw/mapped_cds/%.fa: ../Script/map_cds2aa.R  hmr/Raw/align/%.fa hmr/Raw/cds/%.fa
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^)  hmr/Raw/mapped_cds/
.SECONDARY: hmr/Raw/mapped_cds/%.fa


#update gaps
hmr/Data_6/ups/%.fa: $(SCR)/update_gap.R hmr/Raw/mapped_cds/%.fa hmr/Data_6/ups/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window) $(wall)
.SECONDARY:hmr/Data_6/ups/%.fa


#sliding window method
hmr/Data_6/sws/%.fa: $(SCR)/sw_gap.R hmr/Data_6/ups/%.fa hmr/Data_6/sws/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window) $(wall)
.SECONDARY:hmr/Data_6/sws/%.fa


#############################################
get_sw_Id:$(SCR)/get_sw_Id.sh hmr/Data_6/sws  hmr/Raw/cds/
	bash $< $(word 2,$^) $(word 3,$^) hmr/coatiM/swId.txt hmr/coatiM/Raw/
