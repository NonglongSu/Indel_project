#!/bin/bash
#90 species analysis

default: all
all:

.PHONY:  default all clean
.PHONY:  para100  rms100 em100_5e em100_6e em100_7e plotEM100 plotRmse em90 boot90 plotEM plotCI
.PHONY:  cutFilter nmkb100 nmkb100.1 supp

RSCRIPT = Rscript --vanilla
SCR     = ../Script/90
RM      = rm -i

###################################PART I:EM simulation
#EM 100 simulations
#setup parameter space(true value)
para100:../Script/haha_para_space.R
	$(RSCRIPT) $<


namelist:=$(shell test -f Results/simuList.txt && cat Results/simuList.txt)

#generate a relationship between rms and 100 sample size{5e..6e..7e}.
Results/rms:rms100
rms100:     $(patsubst %, Results/rms/%.json, $(namelist))

#generate 100 simulations for sample size of {5e,6e,7e}
Results/est100:em100_5e em100_6e em100_7e
em100_5e:      $(patsubst %, Results/est100/%.5e.json, $(namelist))
em100_6e:      $(patsubst %, Results/est100/%.6e.json, $(namelist))
em100_7e:      $(patsubst %, Results/est100/%.7e.json, $(namelist))

#fail!
Results/rms/%.json:../Script/haha_61_aitken.R Results/truePar_100.txt
	$(RSCRIPT) $< $@ $(word 2,$^)

#5e
Results/est100/%.5e.json:$(SCR)/dumbEm_simu_100.R Results/truePar_100.txt
	$(RSCRIPT) $< $@ $(word 2,$^)

#6e
Results/est100/%.6e.json:$(SCR)/dumbEm_simu_100.R Results/truePar_100.txt
	$(RSCRIPT) $< $@ $(word 2,$^)

#7e
Results/est100/%.7e.json:$(SCR)/dumbEm_simu_100.R Results/truePar_100.txt
	$(RSCRIPT) $< $@ $(word 2,$^)

#plot true~est
plotEM100:$(SCR)/em_100_est.R Results/truePar_100.txt Results/est100
	for num in 5e 6e 7e; do \
		$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) Figure/EM/est100.$${num}.pdf Figure/EM/qqplot.$${num}.pdf Results/rmse.$${num}.txt $${num}; \
	done

#plot sample size~rms
plotRmse:$(SCR)/plot_rmse.R Results
	$(RSCRIPT) $< $(word 2,$^) Figure/EM/rmse_size.pdf


##############################################PartII 90 species
namelist1:=$(shell test -f Results/specList.txt && cat Results/specList.txt)
Results/est90:	em90
Results/boot90: boot90
em90:	$(patsubst %, Results/est90/%.json,    $(namelist1))
boot90: $(patsubst %, Results/boot90/%.10.json, $(namelist1))

#EM 90 species
Results/est90/%.json:$(SCR)/dumbEm_90.R Raw_data/concat
	$(RSCRIPT) $< $@ $(word 2,$^)

#bootstrapping 90 ests
Results/boot90/%.json:$(SCR)/dumbEm_boot_90.R Raw_data/concat Results/est90
	$(RSCRIPT) $< $@ $(word 2,$^) $(word 3,$^)


#Plot 90 EM ests.
plotEM: $(SCR)/em90_est.R Results/est90
	$(RSCRIPT) $< $(word 2,$^) Figure/EM/est90.pdf Figure/EM/omega_tau.pdf


plotCI: $(SCR)/em90_CI.R Results/boot90 Results/specList.txt Results/CI90
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^)  Figure/EM/omega_tau_CI.pdf

#Artifact.
#Measuring the filtering method.
cutFilter:$(SCR)/cut_filter_90.R Raw_data/concat Raw_data/raw_cds Results/est90
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) Figure/EM/loss_of_filter.pdf


########################################PARI III performance study
#Optimization -- "nmkb"(Nelder-Mead  algorithm  for  derivative-free  optimization)

#generate 100 simulations for sample size of {5e,6e,7e}
Results/nmkb100:nmkb100
nmkb100:$(patsubst %, Results/nmkb100/%.5e.json, $(namelist))

Results/nmkb100/%.5e.json:$(SCR)/nmkb_100.R Results/truePar_100.txt
	$(RSCRIPT) $< $@ $(word 2,$^)

#modified nmkb package
Results/nmkb100.1:nmkb100.1
nmkb100.1:$(patsubst %, Results/nmkb100.1/%.5e.json, $(namelist))

Results/nmkb100.1/%.5e.json:$(SCR)/nmkb_100_1.R Results/truePar_100.txt
	$(RSCRIPT) $< $@ $(word 2,$^)


##########################################PART IV add supp tables
supp:$(SCR)/supp_cha1_tab.R 
	$(RSCRIPT) $< 






