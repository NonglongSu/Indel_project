
default: all
all: $(SUBDIRS)

SUBDIRS := $(wildcard */.)
RSCRIPT  = Rscript --vanilla


.PHONY: all default $(SUBDIRS) 
.PHONY: mafft_phase mafft_phase_in_del mafft_pos_diff

$(SUBDIRS): 
	$(MAKE) -C $@



###################################Four window combined plot 
#Count the proportion of phase 0,1,2 indels         need to test!!!
mafft_phase: ../Script/indel_phase.R
	$(RSCRIPT) $< Results/phase.mafft.pdf

#Plot ins/del phase proportions
mafft_phase_in_del: ../Script/in_del_plot.R 
	$(RSCRIPT) $< Results/Phase.mafft.ins.del.pdf 
	
#Count the proportion of eff phase
mafft_phase_eff: ../Script/indel_phase_eff.R
	$(RSCRIPT) $< Results/phase.mafft.eff.pdf
	
#Positional difference distribution plot
mafft_pos_diff: ../Script/indel_dis_plot.R
	$(RSCRIPT) $< Results/dis.mafft.pdf














