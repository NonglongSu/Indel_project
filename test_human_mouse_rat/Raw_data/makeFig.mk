#make similaity figures

default: all
all: 	 Fig
Fig:     mafft_sim prank_sim

.PHONY: default all Fig mafft_sim prank_sim

RSCRIPT = Rscript --vanilla
SCR     = ../../Script/QC

#############################
fig.1 = Figure/sim.mafft.pdf
fig.2 = Figure/sim.prank.pdf

#Measure the similarity in species pairs
mafft_sim: $(SCR)/sim_dist_plot.R Alignments/mafft_out mapped_cds_mafft $(INPUT)
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) Figure/sim.mafft.pdf

prank_sim: $(SCR)/sim_dist_plot.R Alignments/prank_out mapped_cds_prank $(INPUT)
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(fig.2)

