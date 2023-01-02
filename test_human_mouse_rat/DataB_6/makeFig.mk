default: all
all: 	 Fig
Fig: 	 mafft_phase prank_phase mafft_dis prank_dis dNdS zNzS_rev  

.PHONY:  default all Fig  
.PHONY:  mafft_phase prank_phase mafft_dis prank_dis mafft_phase_eff dNdS zNzS_rev omega_ladder w_zed_corr sim_zed



RSCRIPT  = Rscript --vanilla 

scr1     = ../../Script/?.R
scr1.1   = ../../Script/?.R
scr2     = ../../Script/?.R

scr3     = ../../Script/dnds_genome.R
scr4     = ../../Script/znzs_genome_rev.R

scr5     = ../../Script/omega_ladder.R
scr6     = ../../Script/zed_omega_corr.R
scr7     = ../../Script/sim_zed.R

########################################################PART II
file.1 = Results/phase.mafft.txt
file.2 = Results/phase.prank.txt
file.3 = Results/dis.mafft.txt
file.4 = Results/dis.prank.txt

fig.1  = Figure/phase.mafft.pdf
fig.2  = Figure/phase.prank.pdf
fig.3  = Figure/dis.mafft.pdf
fig.4  = Figure/dis.prank.pdf

fig.5  = Figure/phase.mafft.eff.pdf
fig.6  = Figure/phase.prank.eff.pdf


#Count the proportion of phase 0,1,2 Indels
mafft_phase:$(scr1) Mafft/mapped_cds 
	$(RSCRIPT) $< $(word 2,$^) $(file.1) $(fig.1) $(window)

prank_phase:$(scr1) Prank/mapped_cds 
	$(RSCRIPT) $< $(word 2,$^) $(file.2) $(fig.2) $(window)

#Generate displacement table 
mafft_dis:$(scr2) Mafft/updated_cds Mafft/mapped_cds 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(file.3) $(fig.3) $(window) 

prank_dis:$(scr2) Prank/updated_cds Prank/mapped_cds 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(file.4) $(fig.4) $(window)

#Count the effective phase distribution of phase 1,2 Indels
mafft_phase_eff:$(scr1.1) Mafft/mapped_cds
	$(RSCRIPT) $< $(word 2,$^) $(fig.5)

prank_phase_eff:$(scr1.1) Prank/mapped_cds
	$(RSCRIPT) $< $(word 2,$^) $(fig.6)




#############################################################PART III
file.5 = Results/dNdS.txt
file.6 = Results/zNzS_rev.txt

#Count the dNdS value across the genome (1101 files)
dNdS:$(scr3) Mafft/mapped_anc 
	$(RSCRIPT) $< $(word 2,$^) $(file.5)


#Count the zNzS value across the genome (1101 files)
zNzS_rev:$(scr4) Mafft/mapped_anc 
	$(RSCRIPT) $< $(word 2,$^) $(file.6)


############################################################PART IV

file.7_1  = Results/omega_ladder.txt
file.7_2  = Results/omega_bin.txt
file.8    = Results/omega_zed.txt
fig.7     = Figure/omega_dis.pdf
fig.8     = Figure/omega_zed.pdf

#Count the dNdS value across every single gene and create a table for geneIds:dnds:#ofgaps
omega_ladder:$(scr5) Mafft/mixed_cds
	$(RSCRIPT) $< $(word 2,$^) $(file.7_1) $(file.7_2) $(fig.7)

#Build correlation between omega(w) and zed(znzs)
w_zed_corr:$(scr6) Mafft/mixed_cds/ $(file.7_1) $(file.7_2)
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(file.8) $(fig.8)


#Simulate indels to obtain the expected n/s across phase 0,1,2
sim_zed: $(scr7) Simulation/codon_freq.txt
	$(RSCRIPT) $< $(word 2,$^) 





