#Generate Results/ Figure/

default: all
all: 	 Fig
Fig: 	 mafft_edge mafft_edge_report dnds znzs

.PHONY:  default all Fig
.PHONY:  mafft_edge mafft_edge_report dnds znzs

################################################PART I
RSCRIPT = Rscript --vanilla
sub     = Figure/mafft_edge


scr1    = ../../Script/edge_bias.R
scr2    = ../../Script/edge_report.R
scr3    = ../../Script/omega_ladder.R

#####################################################MAFFT


#Generate edge_bias table
mafft_edge:$(scr1) Mafft/updated_cds Mafft/mapped_cds Figure/mafft_edge/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(window)

#Generate the edge bias HTML report
mafft_edge_report:$(scr2) $(sub)/align.ori.l.txt $(sub)/align.better.l.txt $(sub)/align.ori.r.txt $(sub)/align.better.r.txt
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(word 5,$^) Figure/mafft.edge.report.html



######################################################dnds vs znzs
#Count the dNdS value across the genome 
dnds:../../Script/dnds_pair.R Mafft/mapped_cds
	$(RSCRIPT) $< $(word 2,$^) Results/dnds.txt

znzs:../../Script/znzs_pair.R Mafft/mapped_cds
	$(RSCRIPT) $< $(word 2,$^) Results/znzs.txt 3

