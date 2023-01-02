#!/bin/bash

# Vars
RSCRIPT = Rscript --vanilla 
RM      = rm -fi
Dir     = ../Raw_data_2
scr     = ../../Script

default: all
all:     mafft_upc mafft_map mafft_anc prank_upc prank_map prank_anc mafft_mix 

.PHONY:  all default clean 
.PHONY:  mafft_upc mafft_map mafft_anc prank_upc prank_map prank_anc mafft_mix


#######################################################

NAMELIST := $(shell cat $(Dir)/nameList.txt )


mafft_upc:	 $(patsubst %, Mafft/updated_cds/%.fa,    $(basename $(NAMELIST)))
mafft_map:	 $(patsubst %, Mafft/mapped_cds/%.fa,     $(basename $(NAMELIST)))
prank_upc:	 $(patsubst %, Prank/updated_cds/%.fa,    $(basename $(NAMELIST)))
prank_map:	 $(patsubst %, Prank/mapped_cds/%.fa,     $(basename $(NAMELIST)))

mafft_mix:	 $(patsubst %, Mafft/mixed_cds/%.fa,      $(basename $(NAMELIST)))


Mafft/updated_cds/%.fa:  | $(Dir)/mapped_cds_mafft/%.fa
Mafft/mapped_cds/%.fa:   | Mafft/updated_cds/%.fa 
Mafft/mixed_cds/%.fa:    | Mafft/updated_cds/%.fa 


#######################################################


#Update all gaps
Mafft/updated_cds/%.fa: $(scr)/update_gap_2.R $(Dir)/mapped_cds_mafft/%.fa Mafft/updated_cds/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window) $(wall)
.SECONDARY:Mafft/updated_cds/%.fa

Prank/updated_cds/%.fa: $(scr)/update_gap_2.R $(Dir)/mapped_cds_prank/%.fa Prank/updated_cds/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window) $(wall)
.SECONDARY:Prank/updated_cds/%.fa


#Find the best phase
Mafft/mapped_cds/%.fa: $(scr)/sw_gap_2.R   Mafft/updated_cds/%.fa Mafft/mapped_cds/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window) $(wall)
.SECONDARY:Mafft/mapped_cds/%.fa

Prank/mapped_cds/%.fa: $(scr)/sw_gap_2.R   Prank/updated_cds/%.fa Prank/mapped_cds/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window) $(wall)
.SECONDARY:Prank/mapped_cds/%.fa


#Find the common ancestor and update the focal seqs
mafft_anc: $(scr)/add_com_anc.R  Mafft/mapped_cds Mafft/mapped_anc/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) 

prank_anc: $(scr)/add_com_anc.R  Prank/mapped_cds Prank/mapped_anc/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^)



###########################################################PART III

#Find the best phase and update inside the folder
Mafft/mixed_cds/%.fa: $(scr)/sw_gap_2.1.R   Mafft/updated_cds/%.fa Mafft/mixed_cds/ 
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window) $(wall)
.SECONDARY:Mafft/mixed_cds/%.fa





clean:
	
	@$(RM) Mafft/updated_cds/*
	@$(RM) Mafft/mapped_cds/*
	@$(RM) Mafft/mapped_anc/*
	@$(RM) Prank/updated_cds/*
	@$(RM) Prank/mapped_cds/*
	@$(RM) Prank/mapped_anc/*
