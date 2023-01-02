#!/bin/bash

default:all
all:    mafft_upc mafft_sw

.PHONY: all default clean
.PHONY: mafft_upc mafft_sw prank_upc prank_sw


#Set up vars
RSCRIPT = Rscript --vanilla
RM      = rm -i -f
scr     = ../../Script
Dir     = ../Raw_data

NAMELIST := $(shell cat $(Dir)/nameList.txt)
#####################################################PART I
mafft_upc:  $(patsubst %, Mafft/updated_cds/%.fa,  $(basename $(NAMELIST)))
mafft_sw:   $(patsubst %, Mafft/mapped_cds/%.fa,   $(basename $(NAMELIST)))

Mafft/updated_cds/%.fa:  | $(Dir)/pair_mafft/%.fa
Mafft/mapped_cds/%.fa:   | Mafft/updated_cds/%.fa

prank_upc:  $(patsubst %, Prank/updated_cds/%.fa,   $(basename $(NAMELIST)))
prank_map:  $(patsubst %, Prank/mapped_cds/%.fa,    $(basename $(NAMELIST)))

Prank/updated_cds/%.fa:  | $(Dir)/pair_prank/%.fa
Prank/mapped_cds/%.fa:   | Prank/updated_cds/%.fa


#######################################################PART II
#Update gaps
Mafft/updated_cds/%.fa: $(scr)/update_gap.R $(Dir)/pair_mafft/%.fa Mafft/updated_cds/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window) $(wall)
.SECONDARY:Mafft/updated_cds/%.fa

#Sliding windows 
Mafft/mapped_cds/%.fa:  $(scr)/sw_gap.R Mafft/updated_cds/%.fa Mafft/mapped_cds/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window) $(wall)
.SECONDARY:Mafft/mapped_cds/%.fa


#>>>>>>>>>>>>>>>>>>>>>>>>>>>
#Prank/updated_cds/%.fa: $(scr)/update_gap.R $(Dir)/pair_prank/%.fa Prank/updated_cds/ 
#	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window) $(wall)
#.SECONDARY:Prank/updated_cds/%.fa
#
#Prank/mapped_cds/%.fa:  $(scr)/sw_gap.R Prank/updated_cds/%.fa Prank/mapped_cds_Out/
#	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $(window) $(wall)
#.SECONDARY:Prank/mapped_cds/%.fa



###################################
clean:
	@$(RM) Mafft/updated_cds/*
	@$(RM) Mafft/mapped_cds/*
	

