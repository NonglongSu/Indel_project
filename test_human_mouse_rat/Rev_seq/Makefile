#!/bin/bash

default: all
all:	 rev_cds aa mafft mapped_mafft paired_mafft
.PHONY:  default all clean
.PHONY:  rev_cds aa mafft mapped_mafft paired_mafft 
##################################################

#Set var
RSCRIPT  = Rscript --vanilla
MAFFT    = mafft --maxiterate 1000 --globalpair
RM       = rm -i -f
SCR      = ../../Script/rev_seq

ARRAY    := c"('human', 'mouse', 'rat')"
NAMELIST := $(shell test -f nameList.txt && cat nameList.txt)

########################################################PART I
# Raw_data:	rev_cds
# Raw_data:	aa
# Raw_data:	mafft
# Raw_data:	mapped_mafft
# Raw_data:	paired_mafft

rev_cds:		$(patsubst %, Raw_data/%.rcds.fa,  $(NAMELIST))
aa:				$(patsubst %, Raw_data/%.aa.fa,    $(NAMELIST))
mafft:			$(patsubst %, Raw_data/%.mafft.fa, $(NAMELIST))
mapped_mafft:	$(patsubst %, Raw_data/%.mapped.fa,$(NAMELIST))
paired_mafft:	$(patsubst %, Raw_data/%.paired.fa,$(NAMELIST))


Raw_data/%.rcds.fa:   | ../Raw_data/cds_seq/%.fa
Raw_data/%.aa.fa:     | Raw_data/%.rcds.fa
Raw_data/%.mafft.fa:  | Raw_data/%.aa.fa 
Raw_data/%.mapped.fa: | Raw_data/%.mafft.fa 
Raw_data/%.paired.fa: | Raw_data/%.mapped.fa 



#reverse seq
Raw_data/%.rcds.fa:$(SCR)/rev_cds.R ../Raw_data/cds_seq/%.fa
	$(RSCRIPT) $< $(word 2,$^) $@
.SECONDARY: Raw_data/%.rcds.fa

#translation
Raw_data/%.aa.fa:$(SCR)/dna2aa.R Raw_data/%.rcds.fa
	$(RSCRIPT) $< $(word 2,$^) $@ ${ARRAY}
.SECONDARY: Raw_data/%.aa.fa

#mafft
Raw_data/%.mafft.fa:Raw_data/%.aa.fa
	$(MAFFT) $< > $@
.SECONDARY: Raw_data/%.mafft.fa

#map aa to cds
Raw_data/%.mapped.fa:$(SCR)/map_cds2aa.R  Raw_data/%.mafft.fa Raw_data/%.rcds.fa
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) $@
.SECONDARY: Raw_data/%.mapped.fa

#remove reference
Raw_data/%.paired.fa:$(SCR)/pair_trim_each.R Raw_data/%.mapped.fa
	$(RSCRIPT) $< $(word 2,$^) $@
.SECONDARY: Raw_data/%.paired.fa



#Run bash checkEmpty.sh | checkName.sh folder to fix any broken pipe.
##############################################################
clean:
	@$(RM) Raw_data/*


help:
	@echo "help me, Doctor Ziqi!"
