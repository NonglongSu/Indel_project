#!/bin/bash

#Set var
RSCRIPT = Rscript --vanilla
MAFFT   = mafft --maxiterate 1000 --globalpair
PRANK   = prank -F -f=fasta
RM      = rm -i -f
SCR     = ../../Script

DATABASE   := geneId.txt
NAMELIST   := $(shell test -f $(DATABASE) && cat  $(DATABASE) | cut -f 1 | sed -e 's/"//g')
NAMELIST1  := $(shell test -f nameList.txt && cat nameList.txt)

#########################
default: all
all:

.PHONY: default all clean
.PHONY: homoCall cds test_mul3 test_N test_stop test_preStop filter_cds
.PHONY: addStop aa guideTree mafft prank mapped_mafft mapped_prank pair_mafft pair_prank

###########################################Extract CDS
homoCall: $(SCR)/$(HOMO) $(INPUT)
	$(RSCRIPT) $< $(word 2,$^) geneId.txt

cds_seq:      cds
cds:          $(patsubst %, cds_seq/%.fa, $(NAMELIST))
cds_seq/%.fa: $(SCR)/cds_get.R $(DATABASE)
	$(RSCRIPT) $< $@ $(word 2,$^) cds_seq/
.SECONDARY: cds_seq/%.fa


###########################################Quality control
#Record files not-multiple-of-3.
test_mul3: $(SCR)/QC/test_mul_3.R cds_seq/
	$(RSCRIPT) $< $(word 2,$^) QC/multi_no_3.txt

#Record files that have 'N' nt.
test_N: $(SCR)/QC/test_NNN.R cds_seq/
	$(RSCRIPT) $< $(word 2,$^) QC/NNN.txt

#Record files without stop codons.
test_stop: $(SCR)/QC/test_std.R cds_seq/
	$(RSCRIPT) $< $(word 2,$^) QC/noStop.txt

#Record files with pre stop codons.
test_preStop: $(SCR)/QC/test_preStop.R cds_seq/
	$(RSCRIPT) $< $(word 2,$^) QC/preStop.txt

#combine previous files and filter it out
filter_cds: $(SCR)/QC/filter_cds.sh QC/multi_no_3.txt QC/NNN.txt QC/noStop.txt QC/preStop.txt
	bash $< $(word 2,$^) $(word 3,$^) $(word 4,$^) $(word 5,$^)

#create a NAMELIST1 target (nameList.txt) 
nameList.txt: cds_seq/
	ls $< | sed -e 's/.fa//g' > nameList.txt


########################################################PART II
aa_seq:			aa
Alignments/mafft_out:	mafft
Alignments/prank_out:	prank
mapped_cds_mafft:	mapped_mafft
mapped_cds_prank:	mapped_prank

aa:      	$(patsubst %, aa_seq/%.fa, $(NAMELIST1))
mafft:   	$(patsubst %, Alignments/mafft_out/%.fa, $(NAMELIST1))
prank:   	$(patsubst %, Alignments/prank_out/%.fa, $(NAMELIST1))
mapped_mafft:   $(patsubst %, mapped_cds_mafft/%.fa, $(NAMELIST1))
mapped_prank:   $(patsubst %, mapped_cds_prank/%.fa, $(NAMELIST1))

aa_seq/%.fa:                 | cds_seq/%.fa
Alignments/mafft_out/%.fa:   | aa_seq/%.fa
Alignments/prank_out/%.fa:   | aa_seq/%.fa
mapped_cds_mafft/%.fa:	     | Alignments/mafft_out/%.fa cds_seq/%.fa
mapped_cds_prank/%.fa:	     | Alignments/prank_out/%.fa cds_seq/%.fa


#Translation
aa_seq/%.fa:$(SCR)/dna2aa.R cds_seq/%.fa
	$(RSCRIPT) $< $(word 2,$^) aa_seq/ ${ARRAY}
.SECONDARY: aa_seq/%.fa

#MAFFT
Alignments/mafft_out/%.fa: aa_seq/%.fa
	$(MAFFT) $< > $@
.SECONDARY: Alignments/mafft_out/%.fa

#PRANK
guideTree: $(SCR)/make_guide_tree.R $(INPUT)
	$(RSCRIPT) $< $(word 2,$^) > guide.tree

Alignments/prank_out/%.fa: aa_seq/%.fa
	$(PRANK) -d=$< -t=guide.tree -o=Alignments/prank_out/$*
	test -f Alignments/prank_out/$*.best.fas && \
	samtools faidx Alignments/prank_out/$*.best.fas $(shell cut -f3 $(INPUT)) > $@ && \
	rm Alignments/prank_out/$*.best*
.SECONDARY: Alignments/prank_out/%.fa

#map aa to cds
mapped_cds_mafft/%.fa: $(SCR)/map_cds2aa.R  Alignments/mafft_out/%.fa cds_seq/%.fa
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) mapped_cds_mafft/
.SECONDARY: mapped_cds_mafft/%.fa

mapped_cds_prank/%.fa: $(SCR)/map_cds2aa.R  Alignments/prank_out/%.fa cds_seq/%.fa
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) mapped_cds_prank/
.SECONDARY: mapped_cds_prank/%.fa

#remove reference
pair_mafft:$(SCR)/pair_trim.R mapped_cds_mafft pair_mafft/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^)

pair_prank:$(SCR)/pair_trim.R mapped_cds_prank pair_prank/
	$(RSCRIPT) $< $(word 2,$^) $(word 3,$^) 

#Run bash checkEmpty.sh | checkName.sh folder to fix any broken pipe.
##############################################################
clean:
	@$(RM) cds_seq/*
	@$(RM) aa_seq/*
	@$(RM) Alignments/mafft_out/*
	@$(RM) Alignments/prank_out/*
	@$(RM) mapped_cds_mafft/*
	@$(RM) mapped_cds_prank/*
	@$(RM) pair_mafft/*
	@$(RM) pair_prank/*

help:
	@echo "help me, Doctor Ziqi!"
