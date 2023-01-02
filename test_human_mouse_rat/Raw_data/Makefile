#Usage:make target

#GLobal variables
ARRAY := c"('human', 'mouse', 'rat')"
INPUT := ../Input_species/homo3.1.txt
HOMO  := geneId_get3.R


#Update
touch:
	touch cds_seq/*
	touch aa_seq/*
	touch Alignments/mafft_out/*
	touch Alignments/prank_out/*
	touch mapped_cds_mafft/*
	touch mapped_cds_prank/*
	touch pair_mafft/*
	touch pair_prank/*
.PHONY: touch

#Include more makefiles
include makeRaw.mk makeFig.mk
