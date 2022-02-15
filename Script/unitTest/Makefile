default:all
all:
.PHONY: test_update_gap test_sw_gap test_phase_indel_prob_wz


RSCRIPT = Rscript --vanilla

###########################

test_udpdate_gap: test_update_gap.R
	$(RSCRIPT) $<


test_sw_gap: test_sw_gap.R
	$(RSCRIPT) $<

##############################
test_phase_indel_prob_wz: test_phase_indel_prob_wz.R
	$(RSCRIPT) $<
