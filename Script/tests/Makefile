default:all
all:
.PHONY: test_update_gap test_sw_gap test_phase_indel_prob_wz coverage_report


RSCRIPT = Rscript --vanilla

###########################################unit test

test_update_gap: testthat/test_update_gap.R
	$(RSCRIPT) $<

test_sw_gap: testthat/test_sw_gap.R
	$(RSCRIPT) $<



##############################
test_phase_indel_prob_wz: test_phase_indel_prob_wz.R
	$(RSCRIPT) $<
	

##############################testing coverage report
coverage_report: coverage_report.R
	$(RSCRIPT) $<
