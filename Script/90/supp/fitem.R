library(tidyverse)
library(mixtools)

setwd("~/Dropbox (ASU)/Indel_project/Script/90/supp")
set.seed(403497)
num_gamma_components = 2

files <- fs::dir_ls(path="../../../test_90_species/Raw_data/JCdis_sum", glob="*.max.tsv")

output <- tibble(
    name=character(),
    f1=numeric(), shape1=numeric(), scale1=numeric(),
    f2=numeric(), shape2=numeric(), scale2=numeric()
)

for(f in files) {
    dat      <-  read_tsv(f,show_col_types=FALSE)
    dat_name <-  basename(str_replace(f, ".max.tsv$", ""))
    
    cat(str_glue("Fitting {dat_name}...\n\n"))

    dat <- dat %>% mutate(
        total_len = lenA+lenB,
        diff = abs(lenA-lenB),
        f = (gapA_len+gapB_len)/total_len,
        fdiff = diff/total_len,
        fadj = f-fdiff
        )

    # Calculate Distances
    d <- dat %>% filter(fadj <= 0.5 & fdiff <= 0.1) %>%
        transmute(
            p = mismatch_count/(match_count+mismatch_count),
            t = -0.75*log(1-4/3*p)
        )

    #fit to gamma mixture
    #add tiny offset to manage t=0 data points
    #quiet mixtools's `cat()` calls.
    quiet_em <- quietly(gammamixEM)
    model    <- quiet_em(d$t+1e-8, k=num_gamma_components)       
    model <- model$result

    shape <- model$gamma.pars[1,]
    scale <- model$gamma.pars[2,]
    f <- model$lambda
    output <- bind_rows(output, tibble(
        name = dat_name,
        f1 = f[1], shape1 = shape[1], scale1 = scale[1],
        f2 = f[2], shape2 = shape[2], scale2 = scale[2]
    ))
    write_csv(output, "fitem.csv")
}


