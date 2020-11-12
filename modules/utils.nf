nextflow.enable.dsl=2

process make_phenotypes {
    input:
        path(phenotypes)
        path(to_include)
    output:
        path("phenotypes.txt"), emit: pheno
    script:
        """
        #!/usr/bin/env Rscript
        library(vroom)
        library(tidyverse)
        to_include <- read_delim($to_include),
                                 col_names = FALSE,
                                 delim = " ",
                                 trim_ws = TRUE) %>% 
            mutate(ids = paste0("f.", X1, ".0.0"))
        all_phenotypes <- vroom($phenotypes,
                                col_names = TRUE,
                                num_threads = cores,
                                col_select = c(f.eid,
                                               all_of(to_include\$ids))) %>% 
            rename(FID = f.eid) %>% 
            rename_at(vars(to_include\$ids), ~ to_include\$X2) %>%
            mutate(phenotype = case_when(is.na(age_astma_diagnosed) ~ "control",
                                         TRUE ~ "asthma")) %>%
            mutate(phenotype = as_factor(phenotype))
        
        write_tsv(all_phenotypes, path = "phenotypes.txt")
        """
}