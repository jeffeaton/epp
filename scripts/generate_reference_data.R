devtools::load_all()

test_cases = c("Botswana2017", "Botswana2018", "Mozambique_Maputo_Cidade2018")
for (test_case in test_cases) {
  print(test_case)
  dir.create("tests/testthat/reference", FALSE, TRUE)
  pjnz_path <- system.file("extdata/testpjnz",
                           paste(test_case, ".PJNZ", sep = ""),
                           package = "epp")
  
  hiv_prevalence <- read_epp_data(pjnz_path)
  saveRDS(hiv_prevalence, paste("tests/testthat/reference/", 
                                test_case,  
                                "_hiv_prevalence_data.rds", 
                                sep = ""))
  
  
  population_data <- read_epp_subpops(pjnz_path)
  saveRDS(population_data, paste("tests/testthat/reference/", 
                                 test_case,
                                 "_population_data.rds",
                                 sep = ""))
  
  incidence_prevalence_bf <- read_spt(pjnz_path)
  saveRDS(incidence_prevalence_bf, paste("tests/testthat/reference/", 
                                   test_case,
                                   "_incidence_prevalence_bf_data.rds",
                                   sep = ""))
  
  # spu_data <- read_spu(pjnz_path)
  # saveRDS(spu_data, paste("tests/testthat/reference/", 
  #                         test_case,
  #                         "_spu_data.rds",
  #                         sep = ""))
}
