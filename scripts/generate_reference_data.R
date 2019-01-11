devtools::load_all()

test_cases = c("Botswana2017", "Botswana2018", "Mozambique_Maputo_Cidade2018",
               "DominicanRepublic2017")
for (test_case in test_cases) {
  print(test_case)
  dir.create("tests/testthat/reference", FALSE, TRUE)
  pjnz_path <- system.file("extdata/testpjnz",
                           paste(test_case, ".PJNZ", sep = ""),
                           package = "epp")
  
  epp_data <- read_epp_data(pjnz_path)
  saveRDS(epp_data, paste("tests/testthat/reference/", 
                                test_case,  
                                "_epp_data.rds", 
                                sep = ""))
  
  
  epp_subpop_data <- read_epp_subpops(pjnz_path)
  saveRDS(epp_subpop_data, paste("tests/testthat/reference/", 
                                 test_case,
                                 "_epp_subpop_data.rds",
                                 sep = ""))
  
  spt_data <- read_spt(pjnz_path)
  saveRDS(spt_data, paste("tests/testthat/reference/", 
                                   test_case,
                                   "_spt_data.rds",
                                   sep = ""))
  
  # spu_data <- read_spu(pjnz_path)
  # saveRDS(spu_data, paste("tests/testthat/reference/", 
  #                         test_case,
  #                         "_spu_data.rds",
  #                         sep = ""))
}
