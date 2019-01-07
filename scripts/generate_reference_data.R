devtools::load_all()

test_cases = c("Botswana2017", "Botswana2018", "Mozambique_Maputo_Cidade2018")
for (test_case in test_cases) {
  print(test_case)
  pjnz <- system.file("extdata/testpjnz",
                      paste(test_case, ".PJNZ", sep = ""),
                      package = "epp")
  epp_data <- read_epp_data(pjnz)
  dir.create("tests/testthat/reference", FALSE, TRUE)
  saveRDS(epp_data, paste("tests/testthat/reference/", 
                          test_case,
                          ".rds",
                          sep = ""))
}
