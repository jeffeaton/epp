context("reference")

test_that("Botswana2017 epp data is read correctly", {
  pjnz <- system.file("extdata/testpjnz", "Botswana2017.PJNZ", package="epp")
  epp_data <- epp::read_epp_data(pjnz)
  ref <- readRDS("reference/Botswana2017.rds")
  expect_equal(epp_data, ref)
})

test_that("Botswana2018 epp data is read correctly", {
  pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="epp")
  epp_data <- epp::read_epp_data(pjnz)
  ref <- readRDS("reference/Botswana2018.rds")
  expect_equal(epp_data, ref)
})

test_that("Mozambique_Maputo_Cidade2018 epp data is read correctly", {
  pjnz <- system.file("extdata/testpjnz", "Mozambique_Maputo_Cidade2018.PJNZ", 
                      package="epp")
  epp_data <- epp::read_epp_data(pjnz)
  ref <- readRDS("reference/Mozambique_Maputo_Cidade2018.rds")
  expect_equal(epp_data, ref)
})
