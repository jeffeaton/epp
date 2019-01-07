context("reference")

test_that("Botswana2017 EPP data is read correctly", {
  pjnz <- system.file("extdata/testpjnz", "Botswana2017.PJNZ", package="epp")
  
  hiv_prevalence <- epp::read_epp_data(pjnz)
  hiv_ref <- readRDS("reference/Botswana2017_hiv_prevalence_data.rds")
  expect_equal(hiv_prevalence, hiv_ref)
  
  population_data <- epp::read_epp_subpops(pjnz)
  population_ref <- readRDS("reference/Botswana2017_population_data.rds")
  expect_equal(population_data, population_ref)
  
  spt_data <- epp::read_spt(pjnz)
  spt_ref <- readRDS("reference/Botswana2017_spt_data.rds")
  expect_equal(spt_data, spt_ref)
})

test_that("Botswana2018 EPP data is read correctly", {
  pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="epp")
  
  hiv_prevalence <- epp::read_epp_data(pjnz)
  hiv_ref <- readRDS("reference/Botswana2018_hiv_prevalence_data.rds")
  expect_equal(hiv_prevalence, hiv_ref)
  
  population_data <- epp::read_epp_subpops(pjnz)
  population_ref <- readRDS("reference/Botswana2018_population_data.rds")
  expect_equal(population_data, population_ref)
  
  spt_data <- epp::read_spt(pjnz)
  spt_ref <- readRDS("reference/Botswana2018_spt_data.rds")
  expect_equal(spt_data, spt_ref)
})

test_that("Mozambique_Maputo_Cidade2018 EPP data is read correctly", {
  pjnz <- system.file("extdata/testpjnz", "Mozambique_Maputo_Cidade2018.PJNZ", 
                      package="epp")
  
  hiv_prevalence <- epp::read_epp_data(pjnz)
  hiv_ref <- readRDS(
    "reference/Mozambique_Maputo_Cidade2018_hiv_prevalence_data.rds")
  expect_equal(hiv_prevalence, hiv_ref)
  
  population_data <- epp::read_epp_subpops(pjnz)
  population_ref <- readRDS(
    "reference/Mozambique_Maputo_Cidade2018_population_data.rds")
  expect_equal(population_data, population_ref)
  
  spt_data <- epp::read_spt(pjnz)
  spt_ref <- readRDS("reference/Mozambique_Maputo_Cidade2018_spt_data.rds")
  expect_equal(spt_data, spt_ref)
})
