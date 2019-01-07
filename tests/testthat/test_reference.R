context("reference")

test_that("Botswana2017 EPP data is read correctly", {
  pjnz <- system.file("extdata/testpjnz", "Botswana2017.PJNZ", package="epp")
  
  hiv_prevalence <- epp::read_epp_data(pjnz)
  hiv_ref <- readRDS("reference/Botswana2017_hiv_prevalence_data.rds")
  expect_equal(hiv_prevalence, hiv_ref)
  
  population_data <- epp::read_epp_subpops(pjnz)
  population_ref <- readRDS("reference/Botswana2017_population_data.rds")
  expect_equal(population_data, population_ref)
  
  incidence_prevalence_bf_data <- epp::read_spt(pjnz)
  incidence_prevalence_bf_ref <- readRDS(
    "reference/Botswana2017_incidence_prevalence_bf_data.rds")
  expect_equal(incidence_prevalence_bf_data, incidence_prevalence_bf_ref)
})

test_that("Botswana2018 EPP data is read correctly", {
  pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="epp")
  
  hiv_prevalence <- epp::read_epp_data(pjnz)
  hiv_ref <- readRDS("reference/Botswana2018_hiv_prevalence_data.rds")
  expect_equal(hiv_prevalence, hiv_ref)
  
  population_data <- epp::read_epp_subpops(pjnz)
  population_ref <- readRDS("reference/Botswana2018_population_data.rds")
  expect_equal(population_data, population_ref)
  
  incidence_prevalence_bf_data <- epp::read_spt(pjnz)
  incidence_prevalence_bf_ref <- readRDS(
    "reference/Botswana2018_incidence_prevalence_bf_data.rds")
  expect_equal(incidence_prevalence_bf_data, incidence_prevalence_bf_ref)
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
  
  incidence_prevalence_bf_data <- epp::read_spt(pjnz)
  incidence_prevalence_bf_ref <- readRDS(
    "reference/Mozambique_Maputo_Cidade2018_incidence_prevalence_bf_data.rds")
  expect_equal(incidence_prevalence_bf_data, incidence_prevalence_bf_ref)
})
