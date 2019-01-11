context("reference")

test_that("Botswana2017 EPP data is read correctly", {
  pjnz <- system.file("extdata/testpjnz", "Botswana2017.PJNZ", package="epp")
  
  epp_data <- epp::read_epp_data(pjnz)
  epp_ref <- readRDS("reference/Botswana2017_epp_data.rds")
  expect_equal(epp_data, epp_ref)
  
  epp_subpop_data <- epp::read_epp_subpops(pjnz)
  epp_subpop_ref <- readRDS("reference/Botswana2017_epp_subpop_data.rds")
  expect_equal(epp_subpop_data, epp_subpop_ref)
  
  spt_data <- epp::read_spt(pjnz)
  spt_ref <- readRDS("reference/Botswana2017_spt_data.rds")
  expect_equal(spt_data, spt_ref)
})

test_that("Botswana2018 EPP data is read correctly", {
  pjnz <- system.file("extdata/testpjnz", "Botswana2018.PJNZ", package="epp")
  
  epp_data <- epp::read_epp_data(pjnz)
  epp_ref <- readRDS("reference/Botswana2018_epp_data.rds")
  expect_equal(epp_data, epp_ref)
  
  epp_subpop_data <- epp::read_epp_subpops(pjnz)
  epp_subpop_ref <- readRDS("reference/Botswana2018_epp_subpop_data.rds")
  expect_equal(epp_subpop_data, epp_subpop_ref)
  
  spt_data <- epp::read_spt(pjnz)
  spt_ref <- readRDS("reference/Botswana2018_spt_data.rds")
  expect_equal(spt_data, spt_ref)
})

test_that("Mozambique_Maputo_Cidade2018 EPP data is read correctly", {
  pjnz <- system.file("extdata/testpjnz", "Mozambique_Maputo_Cidade2018.PJNZ", 
                      package="epp")
  
  epp_data <- epp::read_epp_data(pjnz)
  epp_ref <- readRDS("reference/Mozambique_Maputo_Cidade2018_epp_data.rds")
  expect_equal(epp_data, epp_ref)
  
  epp_subpop_data <- epp::read_epp_subpops(pjnz)
  epp_subpop_ref <- readRDS(
    "reference/Mozambique_Maputo_Cidade2018_epp_subpop_data.rds")
  expect_equal(epp_subpop_data, epp_subpop_ref)
  
  spt_data <- epp::read_spt(pjnz)
  spt_ref <- readRDS("reference/Mozambique_Maputo_Cidade2018_spt_data.rds")
  expect_equal(spt_data, spt_ref)
})

test_that("DominicanRepublic2017 EPP data is read correctly", {
  pjnz <- system.file("extdata/testpjnz", "DominicanRepublic2017.PJNZ", 
                      package="epp")
  
  epp_data <- epp::read_epp_data(pjnz)
  epp_ref <- readRDS("reference/DominicanRepublic2017_epp_data.rds")
  expect_equal(epp_data, epp_ref)
  
  epp_subpop_data <- epp::read_epp_subpops(pjnz)
  epp_subpop_ref <- readRDS(
    "reference/DominicanRepublic2017_epp_subpop_data.rds")
  expect_equal(epp_subpop_data, epp_subpop_ref)
  
  spt_data <- epp::read_spt(pjnz)
  spt_ref <- readRDS("reference/DominicanRepublic2017_spt_data.rds")
  expect_equal(spt_data, spt_ref)
})
