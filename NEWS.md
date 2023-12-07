# epp 0.4.4

- Add function `read_epp_workset_options()` to read ANC prevalence adjustment and EPP-ASM model 
  option selections in EPP Curve Fitting page.
  
# epp 0.4.3

- Add function `read_eppxml_results()` to read projection outputs from EPP .xml file including ART coverage.

# epp 0.4.2

- Add type.convert(..., as.is = TRUE) to silence R 4.1 warning.

# epp 0.4.1

- Handle new data structure for household survey data introduced in EPP 2019.

# epp 0.4.0
- Update EPP .xml file parsers to use `xml2` package instead of `XML` package. Improved logic and efficiency for parsers.
- Search XML node trees rather than rely on node position. More robust and code now works to read data from concentrated epidemic structured EPP files
- Parse subpopulation characteristics and turnover from concentrated epidemic files.
- Read incidence data inputs.

# epp 0.3.3

- Only read ANC-RT data if used in EPP fit

# epp 0.3.2

- ~~IMIS offsets likelihood by maximum log-likelihood to avoid underflow.~~  Rolled this back because it introduced a pernicious little bug in mixture weights when offset changes numerical underflow. Should be re-implemented more carefully later.
- Added option to likelihood() to return log values.

# epp 0.3.1

- Updates to IMIS optimization step:
 - Do random sample as well in each iteration
 - Reduce numerical integration stepsize and try again if BFGS encounters non-finite difference error
 - Use last value from center_all to initialise optimization if no positive samples remain

# epp 0.3

- Implement ANC routine testing (ANC-RT) likelihood in EPP 2017 (Spectrum 5.52)
