# epp 0.3

- Implement ANC routine testing (ANC-RT) likelihood in EPP 2017 (Spectrum 5.52)

# epp 0.3.1

- Updates to IMIS optimization step:
 - Do random sample as well in each iteration
 - Reduce numerical integration stepsize and try again if BFGS encounters non-finite difference error
 - Use last value from center_all to initialise optimization if no positive samples remain