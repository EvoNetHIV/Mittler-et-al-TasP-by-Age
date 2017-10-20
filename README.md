# Mittler-et-al-TasP-by-Age


To run example script similar to one used for paper first download evonet and its dependencies:
``` r
if (!require("devtools")) install.packages("devtools")
install_github("EvoNetHIV/TestRepo",subdir="pkg")
install_github( "statnet/tergmLite")
install_github( "statnet/EpiModel", ref ="fast_edgelist")
```

Download from this site: ./scripts/Mittler_et_al_example_script.R and place in working directory (e.g., use getwd() to identify location), then source script:
``` r
source("Mittler_et_al_example_script.R")
```

