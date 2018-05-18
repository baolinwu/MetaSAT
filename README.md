# MetaSAT
 - An R package for meta-analysis of variant-set association test

-----
## FE and HE models
 - We study meta-analysis methods that can account for varying levels of
   heterogeneity of genetic effects both within and across studies.
 - Fixed-effects (FE) model assumes similar effects and thus directly
      sums the summary statistics.
 - Heterogeneous-effects (HE) models can have varying assumptions on effects  
   - Within one study, GBT assumes similar effects across variants; SST 
     (SKAT; variance component test) assumes heterogeneous effects across variants;
   - Across studies, HE takes the sum of SST/GBT over studies. SST type tests
     assume different effects across variants/studies; GBT type tests
     assume similar effects across all variants/studies; and SBT assumes
     similar effects across variants within one study while allows for
     different variant-set level effects across studies, which is
     potentially relevant for cross-ancestry meta-analysis.

