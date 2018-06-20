# MetaSAT
 - An R package for meta-analysis of variant-set association test

-----
## Meta-analysis models
 - We study meta-analysis methods that can account for varying levels of
     heterogeneity of genetic effects both within and across studies.
   - Burden type test (BT): similar variant effects
   - Variance component test (VT): heterogeneous variant effects
   - Adaptive test (AT): adaptively combine BT and VT
 - Fixed-effects (FE) model 
   - Assumes homogeneous effects across studies.
 - Heterogeneous-effects (HE) model
   - heterogeneous effects across studies and across variants within each study
 - Robust heterogeneous-effects (RHE) model
   - heterogeneous effects across studies but homogeneous effects across variants within each study
 - In summary, a total four general tests
   - FE BT: homogeneous effects across and within studies.
   - FE VT: homogeneous effects across studies/heterogeneous effects within studies
   - HE VT: heterogeneous effects across and within studies.
   - RHE BT: heterogeneous effects across studies/homogeneous effects within studies
 - And four ATs implemented
   - FE AT (FAT): combine FE BT and FE VT 
   - HE AT (HAT): combine FE BT and HE VT
   - RHE AT (RAT): combine RHE BT and HE VT
   - BT based AT (BAT): combine FE BT and RHE BT

 - Ref
   - Wu,B. and Zhao,H. (2018). Efficient and powerful meta-analysis of variant-set association tests using MetaSAT. tech rep. 

## Sample R codes
 - The MetaSAT R package also contains the summary statistics and their associated asymptotic covariance matrices 
for 9 variants across 3 studies as discussed in Section 4.
They can be used to reproduce the analysis results in the manuscript.
```R
 data(WZda)
 Us = WZda$Us; Vs = WZda$Vs
 FMSAT(Us,Vs)
 HMSAT(Us,Vs)
 RMSAT(Us,Vs)
 RBAT(Us,Vs)
```