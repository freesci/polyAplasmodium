# polyAplasmodium
Scripts and data related to a manuscript on translation of polyA tracks in Plasmodium species.

## System requirements

R language version 3.5.0 or newer with installed libraries:
- rtracklayer
- data.table (DT)
- biomaRt
- genomeIntervals
- ggplot2
and their dependencies.

### Testing environment

R version 3.5.0 (2018-04-23)
Platform: x86_64-apple-darwin17.5.0 (64-bit)
Running under: macOS  10.14

Matrix products: default
BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
LAPACK: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libLAPACK.dylib

locale:
[1] C

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils    
[7] datasets  methods   base     

other attached packages:
[1] ggplot2_3.0.0          genomeIntervals_1.36.0
[3] intervals_0.15.1       biomaRt_2.36.1        
[5] IRanges_2.14.10        S4Vectors_0.18.3      
[7] BiocGenerics_0.26.0   

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.17           bindr_0.1.1           
 [3] pillar_1.2.3           plyr_1.8.4            
 [5] compiler_3.5.0         GenomeInfoDb_1.16.0   
 [7] XVector_0.20.0         prettyunits_1.0.2     
 [9] bitops_1.0-6           tools_3.5.0           
[11] progress_1.2.0         zlibbioc_1.26.0       
[13] digest_0.6.15          bit_1.1-14            
[15] tibble_1.4.2           gtable_0.2.0          
[17] RSQLite_2.1.1          memoise_1.1.0         
[19] pkgconfig_2.0.1        rlang_0.2.1           
[21] DBI_1.0.0              rstudioapi_0.7        
[23] yaml_2.1.19            bindrcpp_0.2.2        
[25] GenomeInfoDbData_1.1.0 withr_2.1.2           
[27] dplyr_0.7.6            stringr_1.3.1         
[29] httr_1.3.1             hms_0.4.2             
[31] tidyselect_0.2.4       grid_3.5.0            
[33] bit64_0.9-7            glue_1.2.0            
[35] Biobase_2.40.0         R6_2.2.2              
[37] AnnotationDbi_1.42.1   XML_3.98-1.11         
[39] purrr_0.2.5            blob_1.1.1            
[41] magrittr_1.5           scales_0.5.0          
[43] GenomicRanges_1.32.3   assertthat_0.2.0      
[45] colorspace_1.3-2       stringi_1.2.3         
[47] lazyeval_0.2.1         munsell_0.5.0         
[49] RCurl_1.95-4.10        crayon_1.3.4    


## Installation instructions

The main script for the paper is called `riboseq_coverage.R` and can be downloaded
from Github together with accompanying data. There's no specific installation instruction.
After fulfilling the library requirements, the script should run without problems provided 
that you set the working directory to the top level directory of this repository.

## Running time and output

Preferred way of running the script is to rely on Rstudio. 
`riboseq_coverage.R` script runs till the end in less than 30 mins on a modern computer with
fast drive. The script does not save any output - plots are generated in the RStudio panel view.

