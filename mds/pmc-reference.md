## PMC as reference

### Build abstract corpus

> Open Access Common Use PMID list available
> [here](https://ftp.ncbi.nlm.nih.gov/pub/pmc/).

``` r
setwd(ld)
pmc <- read.csv('oa_comm_use_file_list.csv')
pmc0 <- subset(pmc, !is.na(PMID))
```

``` r
corpus <-  PubmedMTK::pmtk_get_records2(pmids = pmc0$PMID, 
                                        cores = 10, 
                                        ncbi_key = key) 

corpus0 <- data.table::rbindlist(corpus)
```

### MeSH annotation frequency

``` r
m0 <- PubmedMTK::pmtk_gather_mesh(x = corpus0)

pmtk_tbl_pmc_ref <- m0[ , list(doc_count = length(unique(pmid))),
           by = list(type, term)]
total_pmid <- length(unique(m0$pmid))
pmtk_tbl_pmc_ref[, doc_prop := doc_count/total_pmid]
```

    ##                  type                             term doc_count     doc_prop
    ##       1: meshHeadings                       algorithms     35621 1.963793e-02
    ##       2: meshHeadings           crystallography,_x-ray      6384 3.519512e-03
    ##       3: meshHeadings data_interpretation,_statistical      3216 1.772987e-03
    ##       4: meshHeadings                 fourier_analysis       472 2.602145e-04
    ##       5: meshHeadings              molecular_structure      8696 4.794122e-03
    ##      ---                                                                     
    ## 1304693: meshHeadings         tricuspid_valve_stenosis         1 5.513020e-07
    ## 1304694: meshHeadings                    ethylestrenol         1 5.513020e-07
    ## 1304695:    chemNames         monoethylglycinexylidide         1 5.513020e-07
    ## 1304696:    chemNames                           savlon         1 5.513020e-07
    ## 1304697:    chemNames            7-propyl_spirolactone         1 5.513020e-07
