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
