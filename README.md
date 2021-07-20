<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/jaytimm/PubmedMTK.svg?branch=main)](https://travis-ci.com/jaytimm/PubmedMTK)
<!-- badges: end -->

# PubmedMTK

An R package for querying the PubMed database & parsing retrieved
records. Toolkit facilitates batch API requests & the creation of custom
corpora for NLP.

-   [Installation](#installation)
-   [Usage](#usage)
    -   [PubMed search](#pubmed-search)
    -   [Retrieve and parse abstract
        data](#retrieve-and-parse-abstract-data)
    -   [Sentence tokenization](#sentence-tokenization)
    -   [KWIC search](#kwic-search)
    -   [Extract MeSH classifications](#extract-mesh-classifications)
    -   [MeSH annotations-based topic
        model](#mesh-annotations-based-topic-model)
-   [Interactive HTML topic summary](#interactive-html-topic-summary)
-   [Tables](#tables)
    -   [MeSH vocabulary](#mesh-vocabulary)
    -   [PMC MeSH annotation
        frequencies](#pmc-mesh-annotation-frequencies)

## Installation

``` r
devtools::install_github("jaytimm/PubmedMTK")
```

## Usage

### PubMed search

The `pmtk_search_pubmed()` function is meant for record-matching
searches typically performed using the [PubMed online
interface](https://pubmed.ncbi.nlm.nih.gov/). The `search_term`
parameter specifies the query term; the `fields` parameter can be used
to specify which fields to query.

``` r
s0 <- PubmedMTK::pmtk_search_pubmed(search_term = 'medical marijuana', 
                                    fields = c('TIAB','MH'))
```

    ## [1] "medical marijuana[TIAB] OR medical marijuana[MH]: 2177 records"

``` r
head(s0)
```

    ##          search_term     pmid
    ## 1: medical marijuana 34258798
    ## 2: medical marijuana 34258360
    ## 3: medical marijuana 34234903
    ## 4: medical marijuana 34232573
    ## 5: medical marijuana 34225825
    ## 6: medical marijuana 34179729

### Retrieve and parse abstract data

For quicker abstract retrieval, be sure to get an [API
key](https://support.nlm.nih.gov/knowledgebase/article/KA-03521/en-us).

``` r
sen_df <- PubmedMTK::pmtk_get_records2(pmids = s0$pmid, 
                                       cores = 6, 
                                       ncbi_key = key) 
```

> Sample record from output:

``` r
sen_df <- data.table::rbindlist(sen_df)

n <- 10
list(pmid = sen_df$pmid[n],
     year = sen_df$year[n],
     articletitle = strwrap(sen_df$articletitle[n], width = 60),
     meshHeadings = strwrap(sen_df$meshHeadings[n], width = 60),
     text = strwrap(sen_df$abstract[n], width = 60)[1:10])
```

    ## $pmid
    ## [1] "34128629"
    ## 
    ## $year
    ## [1] "2021"
    ## 
    ## $articletitle
    ## [1] "Integrative Medicine: Cannabis and Cannabis-Related Drugs."
    ## 
    ## $meshHeadings
    ## [1] "Cannabidiol|Cannabis|Child|Humans|Integrative"         
    ## [2] "Medicine|Medical Marijuana|Pharmaceutical Preparations"
    ## 
    ## $text
    ##  [1] "Cannabis is a genus of flowering herbs in the Cannabaceae"  
    ##  [2] "family. Federal law defines dried plant material"           
    ##  [3] "preparations of the subspecies Cannabis sativa as"          
    ##  [4] "marijuana. The term cannabis refers to all products derived"
    ##  [5] "from Cannabis plants. The active compounds in cannabis are" 
    ##  [6] "cannabinoids, which include delta-9-tetrahydrocannabinol"   
    ##  [7] "(THC) and cannabidiol (CBD). THC is the psychoactive"       
    ##  [8] "component, whereas CBD has no psychoactive effects. There"  
    ##  [9] "are three Food and Drug Administration (FDA)-approved"      
    ## [10] "cannabis-related drugs. Dronabinol and nabilone (Cesamet)"

### Sentence tokenization

``` r
sentences <- PubmedMTK::pmtk_toke_sentences(text = sen_df$abstract,
                                            doc_id = sen_df$pmid)
```

| pmid       | abstract                                                                                                                                                                                                                                                            |
|:-----------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 34258798.1 | Although cannabis is federally prohibited, a majority of U.S. states have implemented medical cannabis laws (MCLs).                                                                                                                                                 |
| 34258798.2 | As more individuals consider the drug for medical treatment, they potentially substitute away from prescription drugs.                                                                                                                                              |
| 34258798.3 | Therefore, an MCL signals competitor entry.                                                                                                                                                                                                                         |
| 34258798.4 | This paper exploits geographic and temporal variation in MCLs to examine the strategic response in direct-to-physician marketing by pharmaceutical firms as cannabis enters the market.                                                                             |
| 34258798.5 | Using office detailing records from 2014-2018 aggregated to the county level, we find weak evidence of a relatively small and delayed response in substitute prescription drug- and opioid-related detailing.                                                       |
| 34258798.6 | While these effects on detailing dollars are more pronounced among smaller pharmaceutical firms, the magnitudes are economically small and likely muted at aggregate levels by the small percent of doctors that actively recommend cannabis for medical treatment. |

### KWIC search

The `pmtk_locate_term()` function allows for quick keyword-in-context
(KWIC) search. A simple wrapper of the `corpus::text_locate` function.

``` r
toks <-  corpus::text_tokens(sen_df$abstract)

egs <- PubmedMTK::pmtk_locate_term(text = toks,
                                   doc_id = sen_df$pmid,
                                   term = c('medical marijuana laws'),
                                   stem = F,
                                   window = 10)

egs$kwic <- paste0('... ', egs$lhs, ' `', egs$instance, '` ', egs$rhs, ' ...')
knitr::kable(egs[1:8, c(1,5)])
```

| doc_id   | kwic                                                                                                                                                                       |
|:---------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 34128629 | … , and driving . physicians can recommend use of marijuana under `medical marijuana laws` but cannot prescribe it , as it is classified as a …                            |
| 33750275 | … moving to reverse marijuana prohibition , most frequently through legalization of `medical marijuana laws` ( mmls ) , and there is concern that marijuana legalization … |
| 33730400 | … recreational marijuana laws ( rml ) , followed by states with `medical marijuana laws` ( mml ) and without legal cannabis use , respectively . …                         |
| 33624387 | … differences-in-differences ( dd ) approach and found that the implementation of `medical marijuana laws` ( mmls ) and recreational marijuana laws ( rmls ) reduced …     |
| 33143941 | … . cannabis legalization was determined by the presence or absence of `medical marijuana laws` ( mml ) and recreational marijuana laws ( rml ) in …                       |
| 33069561 | … , a limited but growing body of literature has found state `medical marijuana laws` ( mmls ) to be associated with lower levels of opioid …                              |
| 32799573 | … management , much research has since focused on the potential for `medical marijuana laws` ( mmls ) to curb the opioid epidemic . nonetheless , …                        |
| 32736294 | … cannabis use disorder are more prevalent in u.s . states with `medical marijuana laws` ( mmls ) , as well as among individuals with elevated …                           |

### Extract MeSH classifications

Subject terms/headings in metadata table include `MeSH` terms, as well
as (some) `keywords` & `chem-names`. The `pmtk_gather_mesh` function
extracts & structures these attributes from metadata.

``` r
m0 <- PubmedMTK::pmtk_gather_mesh(sen_df)
```

### MeSH annotations-based topic model

We can use these MeSH-based abstract representations to explore the
conceptual structure of a particular collection of PubMed records via
topic modeling. Here we implement **Latent Dirichlet allocation**, which
is a topic modeling algorithm that models *each document* in corpus as a
composite of topics, and *each topic* as a composite of terms.

``` r
as.text <- m0[, list(text = paste(term, collapse = " ")), by = pmid]
iter <- text2vec::itoken(as.text$text, ids = as.text$pmid)  
vocab <- text2vec::create_vocabulary(iter)

vocab0 <- text2vec::prune_vocabulary(
  vocab, 
  # doc_proportion_min = 0.0001,
  doc_proportion_max = 0.55,
  doc_count_min = 3) 

vectorizer <- text2vec::vocab_vectorizer(vocab0)
dtm <- text2vec::create_dtm(iter, vectorizer)
```

The `pmtk_summarize_lda` function summarizes and extracts topic
composition from the `text2vec::LDA` output.

``` r
lda <- text2vec::LDA$new(n_topics = 20) 
fit <- lda$fit_transform(dtm, progressbar = F)
```

    ## INFO  [18:54:40.800] early stopping at 120 iteration 
    ## INFO  [18:54:41.156] early stopping at 30 iteration

``` r
tm_summary <- PubmedMTK::pmtk_summarize_lda(
  lda = lda, topic_feats_n = 10)
```

#### Feature composition of first ten topics

| topic_id | topic_features                                                                                                                                                                                                                 |
|---------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|        1 | marijuana \| marijuana_use \| united_states \| female \| adolescent_behavior \| glaucoma \| stress_disorders,\_post-traumatic \| headache \| emergency_service,\_hospital \| adolescents                                       |
|        2 | dronabinol \| nausea \| multiple_sclerosis \| vomiting \| nabiximols \| plant_extracts \| drug_combinations \| cannabis \| pain \| analgesics                                                                                  |
|        3 | cannabinoids \| animals \| epilepsy \| cannabidiol \| endocannabinoids \| anticonvulsants \| treatment_outcome \| brain \| receptors,\_cannabinoid \| mice                                                                     |
|        4 | neoplasms \| chronic_pain \| analgesics \| quality_of_life \| treatment_outcome \| pain_management \| cannabinoids \| cancer_pain \| cancer \| medicinal_cannabis                                                              |
|        5 | legislation,\_drug \| medical_cannabis \| public_health \| colorado \| legalization \| marijuana_smoking \| government_regulation \| public_policy \| health_knowledge,\_attitudes,\_practice \| epidemiology                  |
|        6 | female \| male \| middle_aged \| adult \| cross-sectional_studies \| self_report \| socioeconomic_factors \| motivation \| los_angeles \| longitudinal_studies                                                                 |
|        7 | cannabis \| cannabinoids \| biomedical_research \| hallucinogens \| canada \| san_francisco \| inflammatory_bowel_diseases \| phytotherapy \| research \| occupational_health                                                  |
|        8 | attitude_of_health_personnel \| health_policy \| risk_assessment \| marijuana_smoking \| drug_approval \| treatment_outcome \| australia \| cannabis \| minnesota \| new_zealand                                               |
|        9 | male \| aged \| middle_aged \| adult \| female \| aged,\_80_and_over \| young_adult \| health_knowledge,\_attitudes,\_practice \| surveys_and_questionnaires \| prospective_studies                                            |
|       10 | marijuana_smoking \| canada \| drug_and_narcotic_control \| mental_disorders \| legislation,\_drug \| palliative_care \| health_services_accessibility \| united_kingdom \| physician-patient_relations \| societies,\_medical |

## Interactive HTML topic summary

``` r
tmat <- tidytext::cast_sparse(data = tm_summary$topic_word_dist,
                              row = topic_id,
                              column = feature,
                              value = beta)

set.seed(99)
tsne <- Rtsne::Rtsne(X = as.matrix(tmat), 
                     check_duplicates = T,
                     perplexity = 5)

tsne0 <- data.frame(topic_id = as.integer(rownames(tmat)), tsne$Y)
```

The `pmtk_build_interactive()` will generate an interactive html file
(built on the `flexdashboard` and `crosstalk` packages) that summarizes
(1) topics in two-dimensional, tSNE space and (2) topic composition. An
image of the html is presented below.

``` r
PubmedMTK::pmtk_build_interactive(pmtk_lda = tm_summary,
                                  pmtk_2d = tsne0,
                                  out_dir = '/home/jtimm/Desktop/',
                                  file_name = 'party.html')
```

<img src="README_files/figure-markdown_github/demo.png" width="100%" />

## Tables

### MeSH vocabulary

The package includes as a data frame the MeSH thesaurus &
hierarchically-organized vocabulary – comprised of 2021 versions of
`descriptor` & `trees` files made available via NLM-NIH. [A
workflow](https://github.com/jaytimm/PubmedMTK/blob/main/mds/build-MeSH-df.md)
for re-creating the table from raw data sets.

``` r
library(PubmedMTK)
data(pmtk_tbl_mesh)
pmtk_tbl_mesh[1:5, c(1:3, 5:6)]
```

    ##    DescriptorUI DescriptorName          TermName                cats
    ## 1:      D000001     calcimycin        calcimycin Chemicals and Drugs
    ## 2:      D000001     calcimycin           a-23187 Chemicals and Drugs
    ## 3:      D000001     calcimycin           a 23187 Chemicals and Drugs
    ## 4:      D000001     calcimycin            a23187 Chemicals and Drugs
    ## 5:      D000001     calcimycin antibiotic a23187 Chemicals and Drugs
    ##                     mesh1
    ## 1: Heterocyclic Compounds
    ## 2: Heterocyclic Compounds
    ## 3: Heterocyclic Compounds
    ## 4: Heterocyclic Compounds
    ## 5: Heterocyclic Compounds

### PMC MeSH annotation frequencies

MeSH annotation frequencies for the Open Access Common Use portion of
PMC. Frequencies based on roughly 1.8 million PubMed records. Details
[here](https://github.com/jaytimm/PubmedMTK/blob/main/mds/pmc-reference.md).

``` r
data(pmtk_tbl_pmc_ref)
pmtk_tbl_pmc_ref
```

    ##                  type                   DescriptorName doc_count     doc_prop
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
