# PubmedMTK

***PubMed Mining Toolkit***

An R package for querying the PubMed database & parsing retrieved
records. Toolkit facilitates batch API requests, the creation of custom
corpora for NLP, and the quick exploration & visualization of topic
structure.

-   [Installation](#installation)
-   [Usage](#usage)
    -   [PubMed search](#pubmed-search)
    -   [Extract MeSH classifications](#extract-mesh-classifications)
    -   [MeSH annotations-based topic
        model](#mesh-annotations-based-topic-model)
    -   [Two-dimensional analyses](#two-dimensional-analyses)
-   [Tables](#tables)
    -   [MeSH vocabulary](#mesh-vocabulary)
    -   [PMC MeSH annotation
        frequencies](#pmc-mesh-annotation-frequencies)
    -   [Medline citation counts](#medline-citation-counts)

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

    ## [1] "medical marijuana: 2125 records"

``` r
head(s0)
```

    ##          search_term     pmid
    ## 1: medical marijuana 34044753
    ## 2: medical marijuana 34007062
    ## 3: medical marijuana 33998880
    ## 4: medical marijuana 33998864
    ## 5: medical marijuana 33981161
    ## 6: medical marijuana 33974499

### Retrieve abstract data from PubMed

``` r
sen_df <- PubmedMTK::pmtk_get_records2(pmids = s0$pmid, 
                                       cores = 6, 
                                       ncbi_key = key) 

sen_df <- data.table::rbindlist(sen_df)
```

> Sample record from output:

``` r
n <- 9
list(pmid = sen_df$pmid[n],
     year = sen_df$year[n],
     articletitle = strwrap(sen_df$articletitle[n], width = 60),
     meshHeadings = strwrap(sen_df$meshHeadings[n], width = 60),
     text = strwrap(sen_df$abstract[n], width = 60)[1:10])
```

    ## $pmid
    ## [1] "33933061"
    ## 
    ## $year
    ## [1] "2021"
    ## 
    ## $articletitle
    ## [1] "Opioid use in medical cannabis authorization adult patients"
    ## [2] "from 2013 to 2018: Alberta, Canada."                        
    ## 
    ## $meshHeadings
    ## [1] "Adult|Alberta|Analgesics,"                                  
    ## [2] "Opioid|Cannabis|Female|Humans|Male|Medical Marijuana|Middle"
    ## [3] "Aged|Opioid-Related Disorders|United States"                
    ## 
    ## $text
    ##  [1] "The opioid overdose epidemic in Canada and the United"      
    ##  [2] "States has become a public health crisis - with exponential"
    ##  [3] "increases in opioid-related morbidity and mortality."       
    ##  [4] "Recently, there has been an increasing body of evidence"    
    ##  [5] "focusing on the opioid-sparing effects of medical cannabis" 
    ##  [6] "use (reduction of opioid use and reliance), and medical"    
    ##  [7] "cannabis as a potential alternative treatment for chronic"  
    ##  [8] "pain. The objective of this study is to assess the effect"  
    ##  [9] "of medical cannabis authorization on opioid use (oral"      
    ## [10] "morphine equivalent; OME) between 2013 and 2018 in Alberta,"

### Extract MeSH classifications

Subject terms/headings in metadata table include `MeSH` terms, as well
as (some) `keywords` & `chem-names`. The `pmtk_gather_mesh` function
extracts & structures these attributes from metadata. The resulting
table amounts to a document-term matrix (DTM), in which each PubMed
abstract is represented as a vector of MeSH terms.

``` r
m0 <- PubmedMTK::pmtk_gather_mesh(sen_df)
```

### MeSH annotations-based topic model

We can use these MeSH-based abstract representations to explore the
conceptual structure of a particular collection of PubMed records via
topic modeling. Here we implement **Latent Dirichlet allocation**, which
is a topic modeling algorithm that models *each document* in corpus as a
composite of topics, and *each topic* as a composite of terms. Topic
composition can be interpreted as sets of MeSH terms that frequently
co-occur.

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
composition from the `text2vec::LDA` output. For each possible
topic-feature pair, the model computes the likelihood a given topic
generated a given feature. Output is filtered to the highest scoring
features per topic using the `topic_feats_n`.

``` r
lda <- text2vec::LDA$new(n_topics = 20) 
fit <- lda$fit_transform(dtm, progressbar = F)
```

    ## INFO  [10:17:14.068] early stopping at 90 iteration 
    ## INFO  [10:17:14.271] early stopping at 30 iteration

``` r
tm_summary <- PubmedMTK::pmtk_summarize_lda(
  lda = lda, topic_feats_n = 10)
```

#### Feature composition of first ten topics

| topic_id | topic_features                                                                                                                                                                                                             |
|---------:|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|        1 | adolescent \| marijuana_smoking \| cannabis \| california \| young_adult \| canada \| adolescent_behavior \| substance_use \| commerce \| pharmacists                                                                      |
|        2 | cannabis \| canada \| phytotherapy \| marijuana \| illicit_drugs \| receptor,\_cannabinoid,\_cb1 \| cross-sectional_studies \| endocannabinoid_system \| addiction \| receptor,\_cannabinoid,\_cb2                         |
|        3 | male \| adult \| female \| middle_aged \| surveys_and_questionnaires \| internet \| self_report \| public_opinion \| accidents,\_traffic \| risk_factors                                                                   |
|        4 | cannabidiol \| cannabis \| epilepsy \| dronabinol \| treatment_outcome \| quality_of_life \| plant_extracts \| randomized_controlled_trials_as_topic \| neuralgia \| germany                                               |
|        5 | drug_and_narcotic_control \| united_states \| state_government \| government_regulation \| federal_government \| drug_approval \| health_policy \| marijuana_smoking \| risk_assessment \| legislation,\_medical           |
|        6 | cannabis \| cannabinoids \| attitude_of_health_personnel \| neoplasms \| health_knowledge,\_attitudes,\_practice \| evidence-based_medicine \| risk_assessment \| practice_guidelines_as_topic \| health_personnel \| pain |
|        7 | marijuana \| legislation,\_drug \| united_states \| child \| public_policy \| marijuana_use \| cannabis \| legalization \| commerce \| public_health                                                                       |
|        8 | adult \| aged \| middle_aged \| cross-sectional_studies \| surveys_and_questionnaires \| young_adult \| aged,\_80_and_over \| female \| male \| health_knowledge,\_attitudes,\_practice                                    |
|        9 | legislation,\_drug \| substance-related_disorders \| practice_patterns,\_physicians’ \| phytotherapy \| marijuana_use \| physicians \| politics \| united_states \| analgesics,\_opioid \| socioeconomic_factors           |
|       10 | cannabinoids \| chronic_pain \| medical_cannabis \| pain_management \| israel \| australia \| chronic_disease \| stress_disorders,\_post-traumatic \| fibromyalgia \| cannabinoid_receptor_modulators                      |

### Two-dimensional analyses

``` r
tmat <- tidytext::cast_sparse(data = tm_summary$topic_word_dist,
                              row = topic_id,
                              column = feature,
                              value = beta)

## build 2d data structures --
two_ds <- PubmedMTK::pmtk_2d(mat = tmat)
```

#### Hierarchical clustering

``` r
#two_ds$hc$labels <- tm_summary$topic_summary$topic_features
two_ds$hc %>% ggdendro::ggdendrogram(rotate=TRUE)
```

![](README_files/figure-markdown_github/unnamed-chunk-16-1.png)

#### tSNE

``` r
two_ds$tsne %>%
  ggplot(aes(x = X1,
             y = X2,
             label = topic_id)) +
  ggplot2::geom_point(size = 10, 
                      color = '#a5c8e1',
                      alpha = 0.5) +
  geom_text(size = 3) +
  ggtitle('Topics in 2d space')
```

![](README_files/figure-markdown_github/unnamed-chunk-17-1.png)

## Tables

### MeSH vocabulary

The package includes as a data frame the MeSH thesaurus &
hierarchically-organized vocabulary – comprised of 2021 versions of
`descriptor` & `trees` files made available via NLM-NIH. [A
workflow](https://github.com/jaytimm/PubmedMTK/blob/main/mds/build-MeSH-df.md)
for re-creating the table from raw data sets.

``` r
PubmedMTK::pmtk_tbl_mesh[1:5, c(1:3, 5:6)]
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
PubmedMTK::pmtk_tbl_pmc_ref
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

### Medline citation counts

Table build details available
[here](https://github.com/jaytimm/PubmedMTK/blob/main/mds/medline_citations.md).

``` r
tail(PubmedMTK::pmtk_tbl_citations)
```

    ##          year  total    usa
    ## 1: 2016-01-01 862744 350671
    ## 2: 2017-01-01 848229 343161
    ## 3: 2018-01-01 857322 339199
    ## 4: 2019-01-01 865068 334305
    ## 5: 2020-01-01 871672 329411
    ## 6: 2021-01-01 877300 324518
