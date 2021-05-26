# PubmedMTK

***PubMed Mining Toolkit***

An R package for querying the PubMed database & parsing retrieved
records. Toolkit facilitates batch API requests, the creation of custom
corpora for NLP, and the quick exploration & visualization of topic
structure.

-   [Installation](#installation)
-   [Usage](#usage)
    -   [MeSH vocabulary](#mesh-vocabulary)
    -   [PubMed search](#pubmed-search)
    -   [Extract MeSH classifications](#extract-mesh-classifications)
    -   [MeSH annotations-based topic
        model](#mesh-annotations-based-topic-model)
    -   [Two-dimensional analyses](#two-dimensional-analyses)

## Installation

``` r
devtools::install_github("jaytimm/PubmedMTK")
```

## Usage

### MeSH vocabulary

The package includes as a data frame the MeSH thesaurus &
hierarchically-organized vocabulary – comprised of 2021 versions of
`descriptor` & `trees` files made available via NLM-NIH. [A workflow for
re-creating the table from raw data
sets](https://github.com/jaytimm/PubmedMTK/blob/main/mds/build-MeSH-df.md).

### PubMed search

The `pmtk_search_pubmed()` function is meant for record-matching
searches typically performed using the [PubMed online
interface](https://pubmed.ncbi.nlm.nih.gov/). The `search_term`
parameter specifies the query term; the `fields` parameter can be used
to specify which fields to query.

``` r
s0 <- PubmedMTK::pmtk_search_pubmed(search_term = 'marijuana', fields = c('TIAB','MH'))
```

    ## [1] "marijuana: 21962 records"

``` r
head(s0)
```

    ##    search_term     pmid
    ## 1:   marijuana 34033378
    ## 2:   marijuana 34032489
    ## 3:   marijuana 34028895
    ## 4:   marijuana 34026421
    ## 5:   marijuana 34022843
    ## 6:   marijuana 34018308

### Retrieve abstract data from PubMed

``` r
sen_df <- PubmedMTK::pmtk_get_records2(pmids = s0$pmid, 
                                       cores = 6, 
                                       ncbi_key = key) 

sen_df <- data.table::rbindlist(sen_df)
```

``` r
colnames(sen_df)
```

    ## [1] "pmid"         "year"         "articletitle" "meshHeadings" "chemNames"   
    ## [6] "keywords"     "abstract"

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

    ## INFO  [15:49:43.128] early stopping at 180 iteration 
    ## INFO  [15:49:44.619] early stopping at 20 iteration

``` r
tm_summary <- PubmedMTK::pmtk_summarize_lda(
  lda = lda, topic_feats_n = 10)
```

#### Feature composition of first ten topics

| topic_id | topic_features                                                                                                                                                                            |
|---------:|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|        1 | male \| adult \| marijuana_abuse \| female \| mental_disorders \| cannabis \| psychoses,\_substance-induced \| schizophrenia \| psychotic_disorders \| psychiatric_status_rating_scales   |
|        2 | substance-related_disorders \| united_states \| male \| adolescent \| female \| age_factors \| adult \| attitude \| sex_factors \| data_collection                                        |
|        3 | adolescent \| adolescent_behavior \| child \| substance-related_disorders \| longitudinal_studies \| peer_group \| risk_factors \| social_environment \| violence \| juvenile_delinquency |
|        4 | male \| adult \| female \| treatment_outcome \| depression \| middle_aged \| anxiety \| cannabis \| analgesics,\_opioid \| opioid-related_disorders                                       |
|        5 | adolescent \| male \| students \| surveys_and_questionnaires \| female \| alcohol_drinking \| motivation \| universities \| young_adult \| marijuana_smoking                              |
|        6 | marijuana_abuse \| female \| male \| adult \| follow-up_studies \| alcoholism \| risk_factors \| longitudinal_studies \| cohort_studies \| young_adult                                    |
|        7 | young_adult \| marijuana_use \| substance_use \| sexual_behavior \| hiv_infections \| risk-taking \| alcohol \| male \| adolescents \| adolescent                                         |
|        8 | animals \| dronabinol \| rats \| time_factors \| dose-response_relationship,\_drug \| behavior,\_animal \| brain \| motor_activity \| mice \| drug_tolerance                              |
|        9 | female \| pregnancy \| adult \| child \| prenatal_exposure_delayed_effects \| infant,\_newborn \| cannabis \| prospective_studies \| smoking \| pregnancy_complications                   |
|       10 | medical_marijuana \| phytotherapy \| united_states \| legislation,\_drug \| drug_and_narcotic_control \| cannabis \| pain \| canada \| public_health \| vomiting                          |

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

#### Principal component analysis (PCA)

``` r
two_ds$pca %>%
  ggplot(aes(x = X1,
             y = X2,
             label = topic_id)) +
  ggplot2::geom_point(size = 10, 
                      color = '#a5c8e1',
                      alpha = 0.5) +
  geom_text(size = 3) +
  ggtitle('Topics in 2d PCA space')
```

![](README_files/figure-markdown_github/unnamed-chunk-17-1.png)
