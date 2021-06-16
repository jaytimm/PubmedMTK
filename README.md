# PubmedMTK

***PubMed Mining Toolkit***

An R package for querying the PubMed database & parsing retrieved
records. Toolkit facilitates batch API requests, the creation of custom
corpora for NLP, and the quick exploration & visualization of topic
structure.

-   [Installation](#installation)
-   [Usage](#usage)
    -   [PubMed search](#pubmed-search)
    -   [Retrieve and parse abstract
        data](#retrieve-and-parse-abstract-data)
    -   [KWIC search](#kwic-search)
    -   [Extract MeSH classifications](#extract-mesh-classifications)
    -   [MeSH annotations-based topic
        model](#mesh-annotations-based-topic-model)
    -   [Two-dimensional analyses](#two-dimensional-analyses)
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

    ## [1] "medical marijuana[TIAB] OR medical marijuana[MH]: 2134 records"

``` r
head(s0)
```

    ##          search_term     pmid
    ## 1: medical marijuana 34128629
    ## 2: medical marijuana 34109050
    ## 3: medical marijuana 34082823
    ## 4: medical marijuana 34044753
    ## 5: medical marijuana 34007062
    ## 6: medical marijuana 33998880

### Retrieve and parse abstract data

``` r
sen_df <- PubmedMTK::pmtk_get_records2(pmids = s0$pmid, 
                                       cores = 6, 
                                       ncbi_key = key) 
```

> Sample record from output:

``` r
sen_df <- data.table::rbindlist(sen_df)

n <- 9
list(pmid = sen_df$pmid[n],
     year = sen_df$year[n],
     articletitle = strwrap(sen_df$articletitle[n], width = 60),
     meshHeadings = strwrap(sen_df$meshHeadings[n], width = 60),
     text = strwrap(sen_df$abstract[n], width = 60)[1:10])
```

    ## $pmid
    ## [1] "33974499"
    ## 
    ## $year
    ## [1] "2021"
    ## 
    ## $articletitle
    ## [1] "The Risk of QTc Prolongation with Antiemetics in the"
    ## [2] "Palliative Care Setting: A Narrative Review."        
    ## 
    ## $meshHeadings
    ## [1] "NA"
    ## 
    ## $text
    ##  [1] "Nausea and vomiting are common within the palliative care"  
    ##  [2] "population. Antiemetic agents may help control symptoms,"   
    ##  [3] "but may also place patients at risk for QTc prolongation."  
    ##  [4] "This article reviews pharmacotherapy agents including"      
    ##  [5] "anticholinergics, antihistamines, antidopaminergics, 5-HT3" 
    ##  [6] "receptor antagonists, dronabinol, and medical marijuana and"
    ##  [7] "their associated risk of QTc prolongation. A clinical"      
    ##  [8] "treatment pathway is provided to help guide clinicians in"  
    ##  [9] "choosing the most appropriate antiemetic based upon patient"
    ## [10] "specific factors for QTc prolongation."

### KWIC search

The `pmtk_locate_search()` function allows for quick keyword-in-context
(KWIC) search. A simple wrapper of the `corpus::text_locate` function.

``` r
toks <-  corpus::text_tokens(sen_df$abstract)

egs <- PubmedMTK::pmtk_locate_search(text = toks,
                                     doc_id = sen_df$pmid,
                                     search = c('medical marijuana laws'),
                                     stem = F,
                                     window = 15)
knitr::kable(egs[1:5, ])
```

| doc_id   | lhs                                                                                                               | instance               | rhs                                                                                                   |
|:---------|:------------------------------------------------------------------------------------------------------------------|:-----------------------|:------------------------------------------------------------------------------------------------------|
| 34128629 | impaired executive function , cognition , and driving . physicians can recommend use of marijuana under           | medical marijuana laws | but cannot prescribe it , as it is classified as a schedule i controlled substance .                  |
| 33750275 | background : states are rapidly moving to reverse marijuana prohibition , most frequently through legalization of | medical marijuana laws | ( mmls ) , and there is concern that marijuana legalization may affect adolescent marijuana use       |
| 33730400 | more common in states with recreational marijuana laws ( rml ) , followed by states with                          | medical marijuana laws | ( mml ) and without legal cannabis use , respectively . rml and mml were associated                   |
| 33624387 | insurance . we used a differences-in-differences ( dd ) approach and found that the implementation of             | medical marijuana laws | ( mmls ) and recreational marijuana laws ( rmls ) reduced morphine milligram equivalents per enrollee |
| 33143941 | use were assessed each year . cannabis legalization was determined by the presence or absence of                  | medical marijuana laws | ( mml ) and recreational marijuana laws ( rml ) in each state . difference-in-difference approach     |

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

    ## INFO  [13:46:39.947] early stopping at 120 iteration 
    ## INFO  [13:46:40.235] early stopping at 20 iteration

``` r
tm_summary <- PubmedMTK::pmtk_summarize_lda(
  lda = lda, topic_feats_n = 10)
```

#### Feature composition of first ten topics

| topic_id | topic_features                                                                                                                                                                                                |
|---------:|:--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|        1 | cannabis \| cannabinoids \| germany \| quality_of_life \| palliative_care \| canada \| illicit_drugs \| societies,\_medical \| cannabinoid \| glaucoma                                                        |
|        2 | pain \| vomiting \| analgesics \| marijuana_smoking \| muscle_spasticity \| drug_interactions \| cancer_pain \| israel \| clinical_trials_as_topic \| marijuana_abuse                                         |
|        3 | young_adult \| male \| adult \| female \| adolescent \| middle_aged \| age_factors \| prevalence \| socioeconomic_factors \| sex_factors                                                                      |
|        4 | cannabis \| legislation,\_drug \| cannabinoids \| risk_assessment \| drug_prescriptions \| practice_patterns,\_physicians’ \| marijuana \| neoplasms \| research \| treatment_outcome                         |
|        5 | male \| female \| adult \| treatment_outcome \| attitude_of_health_personnel \| prospective_studies \| anxiety \| self_report \| chronic_disease \| adolescent                                                |
|        6 | legislation,\_drug \| marijuana_use \| marijuana_abuse \| marijuana_smoking \| united_states \| prevalence \| male \| epidemiology \| retrospective_studies \| time_factors                                   |
|        7 | cannabis \| marijuana_smoking \| legislation,\_drug \| health_policy \| commerce \| pregnancy \| marijuana_use \| cannabinoids \| female \| plant_extracts                                                    |
|        8 | united_states \| drug_and_narcotic_control \| nausea \| drug_approval \| evidence-based_medicine \| politics \| united_states_food_and_drug_administration \| minnesota \| palliative_care \| risk_assessment |
|        9 | phytotherapy \| cannabis \| united_states \| plant_preparations \| hiv_infections \| physicians \| cancer \| practice_patterns,\_physicians’ \| physician’s_role \| pharmacists                               |
|       10 | dronabinol \| cannabidiol \| cannabinoids \| animals \| multiple_sclerosis \| nabiximols \| tetrahydrocannabinol \| receptors,\_cannabinoid \| cognition \| drug_combinations                                 |

### Two-dimensional analyses

``` r
tmat <- tidytext::cast_sparse(data = tm_summary$topic_word_dist,
                              row = topic_id,
                              column = feature,
                              value = beta)

## build 2d data structures --
two_ds <- PubmedMTK::pmtk_2d(mat = tmat, seed = 99)
```

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

## Interactive HTML topic summary

``` r
PubmedMTK::pmtk_build_interactive(pmtk_lda = tm_summary,
                                  pmtk_2d = two_ds,
                                  out_dir = '/home/jtimm/Desktop/',
                                  file_name = 'party.html')
```

<img src="demo.png" width="100%" />

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
