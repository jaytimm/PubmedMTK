PubMed Mining Toolkit: an overview
----------------------------------

-   [PubMed Mining Toolkit: an
    overview](#pubmed-mining-toolkit:-an-overview)
-   [Installation](#installation)
-   [Usage](#usage)
    -   [MeSH vocabulary](#mesh-vocabulary)
    -   [Search the PubMed database -
        `pmtk_search_pubmed`](#search-the-pubmed-database---%60pmtk_search_pubmed%60)
        -   [Sample output:](#sample-output:)
        -   [Summary of record counts returned by PubMed
            query:](#summary-of-record-counts-returned-by-pubmed-query:)
    -   [Advanced counting](#advanced-counting)
    -   [Fetch abstract data from
        PubMed](#fetch-abstract-data-from-pubmed)
        -   [Download batch data:
            `pmtk_download_abs()`](#download-batch-data:-%60pmtk_download_abs()%60)
        -   [Load batch data:
            `pmtk_loadr_abs()`](#load-batch-data:-%60pmtk_loadr_abs()%60)
        -   [Corpus details](#corpus-details)
        -   [Record details](#record-details)
    -   [PubMed search results trend
        data](#pubmed-search-results-trend-data)
    -   [MeSH classifications](#mesh-classifications)
        -   [Example vector
            representation](#example-vector-representation)
    -   [MeSH-based topic model](#mesh-based-topic-model)
    -   [Topic model summary: html
        widget](#topic-model-summary:-html-widget)

Installation
------------

``` r
devtools::install_github("jaytimm/PubmedMTK")
```

Usage
-----

``` r
working_dir <- '/home/jtimm/jt_work/GitHub/PubmedMTK/data-raw/'

if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, # quanteda, 
               rentrez, 
               XML, xml2, RCurl,
               reshape2, #text2vec,  
               tokenizers, 
               tm,
               tidytext,
               Matrix.utils,
               janitor,
               ggplot2, knitr,
               magrittr, dplyr, tidyr)
```

``` r
# Set NCBI API key
# ncbi_key <- '4f47f85a9cc03c4031b3dc274c2840b06108'
# rentrez::set_entrez_key(ncbi_key)
```

### MeSH vocabulary

The package includes as a data frame the MeSH thesaurus/
hierarchically-organized vocabulary – comprised of 2021 versions of
`descriptor` & `trees` files made available via NLM-NIH. [A workflow for
re-creating the table from raw data
sets](https://github.com/jaytimm/PubmedMTK/blob/main/build-MeSH-df.md).

``` r
knitr::kable(head(PubmedMTK::pmtk_tbl_mesh))
```

<table>
<colgroup>
<col style="width: 7%" />
<col style="width: 9%" />
<col style="width: 11%" />
<col style="width: 3%" />
<col style="width: 12%" />
<col style="width: 14%" />
<col style="width: 21%" />
<col style="width: 12%" />
<col style="width: 3%" />
<col style="width: 4%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">DescriptorUI</th>
<th style="text-align: left;">DescriptorName</th>
<th style="text-align: left;">TermName</th>
<th style="text-align: left;">code</th>
<th style="text-align: left;">cats</th>
<th style="text-align: left;">mesh1</th>
<th style="text-align: left;">mesh2</th>
<th style="text-align: left;">tree_location</th>
<th style="text-align: left;">tree1</th>
<th style="text-align: left;">tree2</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">D000001</td>
<td style="text-align: left;">calcimycin</td>
<td style="text-align: left;">calcimycin</td>
<td style="text-align: left;">D</td>
<td style="text-align: left;">Chemicals and Drugs</td>
<td style="text-align: left;">Heterocyclic Compounds</td>
<td style="text-align: left;">Heterocyclic Compounds, Fused-Ring</td>
<td style="text-align: left;">D03.633.100.221.173</td>
<td style="text-align: left;">D03</td>
<td style="text-align: left;">D03.633</td>
</tr>
<tr class="even">
<td style="text-align: left;">D000001</td>
<td style="text-align: left;">calcimycin</td>
<td style="text-align: left;">a-23187</td>
<td style="text-align: left;">D</td>
<td style="text-align: left;">Chemicals and Drugs</td>
<td style="text-align: left;">Heterocyclic Compounds</td>
<td style="text-align: left;">Heterocyclic Compounds, Fused-Ring</td>
<td style="text-align: left;">D03.633.100.221.173</td>
<td style="text-align: left;">D03</td>
<td style="text-align: left;">D03.633</td>
</tr>
<tr class="odd">
<td style="text-align: left;">D000001</td>
<td style="text-align: left;">calcimycin</td>
<td style="text-align: left;">a 23187</td>
<td style="text-align: left;">D</td>
<td style="text-align: left;">Chemicals and Drugs</td>
<td style="text-align: left;">Heterocyclic Compounds</td>
<td style="text-align: left;">Heterocyclic Compounds, Fused-Ring</td>
<td style="text-align: left;">D03.633.100.221.173</td>
<td style="text-align: left;">D03</td>
<td style="text-align: left;">D03.633</td>
</tr>
<tr class="even">
<td style="text-align: left;">D000001</td>
<td style="text-align: left;">calcimycin</td>
<td style="text-align: left;">a23187</td>
<td style="text-align: left;">D</td>
<td style="text-align: left;">Chemicals and Drugs</td>
<td style="text-align: left;">Heterocyclic Compounds</td>
<td style="text-align: left;">Heterocyclic Compounds, Fused-Ring</td>
<td style="text-align: left;">D03.633.100.221.173</td>
<td style="text-align: left;">D03</td>
<td style="text-align: left;">D03.633</td>
</tr>
<tr class="odd">
<td style="text-align: left;">D000001</td>
<td style="text-align: left;">calcimycin</td>
<td style="text-align: left;">antibiotic a23187</td>
<td style="text-align: left;">D</td>
<td style="text-align: left;">Chemicals and Drugs</td>
<td style="text-align: left;">Heterocyclic Compounds</td>
<td style="text-align: left;">Heterocyclic Compounds, Fused-Ring</td>
<td style="text-align: left;">D03.633.100.221.173</td>
<td style="text-align: left;">D03</td>
<td style="text-align: left;">D03.633</td>
</tr>
<tr class="even">
<td style="text-align: left;">D000001</td>
<td style="text-align: left;">calcimycin</td>
<td style="text-align: left;">a23187, antibiotic</td>
<td style="text-align: left;">D</td>
<td style="text-align: left;">Chemicals and Drugs</td>
<td style="text-align: left;">Heterocyclic Compounds</td>
<td style="text-align: left;">Heterocyclic Compounds, Fused-Ring</td>
<td style="text-align: left;">D03.633.100.221.173</td>
<td style="text-align: left;">D03</td>
<td style="text-align: left;">D03.633</td>
</tr>
</tbody>
</table>

### Search the PubMed database - `pmtk_search_pubmed`

Find records included in PubMed that match some search term or multiple
search terms. If multiple search terms are specified, independent
queries are performed per term. Output, then, includes PMID results per
search term – which can subsequently be used to fetch full
records/abstracts.

Search terms are by default translated into NCBI syntax; for simplicity,
search is focused on *MeSH headings* (\[MH\]) and *titles & abstracts*
(\[TIAB\]). So: a search for `aging` is translated as
`aging[MH] OR aging[TIAB]`.

``` r
pmed_search <- c('senescence', 
                 'aging', 
                 #'cancer',
                 'beta galactosidase', 
                 'cell cycle', 
                 'p16',
                 'dna damage', 
                 'cellular senescence', 
                 'induced senescence',
                 'secretory phenotype')
```

``` r
search_results1 <- PubmedMTK::pmtk_search_pubmed(pmed_search = pmed_search)
```

#### Sample output:

``` r
search_results1 %>%
  head() %>%
  knitr::kable()
```

| search     | pmid     |
|:-----------|:---------|
| senescence | 33631788 |
| senescence | 33631598 |
| senescence | 33631183 |
| senescence | 33630953 |
| senescence | 33629548 |
| senescence | 33629281 |

#### Summary of record counts returned by PubMed query:

``` r
# ## Total citations per search term are summarized below:
search_results1 %>%
  group_by(search) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  janitor::adorn_totals() %>%
  knitr::kable()
```

| search              |        n|
|:--------------------|--------:|
| cell cycle          |   416905|
| aging               |   368927|
| senescence          |   278659|
| dna damage          |   134703|
| beta galactosidase  |    33254|
| induced senescence  |    25843|
| cellular senescence |    25707|
| p16                 |    14690|
| secretory phenotype |     2703|
| Total               |  1301391|

### Advanced counting

Quick inspection of query results – before fetching record details.

``` r
query_bigrams <- PubmedMTK::pmtk_query_bigrams(search_results1) 
## crosstab_qresults()
```

### Fetch abstract data from PubMed

As a two-step process: using functions `pmtk_download_abs()` and
`pmtk_loadr_abs(()`.

`rentrez` is a lovely package (and maintained by
[rOpenSci](https://github.com/ropensci/rentrez)); however, in my
experience it is not especially well-designed for fetching PubMed
abstracts in bulk or building text corpora. API rate-limits being most
problematic.

While `rentrez` is still employed here, our approach utilizes a
combination of local storage and “more + smaller” API queries to make
the most of rate limits. Each “batch” contains n = 199 records; batch
files are converted from XML to a data frame in RDS format and stored
locally in a user-specified file path.

#### Download batch data: `pmtk_download_abs()`

The `out_file` parameter specifies the file path for local batch file
storage; the `file_prefix` parameter specifies a character string used
to identify batches (along with a batch \#).

``` r
PubmedMTK::pmtk_download_abs(pmids = unique(search_results1$pmid),
                             out_file = paste0(working_dir, 'batches/'),
                             file_prefix = 'sen')
```

#### Load batch data: `pmtk_loadr_abs()`

The `pmtk_loadr_abs()` function loads batch files as two data frames:
the first, a corpus object containing the record id and abstract, and
the second, a metadata object including record id and all other record
details, eg, article name, MeSH terms, Pub Date, etc.

``` r
batch_dir <- paste0(working_dir, 'batches/')
sen_df <- PubmedMTK::pmtk_loadr_abs(in_file = batch_dir, 
                                    file_prefix = 'sen')
```

#### Corpus details

``` r
sen_df$tif %>%
  mutate(includes_abstract = ifelse(is.na(text), 'N', 'Y')) %>%
  count(includes_abstract) %>%
  mutate(tokens = ifelse(includes_abstract == 'Y',
                         n* 210, NA)) %>%
  knitr::kable()
```

| includes\_abstract |      n|   tokens|
|:-------------------|------:|--------:|
| N                  |   1428|       NA|
| Y                  |  24837|  5215770|

#### Record details

``` r
sen_df$meta %>%
  filter(complete.cases(.)) %>% 
  ## !! NA's are not proper stil -- !!!
  slice(1) %>%
  data.table::transpose(keep.names = "var") %>%
  mutate(V1 = gsub('\\|', ' \\| ', V1)) %>%
  knitr::kable()
```

<table>
<colgroup>
<col style="width: 5%" />
<col style="width: 94%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">var</th>
<th style="text-align: left;">V1</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">pmid</td>
<td style="text-align: left;">30053915</td>
</tr>
<tr class="even">
<td style="text-align: left;">doi</td>
<td style="text-align: left;">10.1186/s13073-018-0568-8</td>
</tr>
<tr class="odd">
<td style="text-align: left;">authors</td>
<td style="text-align: left;">Scepanovic P | Alanio C | Hammer C | Hodel F | Bergstedt J | Patin E | Thorball CW | Chaturvedi N | Charbit B | Abel L | Quintana-Murci L | Duffy D | Albert ML | Fellay J</td>
</tr>
<tr class="even">
<td style="text-align: left;">year</td>
<td style="text-align: left;">2018</td>
</tr>
<tr class="odd">
<td style="text-align: left;">articletitle</td>
<td style="text-align: left;">Human genetic variants and age are the strongest predictors of humoral immune responses to common pathogens and vaccines.</td>
</tr>
<tr class="even">
<td style="text-align: left;">journal</td>
<td style="text-align: left;">Genome Med</td>
</tr>
<tr class="odd">
<td style="text-align: left;">volume</td>
<td style="text-align: left;">10</td>
</tr>
<tr class="even">
<td style="text-align: left;">issue</td>
<td style="text-align: left;">1</td>
</tr>
<tr class="odd">
<td style="text-align: left;">pages</td>
<td style="text-align: left;">59</td>
</tr>
<tr class="even">
<td style="text-align: left;">meshHeadings</td>
<td style="text-align: left;">Adult | Aged | Aging | Bacterial Infections | Female | HLA-D Antigens | Humans | Immunity, Humoral | Male | Middle Aged | Polymorphism, Single Nucleotide | Vaccines | Virus Diseases</td>
</tr>
<tr class="odd">
<td style="text-align: left;">chemNames</td>
<td style="text-align: left;">HLA-D Antigens | Vaccines</td>
</tr>
<tr class="even">
<td style="text-align: left;">nctID</td>
<td style="text-align: left;">NCT01699893</td>
</tr>
<tr class="odd">
<td style="text-align: left;">ptype</td>
<td style="text-align: left;">Clinical Trial | Journal Article | Research Support, Non-U.S. Gov’t</td>
</tr>
<tr class="even">
<td style="text-align: left;">keywords</td>
<td style="text-align: left;">Age | GWAS | HLA | Human genomics | Humoral immunity | Immunoglobulins | Infection | Serology | Sex | Vaccination</td>
</tr>
<tr class="odd">
<td style="text-align: left;">revision_date</td>
<td style="text-align: left;">17883</td>
</tr>
</tbody>
</table>

### PubMed search results trend data

``` r
## in theory, this could be used for other things -- 
tr <- subset(sen_df$meta, !grepl('[a-z]', year))
tr$year <- as.Date(paste(tr$year, '01', '01', sep = '-'))
## 
tr1 <- tr[search_results1, on = 'pmid']

## 
meds <- data.table::data.table(PubmedMTK::pmtk_tbl_citations)
tr2 <-  tr1[, list(n = .N), by = list(search, year)]
tr3 <- subset(tr2, year > as.Date('1969', format = '%Y') &
                year < as.Date('2019', format = '%Y') )

tr4 <- meds[tr3, on = 'year']
tr4$per_100k = round(tr4$n / tr4$total * 100000, 3)
```

``` r
## Via ggplot --
tr4 %>%
  ggplot() +
  geom_line(aes(x = year,
                #y = n, 
                y = per_100k,
                group = search,
                color = search),
            size = 1
            ) +
  theme_minimal() +
  ggthemes::scale_color_stata() +
  theme(legend.position = 'right',
        legend.title = element_blank())  +
  ylab('Per 100,000 Medline citations') +
  ggtitle('wrong')
```

![](README_files/figure-markdown_github/unnamed-chunk-20-1.png)

### MeSH classifications

Extract KEYWORDS, MeSH HEADINGS & CHEM-NAMES – output is a
MeSH-comprised text representation.

``` r
## this takes too long --based
meshes <- PubmedMTK::pmtk_gather_mesh(meta_df = sen_df$meta)
txts <- length(unique(meshes$pmid))


## get frequencies -- 
freqs <-  meshes[, list(doc_freq = length(unique(pmid))), 
                 by = list(descriptor_name)]
freqs$doc_prop <- freqs$doc_freq/ txts
freqs1 <- subset(freqs, doc_prop > 0.0001 & doc_prop < 0.02)

meshes1 <- subset(meshes, descriptor_name %in% freqs1$descriptor_name)
meshes1 <- subset(meshes1, nchar(descriptor_name) > 0)
```

#### Example vector representation

``` r
set.seed(999)
meshes1 %>%
  filter(pmid %in% sample(unique(meshes1$pmid), 10)) %>%
  group_by(pmid) %>%
  summarize (d = paste0(descriptor_name, collapse = ' | ')) %>%
  knitr::kable()
```

<table>
<colgroup>
<col style="width: 1%" />
<col style="width: 98%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">pmid</th>
<th style="text-align: left;">d</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">1017451</td>
<td style="text-align: left;">achievement | aptitude | intelligence tests | pattern recognition, visual | psychometrics</td>
</tr>
<tr class="even">
<td style="text-align: left;">1439082</td>
<td style="text-align: left;">chloroplasts | hordeum | plant proteins | plant proteins</td>
</tr>
<tr class="odd">
<td style="text-align: left;">19675177</td>
<td style="text-align: left;">accidents, traffic | automobile driving | depressive disorder | geriatric assessment | health status | physical fitness | probability | risk assessment | task performance and analysis</td>
</tr>
<tr class="even">
<td style="text-align: left;">19816867</td>
<td style="text-align: left;">aorta | cardiac catheterization | diaphragm | dyspnea | echocardiography, transesophageal | heart septal defects, atrial | hemodynamics | hypoxia | magnetic resonance imaging | paralysis | severity of illness index | treatment outcome</td>
</tr>
<tr class="odd">
<td style="text-align: left;">23461035</td>
<td style="text-align: left;">butyrates | cell proliferation | fibroblasts | histone deacetylase inhibitors | map kinase signaling system | mechanistic target of rapamycin complex 1 | multiprotein complexes | phosphorylation | sodium | tor serine-threonine kinases | p38 mitogen-activated protein kinases | butyrates | histone deacetylase inhibitors | multiprotein complexes | sodium | tor serine-threonine kinases | mechanistic target of rapamycin complex 1 | p38 mitogen-activated protein kinases</td>
</tr>
<tr class="even">
<td style="text-align: left;">25555041</td>
<td style="text-align: left;">cognitive function | overactive bladder | treatment | activities of daily living | benzhydryl compounds | benzofurans | cognition disorders | depression | follow-up studies | geriatric assessment | muscarinic antagonists | pyrrolidines | quality of life | treatment outcome | urinary bladder, overactive | benzhydryl compounds | benzofurans | muscarinic antagonists | pyrrolidines</td>
</tr>
<tr class="odd">
<td style="text-align: left;">25653617</td>
<td style="text-align: left;">alzheimer’s | dti | mci</td>
</tr>
<tr class="even">
<td style="text-align: left;">32240105</td>
<td style="text-align: left;">bioinformatics | tumor microenvironment</td>
</tr>
<tr class="odd">
<td style="text-align: left;">32392179</td>
<td style="text-align: left;">creatine | hemolysis | lifespan</td>
</tr>
<tr class="even">
<td style="text-align: left;">8420683</td>
<td style="text-align: left;">ovarian neoplasms | prognosis | randomized controlled trials as topic | retrospective studies</td>
</tr>
</tbody>
</table>

### MeSH-based topic model

> Latent Dirichlet allocation: a topic modeling algorithm that models
> **each document** in corpus as a composite of topics, and **each
> topic** as a composite of terms.

Here, we utilize the MeSH-based abstract representations to build a
topic model. *Exploratory utility*: (1) no real text processing, (2)
information-dense, ie, no fluff, just MeSH. Topic composition can be
interpreted as sets of MeSH terms that frequently co-occur.

``` r
mesh_dtm <- tidytext::cast_sparse(data = meshes1,
                                  row = pmid,
                                  column = descriptor_name,
                                  value = count)

mesh_lda <- text2vec::LDA$new(n_topics = 30) ## This is the model
topic_model_fit <- mesh_lda$fit_transform(mesh_dtm, progressbar = F)
```

The `mtk_summarize_lda` function summarizes and extracts topic
composition from the `text2vec::LDA` output. For each possible
topic-feature pair, the model computes the likelihood a given topic
generated a given feature. Output is filtered to the highest scoring
features per topic using the `topic_feats_n`.

``` r
tm_summary <- PubmedMTK::mtk_summarize_lda(
  lda = mesh_lda, topic_feats_n = 15)
```

<table>
<colgroup>
<col style="width: 2%" />
<col style="width: 97%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: right;">topic_id</th>
<th style="text-align: left;">topic_features</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: right;">1</td>
<td style="text-align: left;">retrospective studies | treatment outcome | prognosis | follow-up studies | risk assessment | severity of illness index | hiv infections | incidence | pain | postoperative complications | predictive value of tests | clinical trials as topic | survival analysis | stroke | comorbidity</td>
</tr>
<tr class="even">
<td style="text-align: right;">2</td>
<td style="text-align: left;">calcium | hippocampus | osteoporosis | bone density | bone and bones | neuronal plasticity | synapses | rats, sprague-dawley | vitamin d | radiography | absorptiometry, photon | lumbar vertebrae | synaptic transmission | action potentials | electrophysiology</td>
</tr>
<tr class="odd">
<td style="text-align: right;">3</td>
<td style="text-align: left;">mice, transgenic | amyloid beta-peptides | phosphorylation | autophagy | peptide fragments | membrane proteins | tau proteins | alzheimer’s disease | parkinson disease | neurodegenerative diseases | amyloid beta-protein precursor | astrocytes | microglia | nerve tissue proteins | alzheimer’s disease</td>
</tr>
<tr class="even">
<td style="text-align: right;">4</td>
<td style="text-align: left;">phenotype | mutation | genotype | genetic predisposition to disease | gene expression profiling | polymorphism, single nucleotide | dna methylation | epigenesis, genetic | case-control studies | gene expression regulation | alleles | genetic variation | apolipoproteins e | polymorphism, genetic | transcriptome</td>
</tr>
<tr class="odd">
<td style="text-align: right;">5</td>
<td style="text-align: left;">cell differentiation | fibroblasts | cell line | cell division | cell survival | stem cells | cell proliferation | gene expression regulation | mice, knockout | gene expression | regeneration | epithelial cells | reverse transcriptase polymerase chain reaction | extracellular matrix | cell movement</td>
</tr>
<tr class="even">
<td style="text-align: right;">6</td>
<td style="text-align: left;">cytokines | t-lymphocytes | spleen | mice, inbred balb c | lymphocyte activation | macrophages | flow cytometry | tumor necrosis factor-alpha | thymus gland | lymphocytes | interleukin-6 | b-lymphocytes | antibodies, monoclonal | hematopoietic stem cells | immunoglobulin g</td>
</tr>
<tr class="odd">
<td style="text-align: right;">7</td>
<td style="text-align: left;">kinetics | erythrocytes | liver | erythrocyte aging | rabbits | proteins | dna | cell membrane | iron | electrophoresis, polyacrylamide gel | in vitro techniques | hemoglobins | molecular weight | isoenzymes | biological transport</td>
</tr>
<tr class="even">
<td style="text-align: right;">8</td>
<td style="text-align: left;">mortality | life expectancy | population dynamics | demography | research | socioeconomic factors | europe | population | health | geriatrics | forecasting | demographic factors | statistics as topic | retirement | developing countries</td>
</tr>
<tr class="odd">
<td style="text-align: right;">9</td>
<td style="text-align: left;">pregnancy | infant, newborn | microscopy, electron | immunohistochemistry | fetus | cell count | axons | spinal cord | gestational age | fertility | oocytes | cats | embryonic and fetal development | reproduction | histocytochemistry</td>
</tr>
<tr class="even">
<td style="text-align: right;">10</td>
<td style="text-align: left;">cell line, tumor | cell proliferation | tumor suppressor protein p53 | cell cycle | gene expression regulation, neoplastic | antineoplastic agents | breast neoplasms | mesenchymal stem cells | dna damage | cell transformation, neoplastic | down-regulation | neoplasms | up-regulation | micrornas | rna, small interfering</td>
</tr>
<tr class="odd">
<td style="text-align: right;">11</td>
<td style="text-align: left;">geriatrics | health services for the aged | nursing homes | long-term care | hospitalization | health services needs and demand | homes for the aged | health knowledge, attitudes, practice | attitude of health personnel | home care services | geriatric nursing | delivery of health care | health promotion | medicare | health policy</td>
</tr>
<tr class="even">
<td style="text-align: right;">12</td>
<td style="text-align: left;">rats, sprague-dawley | rats, wistar | rats, inbred f344 | immunohistochemistry | in vitro techniques | endothelium, vascular | dopamine | brain chemistry | rats, inbred strains | cerebral cortex | nitric oxide | norepinephrine | dose-response relationship, drug | cerebellum | corpus striatum</td>
</tr>
<tr class="odd">
<td style="text-align: right;">13</td>
<td style="text-align: left;">depression | social support | quality of life | adaptation, psychological | stress, psychological | caregivers | interpersonal relations | mental health | self concept | anxiety | family | attitude to health | personal satisfaction | social environment | health status</td>
</tr>
<tr class="even">
<td style="text-align: right;">14</td>
<td style="text-align: left;">skin aging | skin | collagen | ultraviolet rays | face | rejuvenation | retina | cosmetic techniques | patient satisfaction | macular degeneration | treatment outcome | visual acuity | light | skin neoplasms | rhytidoplasty</td>
</tr>
<tr class="odd">
<td style="text-align: right;">15</td>
<td style="text-align: left;">prevalence | health status | activities of daily living | logistic models | chronic disease | socioeconomic factors | multivariate analysis | european continental ancestry group | geriatric assessment | comorbidity | health surveys | age distribution | odds ratio | incidence | disabled persons</td>
</tr>
<tr class="even">
<td style="text-align: right;">16</td>
<td style="text-align: left;">testosterone | organ size | estradiol | testis | rats, inbred strains | estrogens | insulin-like growth factor i | growth hormone | hypothalamus | sex characteristics | luteinizing hormone | androgens | ovary | sexual maturation | hydrocortisone</td>
</tr>
<tr class="odd">
<td style="text-align: right;">17</td>
<td style="text-align: left;">obesity | insulin | blood glucose | adipose tissue | glucose | diabetes mellitus | diabetes mellitus, type 2 | cholesterol | lipids | insulin resistance | body mass index | body composition | energy metabolism | lipid metabolism | triglycerides</td>
</tr>
<tr class="even">
<td style="text-align: right;">18</td>
<td style="text-align: left;">kidney | dose-response relationship, drug | liver | rats, inbred strains | tissue distribution | kinetics | kidney diseases | drug interactions | creatinine | sodium | potassium | glomerular filtration rate | rats, wistar | intestinal mucosa | kidney failure, chronic</td>
</tr>
<tr class="odd">
<td style="text-align: right;">19</td>
<td style="text-align: left;">inflammation | neoplasms | telomere | models, biological | dogs | lung | senescence | telomerase | models, animal | lens, crystalline | osteoarthritis | cancer | cataract | cartilage, articular | leukocytes</td>
</tr>
<tr class="even">
<td style="text-align: right;">20</td>
<td style="text-align: left;">magnetic resonance imaging | cognition disorders | cerebral cortex | image processing, computer-assisted | disease progression | brain mapping | atrophy | hippocampus | cognitive dysfunction | diagnosis, differential | cerebrovascular circulation | frontal lobe | positron-emission tomography | neural pathways | case-control studies</td>
</tr>
<tr class="odd">
<td style="text-align: right;">21</td>
<td style="text-align: left;">exercise | geriatric assessment | activities of daily living | elderly | frail elderly | older adults | accidental falls | quality of life | cognitive dysfunction | walking | independent living | physical fitness | frailty | gait | postural balance</td>
</tr>
<tr class="even">
<td style="text-align: right;">22</td>
<td style="text-align: left;">reproducibility of results | sensitivity and specificity | models, biological | algorithms | reference values | water | computer simulation | stress, mechanical | analysis of variance | models, statistical | materials testing | temperature | surface properties | microscopy, electron, scanning | elasticity</td>
</tr>
<tr class="odd">
<td style="text-align: right;">23</td>
<td style="text-align: left;">molecular sequence data | base sequence | dna-binding proteins | amino acid sequence | transcription factors | transcription, genetic | dna repair | dna | plant leaves | protein binding | saccharomyces cerevisiae | promoter regions, genetic | gene expression regulation, plant | dna damage | nuclear proteins</td>
</tr>
<tr class="even">
<td style="text-align: right;">24</td>
<td style="text-align: left;">mitochondria | antioxidants | reactive oxygen species | oxidation-reduction | superoxide dismutase | hydrogen peroxide | lipid peroxidation | plant extracts | glutathione | free radicals | dna, mitochondrial | catalase | ascorbic acid | oxygen | vitamin e</td>
</tr>
<tr class="odd">
<td style="text-align: right;">25</td>
<td style="text-align: left;">diet | cattle | species specificity | chickens | energy intake | reproduction | seasons | swine | dietary supplements | feeding behavior | dietary proteins | eating | temperature | phospholipids | amino acids</td>
</tr>
<tr class="even">
<td style="text-align: right;">26</td>
<td style="text-align: left;">muscle, skeletal | reference values | muscles | sarcopenia | biomechanical phenomena | sleep | sex characteristics | muscle contraction | muscle strength | exercise | body composition | adaptation, physiological | posture | movement | electromyography</td>
</tr>
<tr class="odd">
<td style="text-align: right;">27</td>
<td style="text-align: left;">behavior, animal | motor activity | stress, physiological | gene expression regulation | caenorhabditis elegans | circadian rhythm | drosophila melanogaster | models, biological | maze learning | caenorhabditis elegans proteins | homeostasis | environment | caloric restriction | mutation | drosophila</td>
</tr>
<tr class="even">
<td style="text-align: right;">28</td>
<td style="text-align: left;">memory | reaction time | analysis of variance | attention | psychomotor performance | mental recall | memory disorders | cognition disorders | learning | electroencephalography | memory, short-term | photic stimulation | executive function | visual perception | pattern recognition, visual</td>
</tr>
<tr class="odd">
<td style="text-align: right;">29</td>
<td style="text-align: left;">cardiovascular diseases | smoking | regression analysis | body mass index | sex characteristics | japan | menopause | follow-up studies | life style | nutritional status | risk | linear models | postmenopause | alcohol drinking | prevalence</td>
</tr>
<tr class="even">
<td style="text-align: right;">30</td>
<td style="text-align: left;">blood pressure | hypertension | myocardium | heart rate | heart | heart failure | hemodynamics | oxygen consumption | electrocardiography | myocardial infarction | aorta | arteriosclerosis | arteries | cardiovascular diseases | antihypertensive agents</td>
</tr>
</tbody>
</table>

### Topic model summary: html widget

``` r
## topic model html widget
#mesh_lda$plot()
mesh_lda$plot(out.dir = "ldavis", open.browser = FALSE)
```

![](README_files/figure-markdown_github/demo-tm-viz.png)
