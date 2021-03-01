PubMed Mining Toolkit - an overview
-----------------------------------

-   [Installation](#installation)
-   [Usage](#usage)
    -   [MeSH vocabulary](#mesh-vocabulary)
    -   [Search the PubMed database -
        `pmtk_search_pubmed()`](#search-the-pubmed-database---%60pmtk_search_pubmed()%60)
    -   [Advanced counting -
        `pmtk_crosstab_query()`](#advanced-counting---%60pmtk_crosstab_query()%60)
    -   [Fetch abstract data from
        PubMed](#fetch-abstract-data-from-pubmed)
    -   [PubMed search results trend
        data](#pubmed-search-results-trend-data)
    -   [Extract MeSH classifications -
        `pmtk_gather_mesh()`](#extract-mesh-classifications---%60pmtk_gather_mesh()%60)
    -   [MeSH-based topic model](#mesh-based-topic-model)
    -   [Topic model summary - html
        widget](#topic-model-summary---html-widget)
    -   [Image summary](#image-summary)

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

### Search the PubMed database - `pmtk_search_pubmed()`

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
pmed_search <- c('human life span',
                 'senescence', 
                 'proteostasis',
                 'dna damage', 
                 'beta galactosidase', 
                 'genomic instability',
                 'telomere attrition',
                 'epigenetic alterations',
                 'mitochondrial dysfunction',
                 'cellular senescence',
                 'stem cell exhaustion')
```

``` r
search_results1 <- PubmedMTK::pmtk_search_pubmed(pmed_search = pmed_search)
```

> Summary of record counts returned by PubMed query

``` r
# ## Total citations per search term are summarized below:
search_results1 %>%
  group_by(search) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  #janitor::adorn_totals() %>%
  knitr::kable()
```

| search                    |       n|
|:--------------------------|-------:|
| senescence                |  278659|
| dna damage                |  134703|
| beta galactosidase        |   33254|
| cellular senescence       |   25707|
| mitochondrial dysfunction |   22063|
| genomic instability       |   21671|
| human life span           |   12145|
| epigenetic alterations    |    5049|
| proteostasis              |    3625|
| telomere attrition        |     853|
| stem cell exhaustion      |     162|

> Unique records

``` r
length(unique(search_results1$pmid))
```

    ## [1] 489299

### Advanced counting - `pmtk_crosstab_query()`

Quick analysis of search term co-occurrence based on results from
`pmtk_search_pubmed()`. Here, *term-A* and *term-B* are said to co-occur
in *abstract-X* if independent PubMed queries for *term-A* and *term-B*
both return *abstract-X*. Ideal for quick exploration.

``` r
search_tab <- PubmedMTK::pmtk_crosstab_query(search_results1) %>%
  mutate(pmi = round(log( (n1n2/1e6) / ( (n1/1e6) * (n2/1e6) )), 3))

search_tab %>% filter(term1 == 'senescence') %>% knitr::kable()
```

| term1      | term2                     |      n1|      n2|   n1n2|     pmi|
|:-----------|:--------------------------|-------:|-------:|------:|-------:|
| senescence | beta galactosidase        |  278659|   33254|   3214|  -1.059|
| senescence | cellular senescence       |  278659|   25707|  16005|   0.804|
| senescence | dna damage                |  278659|  134703|   5685|  -1.887|
| senescence | epigenetic alterations    |  278659|    5049|    201|  -1.946|
| senescence | genomic instability       |  278659|   21671|   1054|  -1.746|
| senescence | human life span           |  278659|   12145|  12003|   1.266|
| senescence | mitochondrial dysfunction |  278659|   22063|   1450|  -1.445|
| senescence | proteostasis              |  278659|    3625|    470|  -0.765|
| senescence | stem cell exhaustion      |  278659|     162|     88|   0.668|
| senescence | telomere attrition        |  278659|     853|    443|   0.623|

``` r
search_tab %>%
  mutate(pmi = ifelse(pmi < 0, 0, pmi),
         term2 = as.character(term2)) %>%
  rowwise() %>%
  
  mutate(ords = paste(sort(c(term1, term2)), collapse = '_')) %>%
  group_by(ords) %>% slice(1) %>% ungroup() %>%
  
  #arrange(term1, term2) %>%
  ggplot(aes(x = term1, y = term2, fill = pmi)) + 
  geom_tile() + 
  geom_text(aes(fill = pmi, label = n1n2), size = 3) + 
  scale_fill_gradient2(low = scales::muted("#d8b365"), 
                       mid = "#f5f5f5", 
                       high = scales::muted('#5ab4ac'),
                       midpoint = 0) +
  theme_minimal() + xlab('') + ylab('') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, hjust = 1))
```

![](README_files/figure-markdown_github/unnamed-chunk-14-1.png)

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

#### Download batch data - `pmtk_download_abs()`

The `out_file` parameter specifies the file path for local batch file
storage; the `file_prefix` parameter specifies a character string used
to identify batches (along with a batch \#).

``` r
PubmedMTK::pmtk_download_abs(pmids = unique(search_results1$pmid),
                             out_file = paste0(working_dir, 'batches/'),
                             file_prefix = 'sen')
```

#### Load batch data - `pmtk_loadr_abs()`

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

| includes\_abstract |       n|    tokens|
|:-------------------|-------:|---------:|
| N                  |   62790|        NA|
| Y                  |  426053|  89471130|

#### Record details

``` r
sen_df$meta %>%
  filter(complete.cases(.)) %>% 
  slice(1) %>%
  data.table::transpose(keep.names = "var") %>%
  mutate(V1 = gsub('\\|', ' \\| ', V1)) %>%
  knitr::kable()
```

<table>
<colgroup>
<col style="width: 3%" />
<col style="width: 96%" />
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
<td style="text-align: left;">31748285</td>
</tr>
<tr class="even">
<td style="text-align: left;">doi</td>
<td style="text-align: left;">10.1136/bmjopen-2018-027984</td>
</tr>
<tr class="odd">
<td style="text-align: left;">authors</td>
<td style="text-align: left;">Ali D | Callan N | Ennis S | Powell R | McGuire S | McGregor G | Weickert MO | Miller MA | Cappuccio FP | Banerjee P</td>
</tr>
<tr class="even">
<td style="text-align: left;">year</td>
<td style="text-align: left;">2019</td>
</tr>
<tr class="odd">
<td style="text-align: left;">articletitle</td>
<td style="text-align: left;">Heart failure with preserved ejection fraction (HFpEF) pathophysiology study (IDENTIFY-HF): does increased arterial stiffness associate with HFpEF, in addition to ageing and vascular effects of comorbidities? Rationale and design.</td>
</tr>
<tr class="even">
<td style="text-align: left;">journal</td>
<td style="text-align: left;">BMJ Open</td>
</tr>
<tr class="odd">
<td style="text-align: left;">volume</td>
<td style="text-align: left;">9</td>
</tr>
<tr class="even">
<td style="text-align: left;">issue</td>
<td style="text-align: left;">11</td>
</tr>
<tr class="odd">
<td style="text-align: left;">pages</td>
<td style="text-align: left;">e027984</td>
</tr>
<tr class="even">
<td style="text-align: left;">meshHeadings</td>
<td style="text-align: left;">Aged | Aged, 80 and over | Aging | Biomarkers | Comorbidity | Diabetes Mellitus | Echocardiography | Exercise Tolerance | Female | Heart Failure | Heart Ventricles | Humans | Hypertension | Male | Observational Studies as Topic | Prospective Studies | Pulse Wave Analysis | Research Design | Stroke Volume | Vascular Stiffness</td>
</tr>
<tr class="odd">
<td style="text-align: left;">chemNames</td>
<td style="text-align: left;">Biomarkers</td>
</tr>
<tr class="even">
<td style="text-align: left;">nctID</td>
<td style="text-align: left;">NCT03186833</td>
</tr>
<tr class="odd">
<td style="text-align: left;">ptype</td>
<td style="text-align: left;">Journal Article</td>
</tr>
<tr class="even">
<td style="text-align: left;">keywords</td>
<td style="text-align: left;">arterial stiffness | comorbidities | heart failure with preserved ejection fraction | pathophysiology</td>
</tr>
<tr class="odd">
<td style="text-align: left;">revision_date</td>
<td style="text-align: left;">18565</td>
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

*FACET on absolute versus relative !!!*

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

![](README_files/figure-markdown_github/unnamed-chunk-21-1.png)

### Extract MeSH classifications - `pmtk_gather_mesh()`

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

> Example vector representations:

``` r
set.seed(999)
meshes1 %>%
  filter(pmid %in% sample(unique(meshes1$pmid), 5)) %>%
  group_by(pmid) %>%
  summarize (mesh_reps = paste0(descriptor_name, collapse = ' | ')) %>%
  knitr::kable()
```

<table>
<colgroup>
<col style="width: 3%" />
<col style="width: 96%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">pmid</th>
<th style="text-align: left;">mesh_reps</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">1439082</td>
<td style="text-align: left;">chloroplasts | hordeum | plant proteins | plant proteins</td>
</tr>
<tr class="even">
<td style="text-align: left;">19675177</td>
<td style="text-align: left;">accidents, traffic | automobile driving | depressive disorder | geriatric assessment | health status | physical fitness | probability | risk assessment | task performance and analysis</td>
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

#### Feature composition of first ten topics

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
</tbody>
</table>

### Topic model summary - html widget

``` r
## topic model html widget
#mesh_lda$plot()
mesh_lda$plot(out.dir = "ldavis", open.browser = FALSE)
```

![](README_files/figure-markdown_github/demo-tm-viz.png)

### Image summary

![](README_files/figure-markdown_github/summary.png)
