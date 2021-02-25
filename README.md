## PubMed Mining Toolkit: an overview

## Installation

``` r
devtools::install_github("jaytimm/PubmedMTK")
```

## Usage

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
`descriptor` & `trees` files made available via NLM-NIH. An account of
how the table was constructed is detailed
[here](https://github.com/jaytimm/PubmedMTK/blob/main/build-MeSH-df.md).

Also, [a useful
reference](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5324252/).

``` r
knitr::kable(head(PubmedMTK::pmtk_tbl_mesh))
```

    ## Warning: replacing previous import 'data.table::melt' by 'reshape2::melt' when
    ## loading 'PubmedMTK'

| DescriptorUI | DescriptorName | TermName           | code | cats                | mesh1                  | mesh2                              | tree_location       | tree1 | tree2   |
|:-------------|:---------------|:-------------------|:-----|:--------------------|:-----------------------|:-----------------------------------|:--------------------|:------|:--------|
| D000001      | calcimycin     | calcimycin         | D    | Chemicals and Drugs | Heterocyclic Compounds | Heterocyclic Compounds, Fused-Ring | D03.633.100.221.173 | D03   | D03.633 |
| D000001      | calcimycin     | a-23187            | D    | Chemicals and Drugs | Heterocyclic Compounds | Heterocyclic Compounds, Fused-Ring | D03.633.100.221.173 | D03   | D03.633 |
| D000001      | calcimycin     | a 23187            | D    | Chemicals and Drugs | Heterocyclic Compounds | Heterocyclic Compounds, Fused-Ring | D03.633.100.221.173 | D03   | D03.633 |
| D000001      | calcimycin     | a23187             | D    | Chemicals and Drugs | Heterocyclic Compounds | Heterocyclic Compounds, Fused-Ring | D03.633.100.221.173 | D03   | D03.633 |
| D000001      | calcimycin     | antibiotic a23187  | D    | Chemicals and Drugs | Heterocyclic Compounds | Heterocyclic Compounds, Fused-Ring | D03.633.100.221.173 | D03   | D03.633 |
| D000001      | calcimycin     | a23187, antibiotic | D    | Chemicals and Drugs | Heterocyclic Compounds | Heterocyclic Compounds, Fused-Ring | D03.633.100.221.173 | D03   | D03.633 |

### Search the PubMed database: `pmtk_search_pubmed()`

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
## tester:: rentrez_search <- rentrez::entrez_search(term = 'cancer', db = 'pubmed')

pmed_search <- c('senescence', 
                 'aging', 
                 'cancer',
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

Sample output:

``` r
search_results1 %>%
  head() %>%
  knitr::kable()
```

| search     | pmid     |
|:-----------|:---------|
| senescence | 33618333 |
| senescence | 33618126 |
| senescence | 33617078 |
| senescence | 33616799 |
| senescence | 33616187 |
| senescence | 33616006 |

Summary of record counts returned by PubMed query:

``` r
# ## Total citations per search term are summarized below:
search_results1 %>%
  group_by(search) %>%
  summarise(n = n()) %>%
  arrange(desc(n)) %>%
  janitor::adorn_totals() %>%
  knitr::kable()
```

| search              |       n |
|:--------------------|--------:|
| cancer              | 3929619 |
| cell cycle          |  416790 |
| aging               |  368768 |
| senescence          |  278509 |
| dna damage          |  134624 |
| beta galactosidase  |   33247 |
| induced senescence  |   25826 |
| cellular senescence |   25673 |
| p16                 |   14681 |
| secretory phenotype |    2705 |
| Total               | 5230442 |

### Advanced counting

Quick inspection of query results – before fetching record details.

``` r
query_bigrams <- PubmedMTK::pmtk_query_bigrams(search_results1) ## crosstab_qresults()
```

### Fetch abstract data from PubMed

As a two-step process: using functions `pmtk_download_abs()` and
`pmtk_loadr_abs(()`.

While `rentrez` is a lovely package (and maintained by
[ROpenSci](https://github.com/ropensci/rentrez)), in my experience it is
not especially well-designed for fetching PubMed abstracts in bulk or
building text corpora. API rate-limits being most problematic.

The approach taken here utilizes a combination of local storage and
“more + smaller” API queries to make the most of rate limits – employing
the `entrez_fetch` function from the `rentrez` package to perform
individual queries. Each “batch” contains n = 199 records; batch files
are converted from XML to a data frame in RDS format and stored locally
in a user-specified file path.

#### Download batch data: `pmtk_download_abs()`

The `out_file` parameter specifies the file path for local batch file
storage is indicated via The `file_prefix` parameter specifies a
character string used to identify batches (along with a batch \#).

``` r
PubmedMTK::pmtk_download_abs(pmids = sen_pmids$pmid,
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
sen_df <- PubmedMTK::pmtk_loadr_abs(in_file = batch_dir, file_prefix = 'sen')
```

#### Record details

``` r
sen_df$meta %>%
  filter(complete.cases(.)) %>% ## !! NA's are not proper stil -- !!!
  slice(1) %>%
  data.table::transpose(keep.names = "var") %>%
  mutate(V1 = gsub('\\|', ' \\| ', V1)) %>%
  knitr::kable()
```

| var           | V1                                                                                                                                                                            |
|:--------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| pmid          | 33608630                                                                                                                                                                      |
| doi           | 10.1038/s42003-020-01619-4                                                                                                                                                    |
| authors       | Kim G \| Kim M \| Kim M \| Park C \| Yoon Y \| Lim DH \| Yeo H \| Kang S \| Lee YG \| Beak NI \| Lee J \| Kim S \| Kwon JY \| Choi WW \| Lee C \| Yoon KW \| Park H \| Lee DG |
| year          | 2021                                                                                                                                                                          |
| articletitle  | Spermidine-induced recovery of human dermal structure and barrier function by skin microbiome.                                                                                |
| journal       | Commun Biol                                                                                                                                                                   |
| volume        | 4                                                                                                                                                                             |
| issue         | 1                                                                                                                                                                             |
| pages         | 231                                                                                                                                                                           |
| meshHeadings  | NA                                                                                                                                                                            |
| chemNames     | NA                                                                                                                                                                            |
| nctID         | NA                                                                                                                                                                            |
| ptype         | Journal Article                                                                                                                                                               |
| keywords      | NA                                                                                                                                                                            |
| revision_date | 18680                                                                                                                                                                         |

### Trend data:

``` r
## pmtk_

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

![](README_files/figure-markdown_github/unnamed-chunk-17-1.png)

### MeSH classifications

Extract KEYWORDS, MeSH HEADINGS & CHEM-NAMES – output is a
MeSH-comprised vector representation –

``` r
## this takes too long -- 
meshes <- PubmedMTK::pmtk_gather_mesh(meta_df = sen_df$meta)
txts <- length(unique(meshes$pmid))


## get frequencies -- 
freqs <-  meshes[, list(doc_freq = length(unique(pmid))), 
                 by = list(descriptor_name)]
freqs$doc_prop <- freqs$doc_freq/ txts
freqs1 <- subset(freqs, doc_prop > 0.0001 & doc_prop < 0.02)

meshes1 <- subset(meshes, descriptor_name %in% freqs1$descriptor_name)
```

Example MeSH-based vector representation:

### MeSH-based topic model

> Latent Dirichlet allocation: a topic modeling algorithm that models
> **each document** in corpus as a composite of topics, and **each
> topic** as a composite of terms.

Here, we utilize the MeSH-based abstract representations to build a
topic model. *Exploratory utility*: (1) no real text processing, (2)
information-dense, ie, no fluff, just MeSH.

``` r
mesh_dtm <- tidytext::cast_sparse(data = meshes1,
                                  row = pmid,
                                  column = descriptor_name,
                                  value = count)

mesh_lda <- text2vec::LDA$new(n_topics = 30) ## This is the model
topic_model_fit <- mesh_lda$fit_transform(mesh_dtm, progressbar = F)
```

    ## INFO  [22:44:47.530] early stopping at 180 iteration 
    ## INFO  [22:45:05.634] early stopping at 20 iteration

The `mtk_summarize_lda` function summarizes and extracts topic
composition from the `text2vec::LDA` output.

``` r
topic_model_summary <- PubmedMTK::mtk_summarize_lda(lda = mesh_lda, 
                                                    topic_feats_n = 15)

summary <- topic_model_summary[ , .(topic_features = paste0(variable, collapse = ' | ')), 
                                by = topic_id]
knitr::kable(summary)
```

| topic_id | topic_features                                                                                                                                                                                                                                                                                                                                            |
|---------:|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|        1 | cell division \| microscopy, electron \| liver \| telomere \| dna \| proteins \| rats, inbred strains \| cell survival \| fibroblasts \| tissue distribution \| cell nucleus \| kinetics \| kidney \| cell count \| mice, inbred strains                                                                                                                  |
|        2 | reference values \| sex characteristics \| bone and bones \| stroke \| tomography, x-ray computed \| osteoarthritis \| diagnosis, differential \| regression analysis \| case-control studies \| osteoporosis \| ultrasonography \| reproducibility of results \| cartilage, articular \| femur \| sensitivity and specificity                            |
|        3 | mortality \| life expectancy \| population dynamics \| demography \| socioeconomic factors \| forecasting \| europe \| population \| health \| retirement \| research \| demographic factors \| fertility \| developing countries \| employment                                                                                                           |
|        4 | molecular sequence data \| base sequence \| dna damage \| amino acid sequence \| dna-binding proteins \| transcription, genetic \| dna repair \| mutation \| transcription factors \| dna \| nuclear proteins \| promoter regions, genetic \| polymerase chain reaction \| protein binding \| histones                                                    |
|        5 | health services for the aged \| nursing homes \| geriatrics \| long-term care \| caregivers \| qualitative research \| homes for the aged \| health knowledge, attitudes, practice \| health services needs and demand \| decision making \| family \| attitude of health personnel \| home care services \| health promotion \| geriatric nursing        |
|        6 | blood pressure \| hypertension \| cardiovascular diseases \| heart rate \| diabetes mellitus \| kidney \| diabetes mellitus, type 2 \| arteriosclerosis \| hemodynamics \| aorta \| coronary disease \| electrocardiography \| atherosclerosis \| arteries \| kidney diseases                                                                             |
|        7 | prevalence \| geriatrics \| japan \| european continental ancestry group \| regression analysis \| china \| socioeconomic factors \| african americans \| asian continental ancestry group \| models, statistical \| educational status \| health surveys \| sex \| research \| history, 20th century                                                     |
|        8 | retrospective studies \| follow-up studies \| prognosis \| incidence \| risk assessment \| treatment outcome \| prevalence \| logistic models \| comorbidity \| predictive value of tests \| multivariate analysis \| hospitalization \| proportional hazards models \| severity of illness index \| odds ratio                                           |
|        9 | cell differentiation \| immunohistochemistry \| stem cells \| parkinson disease \| mice, knockout \| cell proliferation \| mice, transgenic \| gene expression regulation, developmental \| nerve tissue proteins \| mesenchymal stem cells \| regeneration \| astrocytes \| gene expression \| drosophila \| spinal cord                                 |
|       10 | diet \| calcium \| chickens \| dietary supplements \| nutritional status \| energy intake \| feeding behavior \| amino acids \| eating \| dietary proteins \| dietary fats \| nutritional physiological phenomena \| intestinal mucosa \| cattle \| swine                                                                                                 |
|       11 | hippocampus \| rats, sprague-dawley \| cerebral cortex \| rats, wistar \| synapses \| in vitro techniques \| dopamine \| rats, inbred f344 \| neuronal plasticity \| cerebellum \| electric stimulation \| norepinephrine \| axons \| immunohistochemistry \| rats, inbred strains                                                                        |
|       12 | phosphorylation \| blotting, western \| gene expression regulation \| caenorhabditis elegans \| cell line \| mice, knockout \| reverse transcriptase polymerase chain reaction \| membrane proteins \| enzyme activation \| transcription factors \| up-regulation \| carrier proteins \| caenorhabditis elegans proteins \| autophagy \| gene expression |
|       13 | models, biological \| inflammation \| rats, sprague-dawley \| dose-response relationship, drug \| caloric restriction \| neoplasms \| homeostasis \| models, animal \| autophagy \| nitric oxide \| rats, wistar \| neurodegenerative diseases \| senescence \| endothelial cells \| ageing                                                               |
|       14 | postural balance \| reference values \| lung \| biomechanical phenomena \| gait \| retina \| posture \| walking \| movement \| accidental falls \| visual acuity \| pain \| electromyography \| respiration \| oxygen                                                                                                                                     |
|       15 | muscle, skeletal \| obesity \| exercise \| body mass index \| body composition \| bone density \| osteoporosis \| sarcopenia \| adipose tissue \| muscle strength \| physical fitness \| anthropometry \| absorptiometry, photon \| muscle contraction \| oxygen consumption                                                                              |
|       16 | memory \| reaction time \| analysis of variance \| attention \| psychomotor performance \| mental recall \| behavior, animal \| learning \| memory disorders \| maze learning \| photic stimulation \| motor activity \| memory, short-term \| visual perception \| electroencephalography                                                                |
|       17 | antioxidants \| reactive oxygen species \| oxidation-reduction \| superoxide dismutase \| mitochondria \| hydrogen peroxide \| lipid peroxidation \| plant extracts \| rats, wistar \| lens, crystalline \| free radicals \| glutathione \| cataract \| chromatography, high pressure liquid \| catalase                                                  |
|       18 | cytokines \| t-lymphocytes \| spleen \| mice, inbred balb c \| inflammation \| lymphocyte activation \| macrophages \| tumor necrosis factor-alpha \| thymus gland \| lymphocytes \| flow cytometry \| interleukin-6 \| b-lymphocytes \| hematopoietic stem cells \| antibodies, monoclonal                                                               |
|       19 | analysis of variance \| hydrogen-ion concentration \| in vitro techniques \| temperature \| kinetics \| water \| stress, mechanical \| materials testing \| surface properties \| microscopy, electron, scanning \| ethanol \| hot temperature \| solubility \| drug interactions \| staining and labeling                                                |
|       20 | depression \| smoking \| cognition disorders \| chronic disease \| life style \| mental disorders \| psychiatric status rating scales \| exercise \| anxiety \| stress, psychological \| depressive disorder \| alcohol drinking \| psychometrics \| physical activity \| motor activity                                                                  |
|       21 | models, biological \| stress, physiological \| species specificity \| reproduction \| circadian rhythm \| plant leaves \| seasons \| gene expression regulation, plant \| light \| adaptation, physiological \| dogs \| melatonin \| plant proteins \| arabidopsis \| larva                                                                               |
|       22 | cell proliferation \| neoplasms \| cell line, tumor \| tumor suppressor protein p53 \| cell cycle \| gene expression regulation, neoplastic \| breast neoplasms \| antineoplastic agents \| cell transformation, neoplastic \| cancer \| lung neoplasms \| fibroblasts \| tumor cells, cultured \| senescence \| cell survival                            |
|       23 | liver \| mitochondria \| insulin \| blood glucose \| glucose \| cholesterol \| energy metabolism \| muscles \| lipids \| lipid metabolism \| insulin resistance \| triglycerides \| adenosine triphosphate \| fatty acids \| muscle development                                                                                                           |
|       24 | skin aging \| skin \| treatment outcome \| collagen \| ultraviolet rays \| face \| rejuvenation \| fibroblasts \| cosmetic techniques \| extracellular matrix \| wound healing \| patient satisfaction \| double-blind method \| rhytidoplasty \| hyaluronic acid                                                                                         |
|       25 | amyloid beta-peptides \| myocardium \| mice, transgenic \| heart \| peptide fragments \| hiv infections \| tau proteins \| alzheimer’s disease \| alzheimer’s disease \| disease progression \| heart failure \| amyloid beta-protein precursor \| calcium \| organ size \| amyloid                                                                       |
|       26 | magnetic resonance imaging \| cognitive dysfunction \| cognition disorders \| reproducibility of results \| image processing, computer-assisted \| brain mapping \| cerebral cortex \| atrophy \| sleep \| algorithms \| sensitivity and specificity \| executive function \| hippocampus \| disease progression \| mild cognitive impairment             |
|       27 | pregnancy \| infant, newborn \| erythrocytes \| erythrocyte aging \| fetus \| rabbits \| cattle \| sheep \| iron \| gestational age \| kinetics \| oocytes \| hemoglobins \| swine \| blood proteins                                                                                                                                                      |
|       28 | quality of life \| activities of daily living \| geriatric assessment \| health status \| social support \| frail elderly \| adaptation, psychological \| independent living \| frailty \| disabled persons \| mental health \| interpersonal relations \| older adults \| depression \| self concept                                                     |
|       29 | phenotype \| mutation \| genotype \| genetic predisposition to disease \| gene expression profiling \| case-control studies \| polymorphism, single nucleotide \| dna methylation \| epigenesis, genetic \| alleles \| genetic variation \| biological evolution \| polymorphism, genetic \| transcriptome \| apolipoproteins e                           |
|       30 | testosterone \| menopause \| estradiol \| organ size \| testis \| estrogens \| insulin-like growth factor i \| ovary \| growth hormone \| luteinizing hormone \| hypothalamus \| androgens \| hydrocortisone \| sexual maturation \| sex characteristics                                                                                                  |

### Topic model summary: html widget

``` r
## topic model html widget
```
