## A simple restructure of MeSH ontology

### Intro

> A clean & un-adulterated version of MeSH ontology as data frame.
> Included in the R package `PumbedMTK`. Based on two files:
> `desc2021.xml` & `mtrees2021.bin`; available via
> [nlm.nih.gov](https://www.nlm.nih.gov/databases/download/mesh.html).

``` r
if (!require("pacman")) install.packages("pacman")
pacman::p_load(magrittr, dplyr, tidyr, xml2)

setwd(git_dir)
desc <- xml2::read_xml('desc2020.xml')
trees <- read.csv('mtrees2021.bin', header = FALSE, sep =';')
```

### `desc` file

> Extract descriptor details (& concepts & terms) from descriptor file.

``` r
## Descriptor
DescriptorUI <- desc %>% 
  xml2::xml_find_all('.//DescriptorUI') %>% 
  xml2::xml_text()

DescriptorName <- desc %>% 
  xml2::xml_find_all('.//DescriptorName') %>%  
  xml2::xml_text() 

descriptor <- data.frame(DescriptorUI, DescriptorName) %>%
  distinct() %>% 
  arrange(DescriptorUI)

## Concepts
ConceptName <- desc %>% 
  xml2::xml_find_all('.//ConceptName') %>%  
  xml2::xml_text()

ConceptUI <- desc %>% 
  xml2::xml_find_all('.//ConceptUI') %>%  
  xml2::xml_text()

concept <- data.frame(ConceptUI, ConceptName)

## Terms
TermUI <- desc %>% 
  xml2::xml_find_all('.//TermUI') %>%  
  xml2::xml_text() 

TermName <- desc %>% xml2::xml_find_all('.//Term') %>%  
  xml2::xml_find_all('String') %>%
  xml2::xml_text()

term <- data.frame(TermUI, TermName) 
```

> Note: Primary concept == primary term == descriptor name == mesh term.

``` r
concept_term <- term %>%
  
  left_join(concept, by = c('TermName' = 'ConceptName')) %>%
  fill(ConceptUI) %>%
  
  left_join(descriptor, by = c('TermName' = 'DescriptorName')) %>%
  
  mutate(DescriptorName = ifelse(is.na(DescriptorUI), NA, TermName)) %>%
  
  fill(DescriptorUI, DescriptorName) %>%
  select(DescriptorUI, DescriptorName, ConceptUI,
         TermUI, TermName)
```

### `mtrees` file

``` r
tree <- trees %>%
  rename(mesh_heading = V1,
         tree_location = V2) %>%
  select(tree_location, mesh_heading)


# Extract the two highest parent nodes from tree location.  
## For general classification purposes.  

level1 <- tree[nchar(tree$tree_location) == 3, ]
level2 <- tree[nchar(tree$tree_location) == 7, ] %>%
  mutate(join = gsub('\\....', '', tree_location))

top_parents <- level2 %>%
  left_join(level1, by = c('join' = 'tree_location'))

colnames(top_parents) <- c('tree2', 'mesh2', 'tree1', 'mesh1')
top_parents <- top_parents[, c(3:4, 1:2)]
```

> Manully add labels for highest-level node in ontology:

``` r
### 2-4 High-level categories
cats <- 
c('Anatomy', 'Organisms', 'Diseases', 'Chemicals and Drugs',
  'Analytical, Diagnostic and Therapeutic Techniques, and Equipment', 
  'Psychiatry and Psychology', 'Phenomena and Processes', 
  'Disciplines and Occupations', 
  'Anthropology, Education, Sociology, and Social Phenomena', 
  'Technology, Industry, and Agriculture', 'Humanities', 
  'Information Science', 'Named Groups', 'Health Care',
  'Publication Characteristics', 'Geographicals')

code <- c(LETTERS[1:14], 'V', 'Z')

high_tree <- data.frame(code, cats)

meta <- top_parents %>%
  mutate(code = gsub('..$', '', tree1)) %>%
  left_join(high_tree) %>%
  select(code:cats, tree1:mesh2)
```

### Joining two data sets

> Which adds the MeSH tree location to descriptor data via `MeSH`
> heading & `DescriptorName` variables. Note that a single
> `DescriptorName` may be classified (in tree structure) in multiple
> ways.

``` r
clean_col <- function(x) {
  x <- enc2utf8(x)
  x <- trimws(x)
  x <- tolower(x)  }

### 2-5 Join metadata & descriptors/terms
pmtk_tbl_mesh <- concept_term %>%
  
  left_join(tree, by = c('DescriptorName' = 'mesh_heading')) %>%
  mutate (tree2 = substring(tree_location, 1, 7)) %>%
  
  left_join(meta) %>%
  
  # rename(descriptor_id = DescriptorUI,
  #        descriptor_name = DescriptorName,
  #        term_name = TermName) %>%
 
  select(DescriptorUI:DescriptorName, TermName, 
         code, cats,
         mesh1, mesh2,
         tree_location, tree1, tree2) %>%
  
  mutate_at(c('DescriptorName', 'TermName'), clean_col) %>%
  
  filter(complete.cases(.))
```

### Sample records

``` r
knitr::kable(head(PubmedMTK::pmtk_tbl_mesh))
```

| DescriptorUI | DescriptorName | TermName           | code | cats                | mesh1                  | mesh2                              | tree_location       | tree1 | tree2   |
|:-------------|:---------------|:-------------------|:-----|:--------------------|:-----------------------|:-----------------------------------|:--------------------|:------|:--------|
| D000001      | calcimycin     | calcimycin         | D    | Chemicals and Drugs | Heterocyclic Compounds | Heterocyclic Compounds, Fused-Ring | D03.633.100.221.173 | D03   | D03.633 |
| D000001      | calcimycin     | a-23187            | D    | Chemicals and Drugs | Heterocyclic Compounds | Heterocyclic Compounds, Fused-Ring | D03.633.100.221.173 | D03   | D03.633 |
| D000001      | calcimycin     | a 23187            | D    | Chemicals and Drugs | Heterocyclic Compounds | Heterocyclic Compounds, Fused-Ring | D03.633.100.221.173 | D03   | D03.633 |
| D000001      | calcimycin     | a23187             | D    | Chemicals and Drugs | Heterocyclic Compounds | Heterocyclic Compounds, Fused-Ring | D03.633.100.221.173 | D03   | D03.633 |
| D000001      | calcimycin     | antibiotic a23187  | D    | Chemicals and Drugs | Heterocyclic Compounds | Heterocyclic Compounds, Fused-Ring | D03.633.100.221.173 | D03   | D03.633 |
| D000001      | calcimycin     | a23187, antibiotic | D    | Chemicals and Drugs | Heterocyclic Compounds | Heterocyclic Compounds, Fused-Ring | D03.633.100.221.173 | D03   | D03.633 |

### Output

``` r
pmtk_tbl_mesh <- data.table::data.table(pmtk_tbl_mesh)
setwd('/home/jtimm/jt_work/GitHub/PubmedMTK/data')
usethis::use_data(pmtk_tbl_mesh, overwrite=TRUE)
```
