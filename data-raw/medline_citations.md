## Medline citation counts by year of publication

For calculation of relative citation frequency. Count data made
available
[here](https://www.nlm.nih.gov/bsd/medline_cit_counts_yr_pub.html).
Medline counts are complete through 2018; here, we have projected counts
by year for 2019 onward using the R package `forecast`. An imperfect
solution, but a reasonable enough denominator for computing rough
relative frequency counts historically.

``` r
if (!require("pacman")) install.packages("pacman")
pacman::p_load(magrittr, dplyr, ggplot2, rvest, xml2, forecast)
```

### Build table

``` r
url <- 'https://www.nlm.nih.gov/bsd/medline_cit_counts_yr_pub.html'

medl <- url %>%
  xml2::read_html() %>%
  rvest::html_node(xpath = '//*[@id="main-body"]/div/table[1]') %>% 
  rvest::html_table(fill = TRUE) %>%
  select(1:3)

colnames(medl)<- c('year', 'total', 'usa')
medl1 <- medl %>% filter(!grepl('[a-z]', year)) %>%
  mutate(year = as.numeric(gsub('[[:punct:]]', '', year)),
         total = as.numeric(gsub('[[:punct:]]', '', total)),
         usa = as.numeric(gsub('[[:punct:]]', '', usa))) %>%
  filter(year < 2019) %>%
  arrange(year)

## simple forecast for 2019 & 2020
total <- forecast::forecast(ts(medl1$total), h = 3)
usa <- forecast::forecast(ts(medl1$usa), h = 3)

forecasted <- data.frame(year = c(2019:2021),
                         total = round(total$mean[1:3]),
                         usa = round(usa$mean[1:3]))
pmtk_tbl_citations <- rbind(medl1, forecasted) 
pmtk_tbl_citations$year <- as.Date(paste(pmtk_tbl_citations$year , 
                                         1, 1, sep = "-"))
```

``` r
setwd('/home/jtimm/jt_work/GitHub/PubmedMTK/data')
usethis::use_data(pmtk_tbl_citations, overwrite=TRUE)
```

### Hitorical citation counts

``` r
pmtk_tbl_citations %>%
  filter(year < '2021-01-01') %>%
  ggplot(aes(x = year, 
             y = total)) +
  geom_col(show.legend = FALSE, 
           alpha = 0.75,
           fill = 'steelblue') + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_date(date_breaks = '10 years', date_labels = "%Y") +
  labs(title = 'Medline citation counts by year, 1947 - 2020')
```

![](medline_citations_files/figure-markdown_github/unnamed-chunk-4-1.png)
