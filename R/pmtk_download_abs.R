#' Download abstract and meta data for research articles included in PubMed.
#'
#' 
#' @name pmtk_download_abs
#' @param pmids A vector of PMIDs 
#' @param out_file File path for output
#' @param file_prefix String specifying batch name
#' @return A set of data frames as .RDS files  
#' 
#' @importFrom rentrez entrez_fetch
#' @import XML
#' 

#' @export
#' @rdname pmtk_download_abs
#' 
pmtk_download_abs <- function (pmids,
                               out_file,
                               file_prefix) {
  
  batches <- split(pmids, ceiling(seq_along(pmids)/199))  ## -- 
  
  for (i in 1:length(batches)) {
    
    ### tryCatch() -- here -- 
    fetch.pubmed <- rentrez::entrez_fetch(db = "pubmed", 
                                          id = batches[[i]],
                                          rettype = "xml", 
                                          parsed = T)
    
    yy <- as(fetch.pubmed, "character")
    
    xx <- strip_xml(yy)
    encoding <- 'UTF-8'
    Encoding(xx$articletitle) <- encoding
    Encoding(xx$text) <- encoding
    Encoding(xx$authors) <- encoding
    
    ## Also need to make NAs real.  Maybe other function.
    
    xx$text <- ifelse(xx$text %in% 
                        c(' ', 'NA', 'n/a', 'n/a.') | is.na(xx$text),
                      NA, xx$text)
    
    # if(!keep_empty_abs) {
    #   xx <- subset(xx, !is.na(text))}
    
    finally = print(paste0(i, ' / ', length(batches)))
    setwd(out_file)
    saveRDS(xx, file = paste0(file_prefix, i, '.rds'))
  }
}


#' 
#' 
#' 
strip_xml <- function(theFile) {
  newData <- XML::xmlParse(theFile)
  records <- XML::getNodeSet(newData, "//PubmedArticle")
  pmid <- XML::xpathSApply(newData,"//MedlineCitation/PMID", XML::xmlValue)
  doi <- lapply(records, XML::xpathSApply, ".//ELocationID[@EIdType = \"doi\"]", XML::xmlValue)
  doi[sapply(doi, is.list)] <- NA
  doi <- unlist(doi)
  authLast <- lapply(records, XML::xpathSApply, ".//Author/LastName", XML::xmlValue)
  authLast[sapply(authLast, is.list)] <- NA
  authInit <- lapply(records, XML::xpathSApply, ".//Author/Initials", XML::xmlValue)
  authInit[sapply(authInit, is.list)] <- NA
  authors <- mapply(paste, authLast, authInit, collapse = "|")
  
  year <- lapply(records, XML::xpathSApply, ".//PubDate/Year", XML::xmlValue) 
  year[sapply(year, is.list)] <- NA
  year[which(sapply(year, is.na) == TRUE)] <- lapply(records[which(sapply(year, is.na) == TRUE)], 
                                                     XML::xpathSApply, ".//PubDate/MedlineDate", XML::xmlValue)
  year <- gsub(" .+", "", year)
  year <- gsub("-.+", "", year)
  
  
  month <- lapply(records, XML::xpathSApply, ".//PubDate/Month", XML::xmlValue) 
  month[sapply(month, is.list)] <- NA
  month[which(sapply(month, is.na) == TRUE)] <- lapply(records[which(sapply(month, is.na) == TRUE)], 
                                                       XML::xpathSApply, ".//PubDate/MedlineDate", XML::xmlValue)
  
  month <- gsub(" .+", "", month)
  month <- gsub("-.+", "", month)
  
  
  day <- lapply(records, XML::xpathSApply, ".//PubDate/Day", XML::xmlValue) 
  day[sapply(day, is.list)] <- NA
  day[which(sapply(day, is.na) == TRUE)] <- lapply(records[which(sapply(day, is.na) == TRUE)], 
                                                   XML::xpathSApply, ".//PubDate/MedlineDate", XML::xmlValue)
  day <- gsub(" .+", "", day)
  day <- gsub("-.+", "", day)
  
  
  articletitle <- lapply(records, XML::xpathSApply, ".//ArticleTitle", XML::xmlValue) 
  articletitle[sapply(articletitle, is.list)] <- NA
  articletitle <- unlist(articletitle)
  journal <- lapply(records, XML::xpathSApply, ".//ISOAbbreviation", XML::xmlValue) 
  journal[sapply(journal, is.list)] <- NA
  journal <- unlist(journal)
  volume <- lapply(records, XML::xpathSApply, ".//JournalIssue/Volume", XML::xmlValue)
  volume[sapply(volume, is.list)] <- NA
  volume <- unlist(volume)
  issue <- lapply(records, XML::xpathSApply, ".//JournalIssue/Issue", XML::xmlValue)
  issue[sapply(issue, is.list)] <- NA
  issue <- unlist(issue)
  pages <- lapply(records, XML::xpathSApply, ".//MedlinePgn", XML::xmlValue)
  pages[sapply(pages, is.list)] <- NA
  pages <- unlist(pages)
  abstract <- lapply(records, XML::xpathSApply, ".//Abstract/AbstractText", XML::xmlValue)
  abstract[sapply(abstract, is.list)] <- NA
  abstract <- sapply(abstract, paste, collapse = " | ")
  meshHeadings <- lapply(records, XML::xpathSApply, ".//DescriptorName", XML::xmlValue)
  meshHeadings[sapply(meshHeadings, is.list)] <- NA
  meshHeadings <- sapply(meshHeadings, paste, collapse = "|")
  chemNames <- lapply(records, XML::xpathSApply, ".//NameOfSubstance", XML::xmlValue)
  chemNames[sapply(chemNames, is.list)] <- NA
  chemNames <- sapply(chemNames, paste, collapse = "|")
  nctID <- lapply(records, XML::xpathSApply, 
                  ".//DataBank[DataBankName = 'ClinicalTrials.gov']/AccessionNumberList/AccessionNumber", 
                  XML::xmlValue)
  nctID[sapply(nctID, is.null)] <- NA
  nctID <- sapply(nctID, paste, collapse = "|")
  ptype <- lapply(records, XML::xpathSApply, ".//PublicationType", XML::xmlValue)
  ptype[sapply(ptype, is.list)] <- NA
  ptype <- sapply(ptype, paste, collapse = "|")
  
  keywords <- lapply(records, XML::xpathSApply, ".//Keyword", XML::xmlValue)
  keywords[sapply(keywords, is.list)] <- NA
  keywords <- sapply(keywords, paste, collapse = "|")
  
  
  ##just revision -- why -- ?  Only one that is always fully completed --
  pmyear <- lapply(records, XML::xpathSApply,'.//DateRevised/Year', XML::xmlValue)
  pmyear <- gsub(" .+", "", pmyear)
  pmyear <- gsub("-.+", "", pmyear)
  
  pmmonth <- lapply(records, XML::xpathSApply, './/DateRevised/Month', XML::xmlValue) 
  pmmonth <- gsub(" .+", "", pmmonth)
  pmmonth <- gsub("-.+", "", pmmonth)
  
  pmday <- lapply(records, XML::xpathSApply, './/DateRevised/Day', XML::xmlValue) 
  pmday <- gsub(" .+", "", pmday)
  pmday <- gsub("-.+", "", pmday)
  
  ########################

  revision_date <- as.Date(paste(pmyear, pmmonth, pmday, sep = '-'))
    
  data.frame(pmid, doi, authors, year, articletitle, journal, volume, issue, 
             pages, text = abstract, meshHeadings, chemNames, nctID, ptype, keywords,
             revision_date,
             #month, day, pmday, pmmonth, pmyear,
             stringsAsFactors = FALSE) #pmmonth, pmday,
}

