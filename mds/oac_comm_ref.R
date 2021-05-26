

desk_dir <- '/home/jtimm/Desktop/oac/'
git_dir <- '/home/jtimm/jt_work/GitHub/PubmedMTK/data-raw'

setwd(git_dir)
oac <- read.csv('oa_comm_use_file_list.csv') 
oac <- janitor::clean_names(oac)
oac1 <- subset(oac, !is.na(pmid))


PubmedMTK::pmtk_download_abs(pmids = oac1$pmid,
                             out_file = desk_dir,
                             file_prefix = 'ref')

sen_df <- PubmedMTK::pmtk_loadr_abs(in_file = desk_dir, 
                                    file_prefix = 'ref')



meshes <- PubmedMTK::pmtk_gather_mesh(meta_df = sen_df$meta)
txts <- length(unique(meshes$pmid))


meshey <- meshes
meshey1 <- subset(meshey, nchar(descriptor_name) > 0)
meshey1$type <- factor(meshey1$type, 
                       levels = c('mesh_heading', 'chem_name', 'keyword'))
meshey1$descriptor_name <- tolower(meshey1$descriptor_name)

data.table::setorder(meshey1, pmid, type)

meshey2 <- meshey1[, .SD[1], list(pmid, descriptor_name)]


## get frequencies -- 
freqs <-  meshey2[, list(doc_freq = length(unique(pmid))), 
                 by = list(type, descriptor_name)]
freqs$doc_prop <- freqs$doc_freq/ txts
pmtk_tbl_pmc_ref <- subset(freqs, doc_prop > 0.00001)



# meshes1 <- subset(meshes, 
#                   descriptor_name %in% freqs2$descriptor_name)
# 
# meshes2 <- subset(meshes1, nchar(descriptor_name) > 0)
# 
# 
# 
# mesh_dtm <- tidytext::cast_sparse(data = meshes2,
#                                   row = pmid,
#                                   column = descriptor_name,
#                                   value = count)


setwd('/home/jtimm/jt_work/GitHub/PubmedMTK/data/')
usethis::use_data(pmtk_tbl_pmc_ref, overwrite=TRUE)


