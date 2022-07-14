# R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
# Copyright (C) 2022 The R Foundation for Statistical Computing
# Platform: x86_64-pc-linux-gnu (64-bit)

library("stringr")
library("tidyr")
library("dplyr")
library("data.table")

#Read identity mappings
read_identity_mappings <- function(path) {
  x <- data.table::fread(path,
                         sep = "\t",
                         fill = TRUE,
                         check.names = FALSE,
                         stringsAsFactors = FALSE)
  x$V10 <- gsub("d__|p__|c__|o__|f__|g__|s__", "", x$V10) #strip rank tags
  x$V10 <- gsub("[", "", x$V10, fixed=TRUE) #strip rank tags
  x$V10 <- gsub("]", "", x$V10, fixed=TRUE) #strip rank tags
  x$V10 <- gsub("d_bxid|p_bxid|c_bxid|o_bxid|f_bxid|g_bxid|s_bxid", "bxid", x$V10) #strip rank tags
  x.taxonomy <- cbind(x, stringr::str_split_fixed(x$V10, ";", 8))
  colnames(x.taxonomy) <- c(colnames(x), "asv", "domain", "phylum", "class", "order", "family", "genus", "species")
  x.taxonomy[which(x.taxonomy[,4] =="*"), "V4"] <- "51"
  x.taxonomy$V4 <- as.numeric(x.taxonomy$V4)
  x.taxonomy[which(x.taxonomy[,4] <99), "species"] <- ""
  x.taxonomy[which(x.taxonomy[,4] <95.8), "genus"] <- ""
  x.taxonomy[which(x.taxonomy[,4] <89.6), "family"] <- ""
  x.taxonomy[which(x.taxonomy[,4] <86.2), "order"] <- ""
  x.taxonomy[which(x.taxonomy[,4] <83.5), "class"] <- ""
  x.taxonomy[which(x.taxonomy[,4] <80.8), "phylum"] <- ""
  ###
  x.taxonomy[which(x.taxonomy[,"order"] =="Chloroplast"), "species"] <- "Chloroplast_Chloroplast"
  x.taxonomy[which(x.taxonomy[,"order"] =="Chloroplast"), "genus"] <- "Chloroplast"
  x.taxonomy[which(x.taxonomy[,"order"] =="Chloroplast"), "family"] <- "Chloroplast"
  ###
  x.taxonomy[which(x.taxonomy[,"family"] =="Mitochondria"), "species"] <- "Mitochondria_Mitochondria"
  x.taxonomy[which(x.taxonomy[,"family"] =="Mitochondria"), "genus"] <- "Mitochondria"
  ###
  x.taxonomy[apply(x.taxonomy, 1:2, function(i) grepl('unidentified', i))] <- ""
  x.taxonomy[apply(x.taxonomy, 1:2, function(i) grepl("_sp$", i, ignore.case=FALSE))] <- ""
  x.taxonomy[apply(x.taxonomy, 1:2, function(i) grepl("NA", i, ignore.case=FALSE))] <- ""
  x.taxonomy[apply(x.taxonomy, 1:2, function(i) grepl("unclass_", i, ignore.case=FALSE))] <- ""
  x.taxonomy[apply(x.taxonomy, 1:2, function(i) grepl("Incertae_sedis", i, ignore.case=FALSE))] <- ""
  ###
  x.taxonomy[is.na(x.taxonomy)] <- ""
  x.taxonomy$row_no <- stringr::str_split_fixed(x.taxonomy$V9, "_", 2)[,2]
  x.taxonomy$V10 <- x$V10
  x.taxonomy <- x.taxonomy %>% mutate_at(vars(row_no), as.numeric) %>% arrange(., row_no)
  x <- x.taxonomy
}


# Define function to read and sort data from the denovo clustering in UCLUST-format  ##
read_sort_mappings <- function(path, colnames) {
  x <- data.table::fread(path,
                         sep = "\t",
                         fill = TRUE,
                         check.names = FALSE,
                         stringsAsFactors = FALSE)
  #keep only rows with H (hits) and S (singletons) in V1, keep only V9+V10
  x <- x[V1 %in% c("H", "S"),.(V9, V10)]
  #if * in V10, replace with V9
  x <- x[,V10 := ifelse(V10 == "*", V9, V10)]
  #order by FLASV ID
  x <- x[order(as.integer(gsub("[^0-9+$]|\\..*$", "", V9))),]
  colnames(x) <- colnames
  invisible(x)
}

write_tax <- function(tax, file) {
  data.table::fwrite(tax,
                     file,
                     quote = TRUE,
                     sep = ",")
}


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 1) Load UC tables
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
UCout_BEExact.id <- read_identity_mappings(path="bac/UCout_BEExact_seedsonly_usearch.uc")

UCout_NCBI.id <- read_identity_mappings(path="bac/UCout_NCBI_usearch_R-adjusted.uc")

#Join BEExact and NCBI
UCout.mer <- merge(UCout_BEExact.id, UCout_NCBI.id, by = "row_no", all = T)  %>%
                                                    mutate(id.mer = ifelse(V4.x < V4.y , V4.y, V4.x)) %>%
                                                    mutate(species.mer = ifelse(species.x=="", species.y, species.x)) %>%
                                                    mutate(genus.mer = ifelse(genus.x=="", genus.y, genus.x)) %>%
                                                    mutate(family.mer = ifelse(family.x=="", family.y, family.x)) %>%
                                                    mutate(order.mer = ifelse(order.x=="", order.y, order.x)) %>%
                                                    mutate(class.mer = ifelse(class.x=="", class.y, class.x)) %>%
                                                    mutate(phylum.mer = ifelse(phylum.x=="", phylum.y, phylum.x)) %>%
                                                    select(row_no, id.mer, domain.x, phylum.mer, class.mer, order.mer, family.mer, genus.mer, species.mer)

UCout.mer$genus.split <- stringr::str_split_fixed(UCout.mer$species.mer, "_", 2)[,1]
UCout.mer$species.split <- stringr::str_split_fixed(UCout.mer$species.mer, "_", 2)[,2]
UCout.mer <- UCout.mer %>% mutate(g = ifelse(genus.split=="", genus.mer, genus.split))
UCout.mer$g_s <- paste(UCout.mer$g, UCout.mer$species.split, sep="_")
UCout.mer <- UCout.mer %>% mutate(g_s2 = ifelse(species.mer=="", species.mer, g_s))

UCout.mer <- UCout.mer %>% select(row_no, id.mer, domain.x, phylum.mer, class.mer, order.mer, family.mer, g, g_s2)
colnames(UCout.mer) <- c("row_no", "id_mer", "domain", "phylum", "class", "order", "family", "genus", "species")

#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 2) Load cluster_smallmem tables
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ASV_S<- read_sort_mappings(paste0("bac/clustering/ASV_SPECIES.uc"), c("asv", "species"))
S_G <- read_sort_mappings(paste0("bac/clustering/SPECIES_GENUS.uc"), c("species", "genus"))
G_F <- read_sort_mappings(paste0("bac/clustering/GENUS_FAMILY.uc"), c("genus", "family"))
F_O <- read_sort_mappings(paste0("bac/clustering/FAMILY_ORDER.uc"), c("family", "order"))
O_C <- read_sort_mappings(paste0("bac/clustering/ORDER_CLASS.uc"), c("order", "class"))
C_P <- read_sort_mappings(paste0("bac/clustering/CLASS_PHYLUM.uc"), c("class", "phylum"))

#merge each taxonomic level according to the mapping results
denovo_tax <- NULL
denovo_tax <- left_join(ASV_S, S_G, by ="species")
denovo_tax <- left_join(denovo_tax, G_F, by = "genus")
denovo_tax <- left_join(denovo_tax, F_O, by = "family")
denovo_tax <- left_join(denovo_tax, O_C, by = "order")
denovo_tax <- left_join(denovo_tax, C_P, by = "class")

#reorder columns and remove length from FLASV ID's (".xxxx")
denovo_tax <- denovo_tax[,c("asv", "phylum", "class", "order", "family", "genus", "species")]
denovo_tax[,2:7] <- lapply(denovo_tax[,2:7], gsub, pattern = "\\..*$", replacement = "")
#generate denovo names per taxonomic level based on FLASV ID
denovo_tax[["species"]] <- gsub("^[^0-9]+", paste0("bsv"), denovo_tax[["species"]])
denovo_tax[["genus"]] <- gsub("^[^0-9]+", paste0("bsv"), denovo_tax[["genus"]])
denovo_tax[["family"]] <- gsub("^[^0-9]+", paste0("bsv"), denovo_tax[["family"]])
denovo_tax[["order"]] <- gsub("^[^0-9]+", paste0("bsv"), denovo_tax[["order"]])
denovo_tax[["class"]] <- gsub("^[^0-9]+", paste0("bsv"), denovo_tax[["class"]])
denovo_tax[["phylum"]] <- gsub("^[^0-9]+", paste0("bsv"), denovo_tax[["phylum"]])
denovo_tax$row_no <- stringr::str_split_fixed(denovo_tax$asv, "_", 2)[,2]
denovo_tax <- denovo_tax %>% mutate_at(vars(row_no), as.numeric) %>% arrange(., row_no)

#write out
#write_tax(denovo_tax, file = paste0("denovo.csv"))

#denovo_tax$row_no <- NULL
poly_denovo_tax <- NULL
poly_denovo_tax <- list()
poly_denovo_tax$species <- denovo_tax[, .(asv ==first(asv), nParents = uniqueN(genus), genus = paste0(sort(unique(genus)), collapse= ", "), new = first(genus)), by = species][nParents > 1,]
poly_denovo_tax$genus <- denovo_tax[, .(asv ==first(asv), nParents = uniqueN(family), family = paste0(sort(unique(family)), collapse= ", "), new = first(family)), by = genus][nParents > 1,]
poly_denovo_tax$family <- denovo_tax[, .(asv ==first(asv), nParents = uniqueN(order), order = paste0(sort(unique(order)), collapse= ", "), new = first(order)), by = family][nParents > 1,]
poly_denovo_tax$order <- denovo_tax[, .(asv ==first(asv), nParents = uniqueN(class), class = paste0(sort(unique(class)), collapse= ", "), new = first(class)), by = order][nParents > 1,]
poly_denovo_tax$class <- denovo_tax[, .(asv ==first(asv), nParents = uniqueN(phylum), phylum = paste0(sort(unique(phylum)), collapse= ", "), new = first(phylum)), by = class][nParents > 1,]

polyTaxaLog <- unlist(sapply(poly_denovo_tax, function(x) {
  if(nrow(x) > 1) {
    paste0(colnames(x)[1],
           " ",
           x[[1]],
           " has ",
           x[[3]],
           " parents: \"",
           x[[4]],
           "\", and has been assigned the ",
           colnames(x)[4],
           " of ",
           x[[2]],
           ": ",
           x[[5]])
  }
}), use.names = FALSE)
nTaxa <- length(polyTaxaLog)

#issue a warning if one or more taxa had more than one parent, write out logfile
if(nTaxa > 0) {
  warning(paste0(nTaxa,
                 " taxa had more than one parent, see the logfile \"./output/polyphyletics.log\" for details"),
          call. = FALSE)
  #writeLines(polyTaxaLog,
  #           paste0("bac/clustering/polyphyletics.log"))
}

#fix them
denovo_tax[, genus := first(genus), by = species]
denovo_tax[, family := first(family), by = genus]
denovo_tax[, order := first(order), by = family]
denovo_tax[, class := first(class), by = order]
denovo_tax[, phylum := first(phylum), by = class]

denovo_tax$species <- paste(denovo_tax$genus, denovo_tax$species, sep="_")
#write out
#write_tax(merged_tax, file = paste0(outputfolder, "/tax_complete.csv"))


#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# 3) Merge UC + cluster_smallmem tables
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#Merge
UCout.denovo.mer <- merge(UCout.mer, denovo_tax, by = "row_no", all = T) %>%
                                              mutate(domain.fill = ifelse(domain=="", "Bacteria", domain)) %>%
                                              mutate(phylum.fill = ifelse(phylum.x=="", phylum.y, phylum.x)) %>%
                                              mutate(class.fill = ifelse(class.x=="", class.y, class.x)) %>%
                                              mutate(order.fill = ifelse(order.x=="", order.y, order.x)) %>%
                                              mutate(family.fill = ifelse(family.x=="", family.y, family.x)) %>%
                                              mutate(genus.fill = ifelse(genus.x=="", genus.y, genus.x)) %>%
                                              mutate(species.fill = ifelse(species.x=="", species.y, species.x)) %>%
                                              select(row_no, asv, id_mer, domain.fill, phylum.fill, class.fill, order.fill, family.fill, genus.fill, species.fill)

UCout.denovo.mer$species.string <- stringr::str_split_fixed(gsub("d__|p__|c__|o__|f__|g__|s__", "", UCout.denovo.mer$species.fill), "_", 2)[,2]

UCout.denovo.mer$species.fill <- paste(UCout.denovo.mer$genus.fill, UCout.denovo.mer$species.string, sep="_")

UCout.denovo.mer <- UCout.denovo.mer %>% select(row_no, asv, id_mer, domain.fill, phylum.fill, class.fill, order.fill, family.fill, genus.fill, species.fill)
colnames(UCout.denovo.mer) <- c("row_no", "asv", "id_mer", "domain", "phylum", "class", "order", "family", "genus", "species")
poly_mer <- NULL
poly_mer <- list()
poly_mer$species <- UCout.denovo.mer[, .(asv ==first(asv), nParents = uniqueN(genus), genus = paste0(sort(unique(genus)), collapse= ", "), new = first(genus)), by = species][nParents > 1,]
poly_mer$genus <- UCout.denovo.mer[, .(asv ==first(asv), nParents = uniqueN(family), family = paste0(sort(unique(family)), collapse= ", "), new = first(family)), by = genus][nParents > 1,]
poly_mer$family <- UCout.denovo.mer[, .(asv ==first(asv), nParents = uniqueN(order), order = paste0(sort(unique(order)), collapse= ", "), new = first(order)), by = family][nParents > 1,]
poly_mer$order <- UCout.denovo.mer[, .(asv ==first(asv), nParents = uniqueN(class), class = paste0(sort(unique(class)), collapse= ", "), new = first(class)), by = order][nParents > 1,]
poly_mer$class <- UCout.denovo.mer[, .(asv ==first(asv), nParents = uniqueN(phylum), phylum = paste0(sort(unique(phylum)), collapse= ", "), new = first(phylum)), by = class][nParents > 1,]
poly_mer$phylum <- UCout.denovo.mer[, .(asv ==first(asv), nParents = uniqueN(domain), domain = paste0(sort(unique(domain)), collapse= ", "), new = first(domain)), by = phylum][nParents > 1,]

polyTaxaLog <- unlist(sapply(poly_mer, function(x) {
  if(nrow(x) > 1) {
    paste0(colnames(x)[1],
           " ",
           x[[1]],
           " has ",
           x[[3]],
           " parents: \"",
           x[[4]],
           "\", and has been assigned the ",
           colnames(x)[4],
           " of ",
           x[[2]],
           ": ",
           x[[5]])
  }
}), use.names = FALSE)
nTaxa <- length(polyTaxaLog)

#issue a warning if one or more taxa had more than one parent, write out logfile
if(nTaxa > 0) {
  warning(paste0(nTaxa,
                 " taxa had more than one parent, see the logfile \"./output/polyphyletics.log\" for details"),
          call. = FALSE)
  writeLines(polyTaxaLog,
             paste0("bac/clustering/polyphyletics2.log"))
}

#fix them
UCout.denovo.mer[, genus := first(genus), by = species]
UCout.denovo.mer[, family := first(family), by = genus]
UCout.denovo.mer[, order := first(order), by = family]
UCout.denovo.mer[, class := first(class), by = order]
UCout.denovo.mer[, phylum := first(phylum), by = class]
UCout.denovo.mer[, domain := first(domain), by = phylum]

UCout.denovo.mer$tax.vector <- paste(UCout.denovo.mer$domain, UCout.denovo.mer$phylum, UCout.denovo.mer$class, UCout.denovo.mer$order, UCout.denovo.mer$family, UCout.denovo.mer$genus, UCout.denovo.mer$species, sep=";")
UCout.denovo.mer$tax.vector1 <- paste(UCout.denovo.mer$asv, UCout.denovo.mer$domain, UCout.denovo.mer$phylum, UCout.denovo.mer$class, UCout.denovo.mer$order, UCout.denovo.mer$family, UCout.denovo.mer$genus, UCout.denovo.mer$species, sep=";")

UCout.denovo.mer.sub <- UCout.denovo.mer %>% select(asv, tax.vector1)

write.table(UCout.denovo.mer.sub, "bac/clustering/UCout_merged_annotations_bac_seedsonly.txt",row.names=FALSE, sep="\t", quote = FALSE)

ASV.table <- read.table("ASV_bac_counts.txt", header=T, row.names=1, sep="\t")
ASV.table$tax.vector <- UCout.denovo.mer$tax.vector

write.table(ASV.table, file="ASV_bac_counts_annotated.txt", sep="\t", col.names=NA, quote=F) #Kingdom to Species
