library(data.table)
library(dplyr)

in_file <- "./authors.csv"
out_file <- "./authors_formatted.txt"

#--------------------------------------------------------------
## Read file
x <- fread(in_file, data.table = FALSE)

## Fill in blanks in author column
for (k in seq_len(nrow(x))[-1]) {
    if (x$Name[k] == "") x$Name[k] <- x$Name[k-1]
}

## Checks
# if (!all(grepl("\\.$", x$Affiliation))) stop("full stop needed at end of each affiliation")

i <- 0
last_name <- "zzz"
authors <- setNames(vector("list", length(unique(x$Name))), unique(x$Name))
## so the ordering of names here will be correct, now just need to get the institution numbers in correct order
institutions <- list()
for (this_author in names(authors)) {
    for (this_institution in x$Affiliation[x$Name == this_author]) {
        if (!this_institution %in% institutions) institutions <- c(institutions, this_institution)
        authors[[this_author]] <- c(authors[[this_author]], which(institutions == this_institution))
    }
}

## Format author list
af <- c()
for (ai in seq_len(length(authors))) {
    af[ai] <- sprintf("%s%s", names(authors)[ai], paste(authors[[ai]], collapse = ","))
}

## Write
writeLines(text = paste0(paste(af, collapse = ", "), "\n", paste0(seq_len(length(institutions)), ". ", institutions, collapse = "\n"), collapse = ""),
           con = out_file, useBytes = FALSE) # May need useBytes = TRUE
#--------------------------------------------------------------