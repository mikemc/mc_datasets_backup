library(stringr)
library(tidyverse)

data_path <- "~/data/vcu_bias"

## Get the urls for downloading
accession_file <- file.path(data_path, "ncbi_accessions.csv")
if (!file.exists(accession_file)) {
    library(SRAdb)
    # Set this to the path of the SRAdb sqlite file
    sqlfile = file.path("~/data", "SRAmetadb.sqlite")
    file.info(sqlfile)
    sra_con <- dbConnect(SQLite(), sqlfile)
    # Get all the runs for the study SRP050185
    runs <- listSRAfile("SRP050185", sra_con, fileType = "sra") %>% as_tibble
    dbDisconnect(sra_con)
    # Save so we don't have to do this again
    readr::write_csv(runs, accession_file)
}
runs <- readr::read_csv(accession_file)

# For downloading via rsync, we need to replace the "ftp:" at the beginning of
# the FTP path with "rsync:"
runs <- runs %>%
    mutate(rsync = paste0('rsync:', str_sub(ftp, 5)))
runs$ftp %>% head
runs$rsync %>% head

# Download the files with rsync
dir.create(file.path(data_path, "reads"))
commands <- paste("rsync", "--progress", "--copy-links", "--times",
    "--verbose", runs$rsync, file.path(data_path, "reads"))
walk(commands, system)
# walk(commands[1:10], system)

# Convert to fastq
setwd(file.path(data_path, "reads"))
system("fastq-dump --gzip *.sra")
