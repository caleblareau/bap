options(warn=-1)

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tools)))
suppressMessages(suppressWarnings(library(stringr)))

"%ni%" <- Negate("%in%")

options(scipen=999)

# This function / script for munging WASP-style 
# QC stats

args <- commandArgs(trailingOnly = TRUE)
if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}

outfile <- args[i+1] # final resting places
directory <- args[i+2] # file path for the files

lf <- list.files(directory, pattern = "_stats.txt$", full.names = TRUE)

lapply(lf, read.table) %>% rbindlist() %>% data.frame() -> wholeDF
splitDF <- data.frame(str_split_fixed(as.character(wholeDF[,1]), "_", 3), stringsAsFactors = FALSE)
colnames(splitDF) <- c("chromosome", "what", "count")
splitDF$count <- as.numeric(as.character(splitDF$count))

splitDF%>% group_by(what) %>% summarise(total = sum(count)) %>% data.frame()-> sumDF
h1b <- sumDF[sumDF$what == "Haplotype1","total"] - sumDF[sumDF$what == "Haplotype1kept","total"]
h2b <- sumDF[sumDF$what == "Haplotype2","total"] - sumDF[sumDF$what == "Haplotype2kept","total"]

sumDF2 <- rbind(sumDF, data.frame(what = c("Haplotype1biased", "Haplotype2biased"),
                                  total = c(h1b, h2b)))
finalDF <- sumDF2[sumDF2$what %in% c("Discordant", "Haplotype1kept", "Haplotype2kept", "Haplotype1biased", "Haplotype2biased", "Unassigned"),]

write.table(finalDF, file = outfile, row.names = FALSE, col.names = FALSE,
            quote = FALSE, sep = "\t")
