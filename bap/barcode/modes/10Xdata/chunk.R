library(data.table)
library(stringr)

dt <- fread("zcat < 737K-cratac-v1.txt.gz", header = FALSE)
p1 <- sort(unique(substr(dt[[1]], 1, 9)))
p2 <- sort(unique(substr(dt[[1]], 10, 16)))

part_df <- data.frame(parts = c(p1, p2), 
           which = c(rep("1", length(p1)),
                     rep("2", length(p2))))
write.table(part_df, file = "737K-cratac-v1.parts.txt", sep = ",",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

rc <- function(s) chartr("ATGC","TACG",s)
strReverse <- function(x) sapply(lapply(strsplit(x, NULL), rev), paste, 
                                 collapse="") 
a <- "AGTGGAGTGTGCTGGC"
rc(strReverse(a)) %in% dt[[1]]
