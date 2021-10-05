rm(list = ls())

file_path <- commandArgs(trailingOnly = T)

x <- system(sprintf("cut -f2 %s", file_path), intern = T)
id <- which(nchar(x) == 0)
x <- x[2:(id - 1)]
x <- as.numeric(x)

cat(mean(x))

