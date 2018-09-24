library(data.table)

getLongerPolyA <- function(x){
  tmp <- gregexpr("A+", x, perl = TRUE)
  max(attributes(tmp[[1]])$match.length)
}

getFrequencies <- function(my.data.table){
  my.data.table[, lengths := nchar(V6)]
  my.data.table[, lengthsA := sapply(V6, getLongerPolyA)]
  as.data.frame(table(my.data.table$lengthsA))
}


human <- fread("Data/Homo_sapiens_polyA.data")
getFrequencies(human)

tetrahymena <- fread("Data/Tetrahymena_thermophila_polyA.data")
getFrequencies(tetrahymena)

plasmodium <- fread("Data/Plasmodium_falciparum_polyA.data")
getFrequencies(plasmodium)

