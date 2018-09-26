# generate riboseq occupancy based on A-site GWIPS global tracks

library(rtracklayer)
library(data.table)
library(biomaRt)
library(genomeIntervals)


convertToGRanges <- function(polyA.df, before = 50, after = 50){
  GRanges(seqnames = polyA.df$chr, 
          ranges = IRanges(start = polyA.df$start - before,
                           end = polyA.df$start + after,
                           names = polyA.df$chr),
          polyA.length =  nchar(polyA.df$seq),
          direction = polyA.df$strand)
}
# elongating ribosomes, A-site mapping
plasmodium.bw <- "./GWIPS/pf_Global_RiboProElong.08_11_2017.bw"
# elongating ribosomes, A-site mapping
human.bw <- "./GWIPS/All_human.bw"

# plasmodium, stage specific
plasmodium.ET.bw <- "./GWIPS/EarlyTrophozoite_RiboPro.bw"
plasmodium.R.bw <- "./GWIPS/Ring_RiboPro.bw"
plasmodium.LT.bw <- "./GWIPS/LateTrophozoite_RiboPro.bw"
plasmodium.S.bw <- "./GWIPS/Schizont_RiboPro.bw"
plasmodium.M.bw <- "./GWIPS/Merozoite_RiboPro.bw"

# convert to ranges object
human.rs <- rtracklayer::import.bw(human.bw)
plasmodium.rs <- rtracklayer::import.bw(plasmodium.bw)
plasmodium.ET.rs <- rtracklayer::import.bw(plasmodium.ET.bw)
plasmodium.R.rs <- rtracklayer::import.bw(plasmodium.R.bw)
plasmodium.LT.rs <- rtracklayer::import.bw(plasmodium.LT.bw)
plasmodium.S.rs <- rtracklayer::import.bw(plasmodium.S.bw)
plasmodium.M.rs <- rtracklayer::import.bw(plasmodium.M.bw)



#load polyA track information
human.polyA <- data.table::fread("./Data/Homo_sapiens_polyA.data")
names(human.polyA) <- c("transcript", "internal_start", "internal_end", "start", "end", "seq", "empty")

plasmodium.polyA <- data.table::fread("./Data/Plasmodium_falciparum_polyA.data")
names(plasmodium.polyA) <- c("transcript", "internal_start", "internal_end", "start", "end", "seq", "empty")

# use ensembl to identify gene boundaries
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart=ensembl)

map2chr <- getBM(attributes = c("chromosome_name", "strand", "ensembl_transcript_id"), filters = "ensembl_transcript_id", values = human.polyA$transcript, mart = ensembl, uniqueRows = FALSE)

human.polyA.full <- merge(human.polyA, map2chr, by.x = "transcript", by.y = "ensembl_transcript_id")
human.polyA.full$chr <- paste("chr", human.polyA.full$chromosome_name, sep = "")
polyA.human.gr <- unique(convertToGRanges(human.polyA.full, before = 50))
polyA.human.df.psite <- data.frame(location=c(), occupancy=c(), sum=c(), average=c())

#generate data for +/- 50nt from start of polyA

for (i in seq_along(polyA.human.gr)){
  message(i)
  tmp <- as.data.frame(polyA.human.gr[i])
  dat.vec <- rep(0, 101)
  names(dat.vec) <- c(tmp$start:tmp$end)
  tmp.gr <- as.data.frame(subsetByOverlaps(human.rs, polyA.human.gr[i]))
  if (nrow(tmp.gr) > 5){
    for (i in c(tmp$start:tmp$end)){
      for (ii in c(1:nrow(tmp.gr))){
        if (i >= tmp.gr$start[ii] && i <= tmp.gr$end[ii]){
          dat.vec[as.character(i)] <- tmp.gr$score[ii]
        }
      }
    }
    #dat.vec[as.character(tmp.gr$end)] <- tmp.gr$score
    if (tmp$direction < 0){
      dat.vec <- rev(dat.vec)
    }
    
    polyA.human.df.psite <- rbind(polyA.human.df.psite, data.frame(location=c(-50:50), 
                                                                   occupancy=dat.vec,
                                                                   sum=rep(sum(dat.vec), 101),
                                                                   average=rep(mean(dat.vec), 101)))
  }
}


# plot polyA in human
ggplot(polyA.human.df.psite, aes(x=location, y = occupancy, group = location)) + 
  geom_point(position = "jitter", alpha=.1) + theme_light() + 
  stat_summary(fun.data = "mean_cl_normal", geom = "crossbar") + 
  coord_cartesian(xlim=c(-20,40), ylim=c(0,100))

ggplot(subset(polyA.human.df.psite, average >=8.76), aes(x=location, y = occupancy, group = location)) +
  geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim=c(0,200)) + 
  theme(axis.title = element_text(face = "bold", size = rel(2)), 
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "lightgrey"), 
        panel.border = element_rect(fill = NA, color = "grey")) + 
  xlab("NT from start of polyA segment") + 
  ylab("Mean occupancy of elongating ribosomes at A-site")


# read annotation for the genome assembly version used in GWIPS database
gff <- readGff3("PlasmoDB-26_Pfalciparum3D7.gff")
aliases <- getGffAttribute(gff, "Alias")

plasmodium.chrom.names <- c()
plasmodium.strand <- c()
for (i in c(1:nrow(plasmodium.polyA))){
  transcript.name <- stri_match(str = plasmodium.polyA$transcript[i], regex = "(.+):")
  transcript.name <- transcript.name[1,2]
  position <- grep(transcript.name, aliases)
  if(length(position) > 1) message(paste(transcript.name, position, sep = " "))
  seq.name <- annotation(gff[position[1]])
  found.numbers <- stri_extract_all(str = as.character(seq.name$seq_name[1]), regex = "\\d+")
  if (length(found.numbers[[1]]) == 4){
    found.numbers[[1]][3] <- stri_replace(found.numbers[[1]][3], "", regex = "^0")
    plasmodium.chrom.names[i] <- paste("chr", found.numbers[[1]][3], sep = "")
    plasmodium.strand[i] <- ifelse(as.character(seq.name$strand[1]) == "+", 1, -1)
  }
}

plasmodium.chrom.names[which(is.na(plasmodium.chrom.names))] <- "empty"

plasmodium.polyA.full <- plasmodium.polyA
plasmodium.polyA.full$chr <- plasmodium.chrom.names
plasmodium.polyA.full$strand <- plasmodium.strand


polyA.plasmodium.gr <- unique(convertToGRanges(plasmodium.polyA.full, before = 50))


iterateOverRS <- function(bigwigobject, input.gr) {
  output.df <-
    data.frame(
      location = c(),
      occupancy = c(),
      sum = c(),
      average = c(),
      polyA.length = c()
    )
  
  for (i in seq_along(input.gr)) {
    message(i)
    tmp <- as.data.frame(input.gr[i])
    dat.vec <- rep(0, 101)
    names(dat.vec) <- c(tmp$start:tmp$end)
    tmp.gr <-
      as.data.frame(subsetByOverlaps(bigwigobject, input.gr[i]))
    if (nrow(tmp.gr) > 5) {
      for (i in c(tmp$start:tmp$end)) {
        for (ii in c(1:nrow(tmp.gr))) {
          if (i >= tmp.gr$start[ii] && i <= tmp.gr$end[ii]) {
            dat.vec[as.character(i)] <- tmp.gr$score[ii]
          }
        }
      }
      if (tmp$direction < 0) {
        dat.vec <- rev(dat.vec)
      }
      
      output.df <- rbind(
        output.df,
        data.frame(
          location = c(-50:50),
          occupancy = dat.vec,
          sum = rep(sum(dat.vec), 101),
          average = rep(mean(dat.vec), 101),
          polyA.length = rep(tmp$polyA.length, 101)
        )
      )
    }
  }
  output.df
}
# generate data for +/- 50nt around start of polyA
polyA.plasmodium.df.ET.asite <- iterateOverRS(plasmodium.ET.rs, polyA.plasmodium.gr)
polyA.plasmodium.df.R.asite <- iterateOverRS(plasmodium.R.rs, polyA.plasmodium.gr)
polyA.plasmodium.df.LT.asite <- iterateOverRS(plasmodium.LT.rs, polyA.plasmodium.gr)
polyA.plasmodium.df.S.asite <- iterateOverRS(plasmodium.S.rs, polyA.plasmodium.gr)
polyA.plasmodium.df.M.asite <- iterateOverRS(plasmodium.M.rs, polyA.plasmodium.gr)



mean.df <- summary(random.plasmodium.df.asite$occupancy)["Mean"]

ggplot(subset(polyA.plasmodium.df.asite, polyA.length<22 & average >= mean.df), aes(x=location, y = occupancy, group = location)) +
geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim=c(0,200)) + 
  theme(axis.title = element_text(face = "bold", size = rel(2)), 
        panel.background = element_rect(fill = "white"),
        panel.grid = element_line(colour = "lightgrey"), 
        panel.border = element_rect(fill = NA, color = "grey")) + 
  xlab("NT from start of polyA segment") + 
  ylab("Mean occupancy of elongating ribosomes at A-site")



# randomize intervals
n.intervals <- nrow(plasmodium.polyA)
genes <- gff[which(annotation(gff)$type == "gene")]

random.plasmodium.chrom.names <- c()
random.plasmodium.strand <- c()
random.plasmodium.start <- c()
for (i in c(1:nrow(plasmodium.polyA))){

  position <- sample(1:nrow(genes), 1)
  seq.name <- annotation(genes[position])
  found.numbers <- stri_extract_all(str = as.character(seq.name$seq_name[1]), regex = "\\d+")
  if (length(found.numbers[[1]]) == 4){
    found.numbers[[1]][3] <- stri_replace(found.numbers[[1]][3], "", regex = "^0")
    random.plasmodium.chrom.names[i] <- paste("chr", found.numbers[[1]][3], sep = "")
    random.plasmodium.strand[i] <- ifelse(as.character(seq.name$strand[1]) == "+", 1, -1)
  }
  gene.coords <- as.vector(genes[position])
  random.plasmodium.start[i] <- sample(gene.coords[1]:gene.coords[2], 1)
  
}

random.plasmodium.chrom.names[which(is.na(random.plasmodium.chrom.names))] <- "empty"

random.plasmodium.df <- data.frame(chr = random.plasmodium.chrom.names,
                                   strand = random.plasmodium.strand,
                                   start = random.plasmodium.start)
random.plasmodium.df$seq <- rep("NNN", nrow(plasmodium.polyA))

random.plasmodium.gr <- convertToGRanges(random.plasmodium.df)

random.plasmodium.df.asite <- iterateOverRS(plasmodium.rs, random.plasmodium.gr)
random.plasmodium.df.ET.asite <- iterateOverRS(plasmodium.ET.rs, random.plasmodium.gr)
random.plasmodium.df.R.asite <- iterateOverRS(plasmodium.R.rs, random.plasmodium.gr)
random.plasmodium.df.LT.asite <- iterateOverRS(plasmodium.LT.rs, random.plasmodium.gr)
random.plasmodium.df.S.asite <- iterateOverRS(plasmodium.S.rs, random.plasmodium.gr)
random.plasmodium.df.M.asite <- iterateOverRS(plasmodium.M.rs, random.plasmodium.gr)


