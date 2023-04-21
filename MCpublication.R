# This script just prepares the data for the table Terhi needs

projFolder <- "/scratch/project_2002561/MastitisChallenge/RNASeq-Analysis"

fc.files <- list.files(file.path(projFolder, "FeatureCounts", "Annot"), pattern="*novaseq2_fc.txt.summary")

fc.files <- fc.files[-which(grepl("saureus", fc.files))]

fc.files <- fc.files[-which(grepl("6h", fc.files))]

# Get number of samplesw
  length(grep("cntrl1-0h", fc.files)) +
  length(grep("cntrl2-0h", fc.files)) +
  length(grep("cntrl3-0h", fc.files))

  length(grep("cntrl1-3h", fc.files)) +
  length(grep("cntrl2-3h", fc.files)) +
  length(grep("cntrl3-3h", fc.files))

  length(grep("cntrl1-24h", fc.files)) +
  length(grep("cntrl2-24h", fc.files)) +
  length(grep("cntrl3-24h", fc.files))

  length(grep("coli1-3h", fc.files)) +
  length(grep("coli2-3h", fc.files)) +
  length(grep("coli3-3h", fc.files))

  length(grep("coli1-24h", fc.files)) +
  length(grep("coli2-24h", fc.files)) +
  length(grep("coli3-24h", fc.files))
  
  
alignmentFiles <- list.dirs(file.path(projFolder, "BAM"), recursive=FALSE)
alignmentFiles <- alignmentFiles[grep("novaseq2", alignmentFiles)]
alignmentFiles <- alignmentFiles[-grep("saureus", alignmentFiles)]
alignmentFiles <- alignmentFiles[-grep("6h", alignmentFiles)]

alignments <- list()
for(i in 1:length(alignmentFiles)){
  fileIn <- alignmentFiles[i]
  sampleName <- strsplit(fileIn, "/")[[1]][7]
  fileIn <- paste0(fileIn,"/", paste0(sampleName, "_Log.final.out"))
  alignments[[i]] <- readLines(fileIn)
  names(alignments)[i] <- sampleName
}

trimmedReads <- sapply(strsplit(sapply(alignments, "[", 6), "\\t"),"[", 2)
alignedReads <- sapply(strsplit(sapply(alignments, "[", 10), "\\t"),"[", 2)

mean(as.numeric(trimmedReads[grep("cntrl[123]-0h", names(trimmedReads))]))/1000000
mean(as.numeric(trimmedReads[grep("cntrl[123]-3h", names(trimmedReads))]))/1000000
mean(as.numeric(trimmedReads[grep("cntrl[123]-24h", names(trimmedReads))]))/1000000
mean(as.numeric(trimmedReads[grep("coli[123]-3h", names(trimmedReads))]))/1000000
mean(as.numeric(trimmedReads[grep("coli[123]-24h", names(trimmedReads))]))/1000000

options(scipen=999)

alignedReads <- gsub("%", "", alignedReads)
mean(as.numeric(alignedReads[grep("cntrl[123]-0h", names(alignedReads))]))
mean(as.numeric(alignedReads[grep("cntrl[123]-3h", names(alignedReads))]))
mean(as.numeric(alignedReads[grep("cntrl[123]-24h", names(alignedReads))]))
mean(as.numeric(alignedReads[grep("coli[123]-3h", names(alignedReads))]))
mean(as.numeric(alignedReads[grep("coli[123]-24h", names(alignedReads))]))


assigned <- c()
for(i in 1:length(fc.files)){
  fcIn <- readLines(file.path(projFolder, "FeatureCounts", "Annot", fc.files[i]))
  tmp <- as.numeric(sapply(strsplit(fcIn[-1],"\t"), "[",2))
  assigned[i] <- tmp[10]/sum(tmp[-2])
  names(assigned)[i] <- gsub("_fc.txt.summary", "", fc.files[i])
}
mean(as.numeric(assigned[grep("cntrl[123]-0h", names(assigned))]))
mean(as.numeric(assigned[grep("cntrl[123]-3h", names(assigned))]))
mean(as.numeric(assigned[grep("cntrl[123]-24h", names(assigned))]))
mean(as.numeric(assigned[grep("coli[123]-3h", names(assigned))]))
mean(as.numeric(assigned[grep("coli[123]-24h", names(assigned))]))
