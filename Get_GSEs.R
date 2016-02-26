#This package allows interface with GEO using SOFTtext
#asf
#Here we will only be looking at drugs as annotated by Pubmed in the Substances subcategory
#Start by loading all required packages
library(GEOquery)
library(GEOmetadb)
library(reutils)
library(stringr)
require(stringr)

#GET GEO METADB SO CAN SEARCH LOCALLY
if(file.exists("/Users/darneson/Desktop/Box\ Sync/Doug_Macbook_Air/Yang_Lab/PharmOmics_Pipeline/Resources/GEOmetadb.sqlite")) {
  print("The GEO meta details already exist")
} else {
  getSQLiteFile(destdir = "/Users/darneson/Desktop/Box\ Sync/Doug_Macbook_Air/Yang_Lab/PharmOmics_Pipeline/Resources/", destfile = "GEOmetadb.sqlite.gz")
}

#Connect to SQL database
con <- dbConnect(SQLite(),'/Users/darneson/Desktop/Box\ Sync/Doug_Macbook_Air/Yang_Lab/PharmOmics_Pipeline/Resources/GEOmetadb.sqlite')

#Get number of GSEs on GEO --> 56,174 GSEs on GEO
rs <- dbGetQuery(con,'select * from gse')
print(nrow(rs))

#Pubmed IDs are in column 7
pmID <- rs[,7]

#Here we want to get all GEO GSEs that have pubmed IDs
onlyIDs <- vector()
onlyNAs <- vector()

for (i in 1:length(pmID)) {
  if (!is.na(pmID[i]))
    onlyIDs <- c(onlyIDs,pmID[i])
  else onlyNAs <- c(onlyNAs,pmID[i])
}

#32425 GSEs have pubmed IDs
print(length(onlyIDs))

#23749 GSEs do not have pubmed IDs
print(length(onlyNAs))

#Their sum adds to 56,174 which is what we started with, however some of the pubmed IDs are repeated (multiple GSEs for a pubmed publication)

#Find Unique pubmed IDs
uniqueIDs <- unique(onlyIDs)
#There are 24128 Unique pubmed IDs
print(length(uniqueIDs))

#Initialize Start Time
ptm <- proc.time()

#Iterate through all 32425 GSEs that have pubmed IDs (for now doing a subset of 200)
for(entry in 1:length(uniqueIDs)){
  #print every 100 so we get an idea of how long it will take
  if (entry%%100 == 0){
    print(entry)
    print(proc.time() - ptm)
  }
  pmid <- uniqueIDs[entry]
  
  #convert pubmed ID to a string for use in efetch command  
  pmid <- toString(pmid)
  
  #Use pubmed ID to get medline info (including MeSH terms)
  medline <- efetch(pmid, db = "pubmed", "medline")
  
  #Write the fetched pubmed info to a text file b/c the whole entry cannot be accessed directly in R
  write(content(medline), file = "tempSubstances.txt")
  
  #Read the Medline Data back into R
  res <- readLines("tempSubstances.txt")
  
  #Initialize empty vector that will hold pmID and the Substance terms
  id_MHterms <- vector()
  
  #Pull only the PMID, RN(drug) from the Medline info from Pubmed
  for (j in 1:length(res)){
    #line <- gsub(" ","",res[j])
    line <- res[j]
    line <- strsplit(line, " - ")
    line <- unlist(line)
    if (length(line) > 1){
      if (!is.na(line[1]) && (line[1] == "PMID"))
        id_MHterms <- c(id_MHterms,line[2])
      else 
      if (!is.na(line[1]) && (line[1] == "RN ")){ 
        if (length(line) == 2){
          id_MHterms <- c(id_MHterms,gsub("\\*","",line[2]))
        }
        if (length(line) > 2){
          id_MHterms <- c(id_MHterms,gsub("\\*","",(paste((line[2:length(line)]), sep="", collapse="-"))))
        }
      }
    }
  }  
  #Transpose vector so we can append vector as a row to a csv
  transposed <- as.matrix(t(id_MHterms))
  
  #Write the vector (single pubmed entry corresponding to GSE) to a csv
  write.table(transposed, file = "PMID_Substance_terms_all_newest.txt", sep = "\t", 
              col.names = FALSE, row.names = FALSE, append=TRUE)
}

#get running time
print(proc.time() - ptm)