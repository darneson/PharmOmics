#
# Interface with GEO using SOFTtext
# Start by looking studies in GEO that have PubMed annotations
#
# Input:
#   input_1                     description
#
# Optional input:
#   input_2                     description
#   input_3                     description
#
# Written by Douglas Arneson 2016
#
getGDSs <- function(msea) {
  
  cat("\nPharmOmics Version: 02.26.2016\n")
  
  #load required libraries
  library(GEOquery)
  library(GEOmetadb)
  library(reutils)
  library(rentrez)
  library(plyr)
  library(stringr)
  require(stringr)
  
  asdf <- PharmOmics
  
  #FILE PATHS
  pipeline_path <- "/Users/darneson/Desktop/Box\ Sync/Doug_Macbook_Air/Yang_Lab/PharmOmics_Pipeline/"
  soft_files <- list.files(path=paste(pipeline_path,"Resources/gds_files",sep=""), pattern="*.soft.gz", full.names=T, recursive=FALSE)
  
  #get GEO MetaDB so can search locally
  if(file.exists(paste(pipeline_path,"/Resources/GEOmetadb.sqlite", sep=""))) {
    cat("\nThe GEO meta details already exist\n")
  } else {
    getSQLiteFile(destdir = paste(pipeline_path,"/Resources/",sep=""), destfile = "GEOmetadb.sqlite.gz")
  }
  
  #Connect to SQL database
  con <- dbConnect(SQLite(),paste(pipline_path,'/Resources/GEOmetadb.sqlite',sep=""))
  
  #Get number of GSEs on GEO --> 3867 GDSs on GEO
  all_GDSs <- dbGetQuery(con,'select * from gds')
  
  num_GDSs <- nrow(all_GDSs)
  
  cat(paste("\nThere are", num_GDSs, "GDSs currently on GEO\n"))
  
  #Get list of all GDS names
  GDS_names <- all_GDSs[,2]
  
  #start with small subset, delete this line later
  #num_GDSs <- 200
  
  categories <- NULL
  
  #see if the soft files for every GDS has been downloaded, otherwise download them
  if(!length(soft_files) == num_GDSs){
    
    cat(paste("\nStarting download of all necessary GEO GDS SOFT files"))
    
    #not all GDS files from GEOQuery exist, have to use "try" to get the loop to continue --> look into fixing this in the future
    #download all soft files for all GDSs from GEO, also get the different categories for the listed samples to determine which to
    #parse for the local copies
    #Got 3867 files, 20 have 0kb, so 3847 files to work with
    for(entry in 1:num_GDSs){
      try(current_gds <- getGEO(GDS_names[entry], destdir=paste(pipeline_path,"/Resources/gds_files",sep="")))
      try(current_gds <- getGEO(filename=paste(pipeline_path,"/Resources/gds_files/",GDS_names[entry],".soft.gz",sep="")))
      
      to_append <- NULL
      
      col_headers <- colnames(Columns(current_gds))
      to_append <- col_headers[ ! col_headers %in% categories]
      categories <- c(categories, to_append)
      
      if (entry%%10 == 0){
        cat(paste("\n",entry,"GEO GDS SOFT files have been downloaded out of",num_GDSs,"\n"))
      }
    }
    #write all the different descriptive categories of samples to a text file
    write(categories, file = paste(pipeline_path,'Resources/sample_categories.txt',sep=""), sep = '\t')
  }else{
    cat(paste("\nYou already have downloaded all necessary GEO GDS SOFT files, proceeding to the next step\n"))
  }
  
  GDS_substance <- paste(pipeline_path,"Resources/GDS_with_substance.txt",sep="")
  
  if (file.exists(GDS_substance)){
    cat(paste("\nAll GDSs have already been related to their corresponding substances\n"))
  }else{
    cat(paste("\nRelating all GDSs to their corresponding substances and writing to file\n"))
    #Initialize a start time to track how long the iterations will take
    ptm <- proc.time()
    for (i in 1:length(soft_files)){
      #print every 10 so we get an idea of how long it will take
      if (i%%10 == 0){
        time_elapsed <- (proc.time() - ptm)[3]
        remaining_time <- (length(soft_files)-i)*(time_elapsed/i)
        cat("\n", i, "GDSs have been processed in", time_elapsed, "seconds\n",
            "Approximately", remaining_time, "seconds still remain")
      }
        
      #do not conisder the GDS files that do not exist (i.e. GEOQuery had them, but when downloaded from GEO, they did not exist so they
      #have a size of 0KB
      #based on all the categories obtained from the categories file, we are only interested in the "agent" columns
      if (file.size(soft_files[i]) != 0){
        gds_file <- getGEO(filename=soft_files[i])
        
        #get what each sample was treated with
        sample_treatment <- unique(Columns(gds_file)$agent)
        
        #only interested in samples which have some sort of drug treatment
        if (!is.null(sample_treatment)){
          gds_name <- unlist(strsplit(soft_files[i], "/"))[11]
          gds_name <- unlist(strsplit(gds_name, "\\."))[1]
          sample_treatment <- data.frame(lapply(sample_treatment, as.character), stringsAsFactors=FALSE)
          to_write <- c(gds_name,sample_treatment)
    
          #Write the vector (GDS and corresponding substances) to a txt
          write.table(to_write, file = GDS_substance, sep = "\t", 
                      col.names = FALSE, row.names = FALSE, append=TRUE) 
        }
      }
    }
  }
}