setwd("/Users/darneson/Desktop/Box Sync/Doug_Macbook_Air/Yang_Lab/PharmOmics_Pipeline/Scripts")
source("PharmOmics_2.R")

PharmOmics <- list()
GDS_substances <- getGDSs(PharmOmics)
UMLS_java <- writeUMLSjava(GDS_substances)

setwd("/Users/darneson/Desktop/Box Sync/Doug_Macbook_Air/Yang_Lab/PharmOmics_Pipeline/Resources/UMLS/utsapi2_0")
system("javac -g GDS_UMLS.java")
system("java GDS_UMLS darneson d4Av#otr 2015AB")
setwd("/Users/darneson/Desktop/Box Sync/Doug_Macbook_Air/Yang_Lab/PharmOmics_Pipeline/Scripts")

