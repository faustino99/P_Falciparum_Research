#Gene Set Formulation
#This script will create gene sets as tab delimited text files and will also create GMT format files for use in gene set enrichment analysis.

#step 1: Create a tab delimited text file of the csv table of Plasmodium falciparum's genes taken from APILOC's Github repository.
a<-read.table(file="Plasmodium_falciparum.csv", sep="\t", stringsAsFactors=FALSE, header=TRUE, strip.white = TRUE)


###########Updating Gene Symbols############
#make list of just first column
p<-a[[2]]

#read gene aliases table taken from plasmodb
q<-read.table(file = "Pf3D7_gene_aliases.txt", sep = "\t", stringsAsFactors = FALSE,fill = TRUE)

#Change gene names and add it to the original table
for (y in 2:4) {
  for (n in 1:length(q[[1]])) {
    if (q[n,y] != "") {
      p[grep(q[n,y],p,ignore.case = TRUE)]<-unlist(q[n,1])
    }  
  }
}
a[[2]]<-p

###############################
# make .chip file w/ updated gene symbols
############################### 
v<-read.csv(file = "Pf_microarray_genesymbols.csv",header = TRUE,skip = 23)
w<-v[,c(1,15)]
z<-as.vector(w[[2]])
for (y in 2:4) {
  for (n in 1:length(q[[1]])) {
    if (q[n,y] != "") {
      z[grep(q[n,y],z)]<-unlist(q[n,1])
    }
  }
}
w[[2]]<-z
w[[3]]<-v[[14]]
names(w)<-c("Probe Set ID","Gene Symbol","Gene Title")
write.table(w, file = "pfchip.chip", sep = "\t",quote = FALSE,row.names = FALSE)

###   make sure to replace "---" with NA   ###

###############################################

#step 2:create gene sets based on localization

#2a) narrow table down to the localization columns
localization_column<-a[[4]]

#2c) add plasmoDB ID column to localization list to make dataframe
localization_dataframe1<-data.frame(a[[2]], localization_column)
names(localization_dataframe1)<-c("PlasmoDB_ID", "Localization") #label columns

#2d) get rid of blank rows
localization_list<- localization_dataframe1[localization_dataframe1$PlasmoDB_ID!="", ]
localization_list<- localization_list[localization_list$Localization!="", ]

#2e) trim trailing white space on PlasmoDB IDs
unlist<-trimws(unlist(localization_list[1]), which = "right")
localization_list[1]<-list(unlist)



######################
#initial strategy for creating genesets, less efficient than current code
#factor the list
##b<-factor(unfactored_localization_list)
#create levels list 
##levels<-levels(localization_list[[2]])
##levels<-data.frame(levels)
#make gene set for apical based off levels
##apical_geneset<-localization_list[as.integer(localization_list[[2]])<=25 & as.integer(localization_list[[2]])>=5, ]
######################


#2f)better strategy for making genesets using grep

#Apical
#select all entries with "apical" in their localization description
apical_geneset<-localization_list[grep("apical|rhoptry bulb|merozoite surface|rhoptry neck", localization_list$Localization, ignore.case = TRUE), ] 
#get rid of "not apical" entries
apical_geneset<-apical_geneset[grep("not apical|not merosoite surface|not rhoptry neck", apical_geneset$Localization, ignore.case = TRUE, invert=TRUE), ] 
apical_geneset<-apical_geneset[grep("no matching", apical_geneset$PlasmoDB_ID, invert = TRUE), ]
#remove duplicaates
apical_geneset<-apical_geneset[!duplicated(apical_geneset[1]), ]
#export to tsv file
write.table(apical_geneset[[1]], file = "localization_gene_data/Apical_Geneset", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

#Endoplasmic reticulum & Golgi Aparatus
ER_golgi_geneset<-localization_list[grep("Er |ER |endoplasmic reticulum|golgi|Golgi", localization_list$Localization), ]
ER_golgi_geneset<-ER_golgi_geneset[grep("not ER |not golgi|Not golgi|not Golgi", ER_golgi_geneset$Localization, invert = TRUE), ]
ER_golgi_geneset<-ER_golgi_geneset[!duplicated(ER_golgi_geneset[1]), ]
write.table(ER_golgi_geneset[[1]], file = "localization_gene_data/ER_Golgi_Geneset", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

#Nucleus
nucleus_geneset<-localization_list[grep("nucleus|nuclear", localization_list$Localization, ignore.case = TRUE), ]
nucleus_geneset<-nucleus_geneset[grep("not nucleus|near nucleus|not nuclear", nucleus_geneset$Localization, ignore.case = TRUE, invert = TRUE), ]
nucleus_geneset<-nucleus_geneset[!duplicated(nucleus_geneset[1]), ]
write.table(nucleus_geneset[[1]], file = "localization_gene_data/Nucleus_Geneset", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

#Parasitophorous vacuole (including PVM)
PV_geneset<- localization_list[grep("parasitophorous vacuole|pv", localization_list$Localization, ignore.case = TRUE), ]
PV_geneset<- PV_geneset[grep("not pv",PV_geneset$Localization, ignore.case = TRUE, invert = TRUE), ]
PV_geneset<- PV_geneset[!duplicated(PV_geneset[1]), ]
write.table(PV_geneset[[1]], file = "localization_gene_data/PV_Geneset", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

#Mitochondrion
mitochondrion_geneset<- localization_list[grep("mitochondri",localization_list$Localization, ignore.case = TRUE), ]
mitochondrion_geneset<- mitochondrion_geneset[grep("not mitochondri",mitochondrion_geneset$Localization, ignore.case = TRUE, invert = TRUE), ]
mitochondrion_geneset<- mitochondrion_geneset[!duplicated(mitochondrion_geneset[1]), ]
write.table(mitochondrion_geneset[[1]], file = "localization_gene_data/Mitochondrion_Geneset", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

#Parasite Plasma Membrane
PPM_geneset<- localization_list[grep("PPM|PM|plasma membrane|Sporozoite surface|surrounding intracellular|parasite membrane",localization_list$Localization, ignore.case = TRUE), ]
PPM_geneset<- PPM_geneset[grep("not PPM", PPM_geneset$Localization, ignore.case = TRUE, invert = TRUE), ]
PPM_geneset<- PPM_geneset[grep("lots", PPM_geneset$PlasmoDB_ID, invert = TRUE), ]
PPM_geneset<- PPM_geneset[!duplicated(PPM_geneset[1]), ]
write.table(PPM_geneset[[1]], file = "localization_gene_data/PPM_Geneset", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

#Cytoplasm
cytoplasm_geneset<- localization_list[grep("cytoplasm|cytosol|peripheral|intracellular inclusion", localization_list$Localization, ignore.case = TRUE), ]
cytoplasm_geneset<- cytoplasm_geneset[grep("not cytoplasm|not erythrocyte cytoplasm|not RBC cytoplasm|not cytosol",cytoplasm_geneset$Localization, ignore.case = TRUE, invert = TRUE), ]
cytoplasm_geneset<- cytoplasm_geneset[!duplicated(cytoplasm_geneset[1]), ]
write.table(cytoplasm_geneset[[1]], file = "localization_gene_data/Cytoplasm_Geneset", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

#Apicoplast
apicoplast_geneset<- localization_list[grep("apicoplast", localization_list$Localization, ignore.case = TRUE ), ]
apicoplast_geneset<- apicoplast_geneset[grep("not apicoplast", apicoplast_geneset$Localization, ignore.case = TRUE, invert = TRUE), ]
apicoplast_geneset<- apicoplast_geneset[!duplicated(apicoplast_geneset[1]), ]
write.table(apicoplast_geneset[[1]], file = "localization_gene_data/Apicoplast_Geneset", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

#Inner Membrane Complex #at 15 genes
IMC_geneset<- localization_list[grep("inner membrane|IMC",localization_list$Localization, ignore.case = TRUE), ]
IMC_geneset<- IMC_geneset[grep("not IMC", IMC_geneset$Localization, ignore.case = TRUE, invert = TRUE), ]
IMC_geneset<- IMC_geneset[!duplicated(IMC_geneset[1]), ]
write.table(IMC_geneset[[1]], file = "localization_gene_data/IMC_Geneset", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

#Food Vacuole
food_vacuole_geneset<- localization_list[grep("food vacuole|FV", localization_list$Localization, ignore.case = TRUE), ]
food_vacuole_geneset<- food_vacuole_geneset[grep("not food vacuole|not FV", food_vacuole_geneset$Localization, ignore.case = TRUE, invert = TRUE), ]
food_vacuole_geneset<- food_vacuole_geneset[!duplicated(food_vacuole_geneset[1]), ]
write.table(food_vacuole_geneset[[1]], file = "localization_gene_data/Food_Vacuole_Geneset", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

#Exported
exported_geneset<- localization_list[grep("exported|erythrocyte|cleft like|maurer's cleft", localization_list$Localization, ignore.case = TRUE), ]
exported_geneset<- exported_geneset[grep("not erythrocyte|not maurer's cleft", exported_geneset$Localization, invert = TRUE, ignore.case = TRUE), ]
exported_geneset<- exported_geneset[!duplicated(exported_geneset[1]), ]
write.table(exported_geneset[[1]], file = "localization_gene_data/Exported_Geneset", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)


#3 create a GMT file of all the genesets

#a) download GSEA base package if not already installed

#how to download:
#source("https://bioconductor.org/biocLite.R")
#biocLite("GSEABase")

library(GSEABase)

#b) Create Genesets for all localizations
a<-GeneSet(setName= "apical", apical_geneset[,1])
b<-GeneSet(setName= "apicoplast", apicoplast_geneset[,1])
c<-GeneSet(setName= "cytoplasm", cytoplasm_geneset[,1])
d<-GeneSet(setName= "ER_golgi", ER_golgi_geneset[,1])
e<-GeneSet(setName= "exported", exported_geneset[,1])
f<-GeneSet(setName= "food_vacuole", food_vacuole_geneset[,1])
g<-GeneSet(setName= "IMC", IMC_geneset[,1])
h<-GeneSet(setName= "mitochondrion", mitochondrion_geneset[,1])
i<-GeneSet(setName= "nucleus", nucleus_geneset[,1])
j<-GeneSet(setName= "PPM", PPM_geneset[,1])
k<-GeneSet(setName= "PV", PV_geneset[,1])


#c) make gene set collection
gsc<-GeneSetCollection(a,b,c,d,e,f,g,h,i,j,k)

#d) export GMT file
toGmt(gsc, con = "localization_gene_data/Location.gmt")

# Step 3: create tab delimited file of transmission expression data of the Plasmodium Falciparum during its 48 hour intraerythrocytic developmental cyle

# read 3d7_smoothed.txt taken from Derisi microarray expression data in PlamsoDB
a<-read.table(file="transmission_expression.txt", sep = "\t", quote = "",header = TRUE, stringsAsFactors = TRUE)

# a) initialize then create columns for standard deviation and mean + 1.25 SD values for each protein
a$"SD"<- NA
a$"Mean+1.25_SD"<- NA
for (gene in 1:(length(a[[1]]))) {
  a[gene, 55]<- sd(a[gene, 2:54]) 
  a[gene, 56]<- a[gene, 55] * 1.25 + mean(unlist(a[gene, 2:54]))
}

#b) histogram that plots SD values for all the proteins
png(filename = "Images/Derisi_Histogram.png")
hist(a[,55], breaks = 100, xlab = "standard deviation", main = "Histogram of Standard Deviation", col = "dark green")
dev.off()

#c) for loop goes through each protein's timepoint and determines whether it has a peak expression value (x => mean+1.25 SD)
# any non-peak values are deleted from the table
for (x in 1:length(a[[1]])) {
  for (y in 2:(length(a[1,])-2)) {
    if (!isTRUE(a[x,y]>=a[x,56])) {
      a[x,y]<- NA
    } 
  }  
}

# d) save all of the peak expression protein lists for each time point as a GeneSet and as computer tsv files
#e) Create list of the lengths of each time point protein list

y<- NULL
Derisi_transmission_genesets_1.25SD<- NULL
library(GSEABase) # load GSEA Base package (if not already loaded)

for (x in 2:54) {
  y<- append(y, length(a[!is.na(a[,x]),1]))
  Derisi_transmission_genesets_1.25SD<- append(Derisi_transmission_genesets_1.25SD,GeneSet(setName = paste0("TP ",x-1), as.character(a[!is.na(a[,x]),1])))
  write.table(a[!is.na(a[,x]),1], file = paste0("gene_expression_data/Derisi/Derisi_1.25_SD_TP",x-1), sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
}

#f) create plot of the gene list lengths over time in order to determine moments of peak activity within the cell cycle
x<- c(1:53)
plot(x,y,ylab = "# of genes", xlab = "Hours", type = "b")

#g) save GeneSets to a GMT file
Derisi_1.25SD_GSC<-GeneSetCollection(Derisi_transmission_genesets_1.25SD)
toGmt(Derisi_1.25SD_GSC, con = "gene_expression_data/Derisi/Derisi_transmission_Expression_1.25_SD_Timepoints.gmt")

##########
### Step 3 1/2: do same as above only with peak expression cutoff being mean + 1.5 SD ###
##########

b<-read.table(file="transmission_expression.txt", sep = "\t", quote = "",header = TRUE, stringsAsFactors = TRUE)
b$"SD"<- NA
b$"Mean+1.5_SD"<- NA

for (gene in 1:(length(b[[1]]))) {
  b[gene, 55]<- sd(b[gene, 2:54]) 
  b[gene, 56]<- b[gene, 55] * 1.5 + mean(unlist(b[gene, 2:54]))
}

for (x in 1:length(b[[1]])) {
  for (y in 2:(length(b[1,])-2)) {
    if (!isTRUE(b[x,y]>=b[x,56])) {
      b[x,y]<- NA
    } 
  }  
}

y2<- NULL
Derisi_transmission_genesets_1.5SD<- NULL

for (x in 2:54) {
  y2<- append(y2, length(b[!is.na(b[,x]),1]))
  Derisi_transmission_genesets_1.5SD<- append(Derisi_transmission_genesets_1.5SD,GeneSet(setName = paste0("TP ",x-1), as.character(b[!is.na(b[,x]),1])))
  write.table(b[!is.na(b[,x]),1], file = paste0("gene_expression_data/Derisi/Derisi_1.5_SD_TP",x-1), sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
}

Derisi_1.5SD_GSC<-GeneSetCollection(Derisi_transmission_genesets_1.5SD)
toGmt(Derisi_1.5SD_GSC, con = "gene_expression_data/Derisi/Derisi_transmission_Expression_1.5_SD_Timepoints.gmt")

#create plot comparing two different peak expression genesets
y<- NULL
for (x in 2:54) {
  y<- append(y, length(a[!is.na(a[,x]),1]))
}

x<- c(1:53)
png(filename = "Images/Derisi_peak_expression_sizes.png")
plot(x,y2,ylab = "# of genes", xlab = "Hours", col = "blue", pch=18, type= "b",ylim = c(0,1250), main = "Number of Peak Expression Genes Over Time")
points(x,y, pch=20,  type = "b")
abline(h=500,lty=2)
abline(h=15, lty=2)
legend("topright", bty = "o", legend = c("Mean+1.25 SD", "Mean+1.5 SD"), lty = c(1,1), col = c("black","blue"), title = expression(bold("Peak Expression Value")))
dev.off()

#############Results from the next experiment's data did not come out as expected #############
####Commenting out this experiment since the results are not what we want ######

# #4 Do the same as step 3 but for a different data table from a different experiment doing the same thing as the previous one
# a <- read.table(file = "profiles_smoothed.txt", quote = "",stringsAsFactors = TRUE, header = TRUE, sep = "\t")
# 
# a$"SD"<- NA
# a$"Mean+1.25_SD"<- NA
# 
# for (gene in 1:(length(a[[1]]))) {
#   a[gene, 50]<- sd(a[gene, 2:49])
#   a[gene, 51]<- a[gene, 50] * 1.25 + mean(unlist(a[gene, 2:49]))
# }
# 
# hist(a[,50], breaks = 100, xlab = "standard deviation", main = "Histogram of Standard Deviation", col = "blue")
# 
# for (x in 1:length(a[[1]])) {
#   for (y in 2:(length(a[1,])-2)) {
#     if (!isTRUE(a[x,y]>=a[x,51])) {
#       a[x,y]<- NA
#     } 
#   }  
# }
# 
# y<- NULL
# Transmission_genesets2<- NULL
# library(GSEABase) # load GSEA Base package (if not already loaded)
# 
# for (x in 2:49) {
#   assign((paste0("TP",x)),a[!is.na(a[,x]),1])
#   y<- append(y, length(a[!is.na(a[,x]),1]))
#   Transmission_genesets2<- append(Transmission_genesets2,GeneSet(setName = paste0("TP",x), as.character(a[!is.na(a[,x]),1])))
#   write.table(a[!is.na(a[,x]),1], file = paste0("gene_expression_data2/TP",x), sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
# }
# 
# x<- c(1:48)
# plot(x,y,ylab = "# of genes", xlab = "Hours")
# 
# 
# GSC2<-GeneSetCollection(Transmission_genesets2)
# toGmt(GSC2, con = "gene_expression_data2/Transmission_Expression_Timepoints.gmt")

#################################
# Step 5: Creating TP genesets of peak expression from experiments that
# use RNA-seq to find gene expression every 5 hours

# read Bartfai Time Series, profiles.diff file from PlasmoDB:
d<- read.table(file = "profiles.diff.txt",header = TRUE, stringsAsFactors = TRUE, quote = "", sep = "\t")
e<- read.table(file = "profiles.diff.txt",header = TRUE, stringsAsFactors = TRUE, quote = "", sep = "\t")

# make SD,mean+1.25SD, mean+1.5SD columns
d$"SD"<- NA
e$"SD"<- NA
d$"Mean+1.25_SD"<- NA
e$"Mean+1.5_SD"<- NA

# compute the three columns
for (x in 1:length(d[[1]])) {
  d[x,10]<- sd(d[x,2:9])
  e[x,10]<- sd(e[x,2:9])
  
  d[x, 11]<- d[x, 10] * 1.25 + mean(unlist(d[x, 2:9]))
  e[x, 11]<- e[x, 10] * 1.5 + mean(unlist(e[x, 2:9]))
}

## hist(d[[10]], breaks = 100, xlab = "standard deviation", main = "Histogram of Standard Deviation", col = "blue",xlim = c(0,1000),ylim = c(0,10000))

for (x in 1:length(d[[1]])) {
  for (y in 2:(length(d[1,])-2)) {
    if (d[x,y]<d[x,11] || d[x,y]==0) {
      d[x,y]<- NA
    } 
    if (e[x,y]<e[x,11] || e[x,y]==0) {
      e[x,y]<- NA
    } 
  }  
}

y<- NULL
Bartfai_transmission_1.25SD_genesets<-NULL
Bartfai_transmission_1.5SD_genesets<-NULL
y2<- NULL
for (x in 2:9) {
  y<- append(y, length(d[!is.na(d[,x]),1]))
  Bartfai_transmission_1.25SD_genesets<- append(Bartfai_transmission_1.25SD_genesets, GeneSet(setName = paste0((x-1)*5," hours"), as.character(d[!is.na(d[,x]),1])))
  
  y2<- append(y2, length(e[!is.na(e[,x]),1]))
  Bartfai_transmission_1.5SD_genesets<- append(Bartfai_transmission_1.5SD_genesets, GeneSet(setName = paste0((x-1)*5," hours"), as.character(e[!is.na(e[,x]),1])))
  
  write.table(d[!is.na(d[,x]),1], file = paste0("gene_expression_data/Bartfai/Bartfai_1.25SD_TP",x-1), sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(e[!is.na(e[,x]),1], file = paste0("gene_expression_data/Bartfai/Bartfai_1.5SD_TP",x-1), sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
}

Bartfai_1.25SD_GSC<-GeneSetCollection(Bartfai_transmission_1.25SD_genesets)
Bartfai_1.5SD_GSC<-GeneSetCollection(Bartfai_transmission_1.5SD_genesets)

toGmt(Bartfai_1.25SD_GSC, con = "gene_expression_data/Bartfai/Bartfai_1.25SD_Transmission_Expression_Timepoints.gmt")
toGmt(Bartfai_1.5SD_GSC, con = "gene_expression_data/Bartfai/Bartfai_1.5SD_Transmission_Expression_Timepoints.gmt")

x<-c(5,10,15,20,25,30,35,40)

png(filename = "Images/Bartfai_peak_expression_sizes.png")
plot(x,y,ylab = "# of genes", type = "b", xlab = "Hours", ylim = c(150,650), main = "Number of Peak Expression Genes Over Time")
points(x,y2, col="blue", type = "b", pch = 18)
abline(h=500, lty = 2)
legend(5,650, cex = .9, bty = "o", legend = c("Mean+1.25SD", "Mean+1.5SD"), lty = c(1,1), col = c("black","blue"), title = expression(bold("Peak Expression Value")))
dev.off()

########################
########################
###Step 5b) anallyze Stunnenberg experiment

# read Stunnenberg Time Series, profiles.diff2 file:
f<- read.table(file = "profiles.diff2.txt",header = TRUE, stringsAsFactors = TRUE, quote = "", sep = "\t")
g<- read.table(file = "profiles.diff2.txt",header = TRUE, stringsAsFactors = TRUE, quote = "", sep = "\t")

# make SD,mean+1.25SD, mean+1.5SD columns
f$"SD"<- NA
g$"SD"<- NA
f$"Mean+1.25_SD"<- NA
g$"Mean+1.5_SD"<- NA

# compute the three columns
for (x in 1:length(f[[1]])) {
  f[x,10]<- sd(f[x,2:9])
  g[x,10]<- sd(g[x,2:9])
  
  f[x, 11]<- f[x, 10] * 1.25 + mean(unlist(f[x, 2:9]))
  g[x, 11]<- g[x, 10] * 1.5 + mean(unlist(g[x, 2:9]))
}

## hist(d[[10]], breaks = 100, xlab = "standard deviation", main = "Histogram of Standard Deviation", col = "blue",xlim = c(0,1000),ylim = c(0,10000))

for (x in 1:length(f[[1]])) {
  for (y in 2:(length(f[1,])-2)) {
    if (f[x,y]<f[x,11] || f[x,y]==0) {
      f[x,y]<- NA
    } 
    if (g[x,y]<g[x,11] || g[x,y]==0) {
      g[x,y]<- NA
    } 
  }  
}

y<- NULL
Stunnenberg_transmission_1.25SD_genesets<-NULL
Stunnenberg_transmission_1.5SD_genesets<-NULL
y2<- NULL
for (x in 2:9) {
  y<- append(y, length(f[!is.na(f[,x]),1]))
  Stunnenberg_transmission_1.25SD_genesets<- append(Stunnenberg_transmission_1.25SD_genesets, GeneSet(setName = paste0((x-1)*5," hours"), as.character(f[!is.na(f[,x]),1])))
  
  y2<- append(y2, length(g[!is.na(g[,x]),1]))
  Stunnenberg_transmission_1.5SD_genesets<- append(Stunnenberg_transmission_1.5SD_genesets, GeneSet(setName = paste0((x-1)*5, " hours"), as.character(g[!is.na(g[,x]),1])))
  
  write.table(f[!is.na(f[,x]),1], file = paste0("gene_expression_data/Stunnenberg/Stunnenberg_1.25SD_TP",x-1), sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(g[!is.na(g[,x]),1], file = paste0("gene_expression_data/Stunnenberg/Stunnenberg_1.5SD_TP",x-1), sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
}

Stunnenberg_1.25SD_GSC<-GeneSetCollection(Stunnenberg_transmission_1.25SD_genesets)
Stunnenberg_1.5SD_GSC<-GeneSetCollection(Stunnenberg_transmission_1.5SD_genesets)

toGmt(Stunnenberg_1.25SD_GSC, con = "gene_expression_data/Stunnenberg/Stunnenberg_1.25SD_Transmission_Expression_Timepoints.gmt")
toGmt(Stunnenberg_1.5SD_GSC, con = "gene_expression_data/Stunnenberg/Stunnenberg_1.5SD_Transmission_Expression_Timepoints.gmt")


x<-c(5,10,15,20,25,30,35,40)
png(filename = "Images/Stunnenberg_peak_expression_sizes.png")
plot(x,y,ylab = "# of genes", type = "b", xlab = "Hours", ylim = c(350,700),main = "Number of Peak Expression Genes Over Time")
points(x,y2, col="blue", type = "b", pch = 18)
abline(h=500, lty = 2)
legend("bottomright", cex=.8, legend = c("Mean + 1.25 SD", "Mean + 1.5 SD"), lty = c(1,1), col = c("black","blue"), bty = "o", title = expression(bold("Peak Expression Definition")))
dev.off()

#########
## Step 6: Create peak expression genesets for different stages of the cell
#########

# read profiles.diff3.txt taken from Su's seven stages rna-seq experiment on PlasmoDB
h<- read.table(file = "profiles.diff3.txt",header = TRUE, stringsAsFactors = TRUE, quote = "", sep = "\t")
i<- read.table(file = "profiles.diff3.txt",header = TRUE, stringsAsFactors = TRUE, quote = "", sep = "\t")

# make SD,mean+1.25SD, mean+1.5SD columns
h$"SD"<- NA
i$"SD"<- NA
h$"Mean+1.25_SD"<- NA
i$"Mean+1.5_SD"<- NA

# compute the three columns
for (x in 1:length(h[[1]])) {
  h[x,9]<- sd(h[x,2:8])
  i[x,9]<- sd(i[x,2:8])
  
  h[x, 10]<- h[x, 9] * 1.25 + mean(unlist(h[x, 2:8]))
  i[x, 10]<- i[x, 9] * 1.5 + mean(unlist(i[x, 2:8]))
}

## hist(d[[10]], breaks = 100, xlab = "standard deviation", main = "Histogram of Standard Deviation", col = "blue",xlim = c(0,1000),ylim = c(0,10000))

for (x in 1:length(h[[1]])) {
  for (y in 2:(length(h[1,])-2)) {
    if (h[x,y]<h[x,10] || h[x,y]==0) {
      h[x,y]<- NA
    } 
    if (i[x,y]<i[x,10] || i[x,y]==0) {
      i[x,y]<- NA
    } 
  }  
}

y<- NULL
Su_transmission_1.25SD_genesets<-NULL
Su_transmission_1.5SD_genesets<-NULL
y2<- NULL

r<-c("Ring", "Early Trophozite", "Late Trophozoite", "Schizont", "Gametocyte II", "Gametocyte V", "Ookinete")
for (x in 2:8) {
  q<-x-1
  y<- append(y, length(h[!is.na(h[,x]),1]))
  Su_transmission_1.25SD_genesets<- append(Su_transmission_1.25SD_genesets, GeneSet(setName = r[q], as.character(h[!is.na(h[,x]),1])))
  
  y2<- append(y2, length(i[!is.na(i[,x]),1]))
  Su_transmission_1.5SD_genesets<- append(Su_transmission_1.5SD_genesets, GeneSet(setName = r[q], as.character(i[!is.na(i[,x]),1])))
  
  write.table(h[!is.na(h[,x]),1], file = paste0("gene_expression_data/Su/Su_1.25SD_TP",x-1), sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(i[!is.na(i[,x]),1], file = paste0("gene_expression_data/Su/Su_1.5SD_TP",x-1), sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
}

Su_1.25SD_GSC<-GeneSetCollection(Su_transmission_1.25SD_genesets)
Su_1.5SD_GSC<-GeneSetCollection(Su_transmission_1.5SD_genesets)

toGmt(Su_1.25SD_GSC, con = "gene_expression_data/Su/Su_1.25SD_Transmission_Expression_Timepoints.gmt")
toGmt(Su_1.5SD_GSC, con = "gene_expression_data/Su/Su_1.5SD_Transmission_Expression_Timepoints.gmt")

x<-c(1:7)
png(filename = "Images/Su_peak_expression_sizes.png")
par(mar=c(9,5,4,4))
plot(x,y,ylab = expression(bold("# of genes")), xlab = "", type = "b", main = "Number of Peak Expression Genes During Cell Cycle", xaxt = "n")
mtext(expression(bold("Stages")),side=1,line=7)
axis(1, at=1:7, labels=c("Ring", "Early Trophozite", "Late Trophozoite", "Schizont", "Gametocyte II", "Gametocyte V", "Ookinete"), las = 2)
points(x,y2, col="blue", type = "b", pch = 18)
abline(h=500, lty = 2)
abline(h=15, lty= 2)
legend("topleft", legend = c("Mean + 1.25 SD", "Mean + 1.5 SD"), lty = c(1,1), col = c("black","blue"), bty = "o", title = expression(bold("Peak Expression Definition")))
dev.off()


####### Step 6b: Make genesets for Winzeler's Gametocye experiment

# read 3D7_Winzeler_Gametocyte.txt taken from Su's seven stages rna-seq experiment on PlasmoDB
j<- read.table(file = "3D7_Winzeler_Gametocyte.txt",header = TRUE, stringsAsFactors = TRUE, quote = "", sep = "\t")
k<- read.table(file = "3D7_Winzeler_Gametocyte.txt",header = TRUE, stringsAsFactors = TRUE, quote = "", sep = "\t")

# make SD,mean+1.25SD, mean+1.5SD columns
j$"SD"<- NA
k$"SD"<- NA
j$"Mean+1.25_SD"<- NA
k$"Mean+1.5_SD"<- NA

# compute the three columns
for (x in 1:length(j[[1]])) {
  j[x,8]<- sd(j[x,2:7])
  k[x,8]<- sd(k[x,2:7])
  
  j[x, 9]<- j[x, 8] * 1.25 + mean(unlist(j[x, 2:7]))
  k[x, 9]<- k[x, 8] * 1.5 + mean(unlist(k[x, 2:7]))
}

## hist(d[[10]], breaks = 100, xlab = "standard deviation", main = "Histogram of Standard Deviation", col = "blue",xlim = c(0,1000),ylim = c(0,10000))

for (x in 1:length(j[[1]])) {
  for (y in 2:(length(j[1,])-2)) {
    if (j[x,y]<j[x,9] || j[x,y]==0) {
      j[x,y]<- NA
    } 
    if (k[x,y]<k[x,9] || k[x,y]==0) {
      k[x,y]<- NA
    } 
  }  
}

y<- NULL
Winzeler_transmission_1.25SD_genesets<-NULL
Winzeler_transmission_1.5SD_genesets<-NULL
y2<- NULL
r<-c("1 day","2 days","3 days","6 days","8 days","12 days")

for (x in 2:7) {
  q<-x-1
  y<- append(y, length(j[!is.na(j[,x]),1]))
  Winzeler_transmission_1.25SD_genesets<- append(Winzeler_transmission_1.25SD_genesets, GeneSet(setName = r[q], as.character(j[!is.na(j[,x]),1])))
  
  y2<- append(y2, length(k[!is.na(k[,x]),1]))
  Winzeler_transmission_1.5SD_genesets<- append(Winzeler_transmission_1.5SD_genesets, GeneSet(setName = r[q], as.character(k[!is.na(k[,x]),1])))
  
  write.table(j[!is.na(j[,x]),1], file = paste0("gene_expression_data/Winzeler/Winzeler_1.25SD_TP",x-1), sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(k[!is.na(k[,x]),1], file = paste0("gene_expression_data/Winzeler/Winzeler_1.5SD_TP",x-1), sep = "\t",quote = FALSE, row.names = FALSE, col.names = FALSE)
}

Winzeler_1.25SD_GSC<-GeneSetCollection(Winzeler_transmission_1.25SD_genesets)
Winzeler_1.5SD_GSC<-GeneSetCollection(Winzeler_transmission_1.5SD_genesets)

toGmt(Winzeler_1.25SD_GSC, con = "gene_expression_data/Winzeler/Winzeler_1.25SD_Transmission_Expression_Timepoints.gmt")
toGmt(Winzeler_1.5SD_GSC, con = "gene_expression_data/Winzeler/Winzeler_1.5SD_Transmission_Expression_Timepoints.gmt")

x<-c(1,2,3,6,8,12)
png(filename = "Images/Winzeler_peak_expression_sizes.png")
plot(x,y,ylab = expression(bold("# of genes")), xlab = expression(bold("Days as a Gametocyte")), type = "b", main = "Number of Peak Expression Genes for Gametocyte", ylim = c(0,1250))
points(x,y2, col="blue", type = "b", pch = 18)
abline(h=500, lty = 2)
abline(h=15, lty= 2)
legend("topleft", legend = c("Mean + 1.25 SD", "Mean + 1.5 SD"), lty = c(1,1), col = c("black","blue"), bty = "o", title = expression(bold("Peak Expression Definition")))
dev.off()


####
## Step 7: make .gmt file of gene pathways
####

w<-read.csv(file = "Table_S3.csv", header = TRUE,stringsAsFactors = FALSE)

q<-read.table(file = "Pf3D7_gene_aliases.txt", sep = "\t", stringsAsFactors = FALSE,fill = TRUE)
for (z in 2:100) {
for (y in 2:4) {
  for (n in 1:length(q[[1]])) {
    if (q[n,y] != "") {
      w[grep(q[n,y],w[[z]],ignore.case = TRUE),z]<-unlist(q[n,1])
    }  
  }
}
}
w_2<- as.data.frame(append(w,list(C=NA),after = 1))
write.table(w_2[[1]],file = "PathwayGenesetNames",row.names = F,col.names = F)
write.table(w_2,file="Table_S3.gmt",row.names = FALSE,row.columns=F,quote = FALSE, sep = "\t")
  