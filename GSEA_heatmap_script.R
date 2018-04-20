# GSEA analysis via heatmap
## For each experiment w/ a .gmt file put through GSEA, this script will
## create a heatmap of the enrichment values for each geneset at each phenotype label tested



# GSEApath is a string containing the filepath to the directory that containins all of the GSEA output files.
## ex: "gsea_home/output/jul27"

# genesetlist is a character vector IDENTICAL TO the names of each geneset in the geneset collection run through GSEA
## GSEA normally puts everything in capital letters
## ex: c("TP1","TP2","TP3")


# heatmap_path is the path to the directory that the heatmap should be saved under

# heatmap_name is the name of the heatmap file (leave out extensions)

# filenames is a character vector of all the filenames of the GSEA results folders in the GSEApath directory

### IMPORTANT: ALL FILE NAMES MUST HAVE THE FOLLOWING FORMAT: ###
# geneset.phenotype1.phenotype2.Gsea.#
## where phenotype1 & phenotype2 are the two phenotypes tested in the file. 
## geneset can be anything, used to differentiate between genesets
### Phenotype labels must be in the correct order (if the phenotype set tested is "abc" vs. "def" then phenotype 1 must be "abc")
## Gsea.# (# is a random number) should be automatically added to the end of the name after the GSEA program runs
## ex: c("Derisi.106.106CQ.Gsea.1469642172842","Derisi.106.176I.Gsea.1469642286768", ....)

## make sure to install plyr & RColorBrewer & gplots if not already installed
# install.packages("plyr")
# install.packages("ggplots")
# install.packages("RColorBrewer")


detach("GSEABase", unload=TRUE) #may cause plyr package to not work
library(plyr)
library(gplots)
library(RColorBrewer)



GSEA_heatmap<- function(GSEAhome, filenames, genesetlist,heatmap_path, heatmap_name) {
  b<-as.data.frame(genesetlist)
  names(b)<-"NAME"
  l<-b
  for (n in filenames) {
    a1<-read.table(file = paste0(GSEAhome,"/",n,"/gsea_report_for_",strsplit(x=n, "[.]")[[1]][2],"_", strsplit(x = n,"[.]")[[1]][5],".xls"),sep = "\t",header = TRUE, stringsAsFactors = T)
    a2<-read.table(file = paste0(GSEAhome,"/",n,"/gsea_report_for_",strsplit(x=n, "[.]")[[1]][3],"_", strsplit(x = n,"[.]")[[1]][5],".xls"),sep = "\t",header = TRUE, stringsAsFactors = T)
    a<- rbind(a1,a2)
    a<-a[order(a[[1]]),]
    a<-join(b,a,"NAME")
    for (x in 1:length(a[[1]])) {
      if (a[x,8] >= 0.05 || is.na(a[x,8])) {
        a[x,6] <- 0
      }
    }
    a[order(a[[1]]),]
    l[order(l[[1]]),]
    l<-cbind(l,a[, 6])
    colnames(l)[names(l)== "a[, 6]"] <- paste0(strsplit(n,"[.]")[[1]][2]," vs. ",strsplit(n,"[.]")[[1]][3])
  }
  rownames(l) <- l[,1]
  l<-l[,2:ncol(l)]
  lm<-as.matrix(l)
  crp<-colorRampPalette(brewer.pal(11,"PiYG"))
  png(filename = paste0(heatmap_path,"/",heatmap_name,".png"),res = 300,width = 9*300, height = 9*300)
  heatmap.2(lm, Colv = "NA", Rowv = "NA", col = crp, dendrogram = "none", key.title = "Key", density.info = "none",trace = "none", margins = c(12,14),revC = T,main = heatmap_name,cexCol = 1.5)
  dev.off()
  return(paste0("Heatmap saved under ",heatmap_path,"/",heatmap_name,".png"))
}
for (n in 2:4) {
  o<-l[11,n]
  p<-l[12,n]
  r<- (o+p)/2
}

#variables for each geneset
#Su
genesetlist<- c("RING", "EARLY TROPHOZITE", "LATE TROPHOZITE", "SCHIZONT", "GAMETOCYTE II", "GAMETOCYTE V", "OOKINETE")
GSEAhome<- "gsea_home/output/aug01"
filenames <- c("Su.352K.176I.Gsea.1470107601473",
               "Su.352K.106.Gsea.1470107554025",
               "Su.176I.106.Gsea.1470107512079",
               "Su.352K.352KCQ.Gsea.1470107465438",
               "Su.176I.176ICQ.Gsea.1470107418962",
               "Su.106.106CQ.Gsea.1470107359677")
filenames<- c("Su.treated.control.Gsea.1470163559049","Su.control.treated.Gsea.1470163543863")
heatmap_name <- "Su"
GSEA_heatmap(GSEAhome,filenames,genesetlist,heatmap_path,heatmap_name)

#Winzeler
genesetlist<- c("1 DAY", "2 DAYS", "3 DAYS", "6 DAYS", "12 DAYS")
filenames<- c("Winzeler.352K.176I.Gsea.1470150206495",
              "Winzeler.352K.106.Gsea.1470150178888",
              "Winzeler.176I.106.Gsea.1470150159290",
              "Winzeler.352K.352KCQ.Gsea.1470150123369",
              "Winzeler.176I.176ICQ.Gsea.1470150099291",
              "Winzeler.106.106CQ.Gsea.1470150064061")
filenames<-c("Winzeler.treated.control.Gsea.1470163471396","Winzeler.control.treated.Gsea.1470163450691")
GSEAhome<- "gsea_home/output/aug02"
heatmap_name <- "Winzeler 1"
GSEA_heatmap(GSEAhome,filenames,genesetlist,heatmap_path,heatmap_name)

#Bartfai
genesetlist<-NULL
for (n in 1:8) {
  genesetlist<-append(genesetlist,paste0(n*5," HOURS"))
}
GSEAhome<- "gsea_home/output/aug02"
filenames<- c("Bartfai.352K.176I.Gsea.1470151293247",
              "Bartfai.352K.106.Gsea.1470151272029",
              "Bartfai.176I.106.Gsea.1470151241819",
              "Bartfai.352K.352KCQ.Gsea.1470151203155",
              "Bartfai.176I.176ICQ.Gsea.1470151145651",
              "Bartfai.106.106CQ.Gsea.1470151126935")
filenames <- c("Bartfai.control.treated.Gsea.1470163348639","Bartfai.treated.control.Gsea.1470163362030")
heatmap_name <- "Bartfai"
GSEA_heatmap(GSEAhome,filenames,genesetlist,heatmap_path,heatmap_name)

#Derisi
genesetlist<-NULL
for (n in 1:53) {
  genesetlist<-append(genesetlist,paste0("TP ",n))
}
GSEAhome<- "gsea_home/output/aug02"
GSEAhome<- "gsea_home/output/nov07"

heatmap_name <- "Derisi 1.1"
filenames<- c("Derisi.176I.106.Gsea.1470151757902",
              "Derisi.352K.106.Gsea.1470151776876",
              "Derisi.352K.176I.Gsea.1470151789168",
              "Derisi.352K.352KCQ.Gsea.1470151740610",
              "Derisi.176I.176ICQ.Gsea.1470151721968",
              "Derisi.106.106CQ.Gsea.1470151700924")
filenames<- c("Derisi.352K.352KCQ.Gsea.1478577736851",
              "Derisi.352KCQ.352K.Gsea.1478577757461")
filenames<- c("Derisi.control.treated.Gsea.1470163171163","Derisi.treated.control.Gsea.1470163186828")
GSEA_heatmap(GSEAhome,filenames,genesetlist,heatmap_path,heatmap_name)

#stunenberg
GSEAhome<- "gsea_home/output/aug02"
genesetlist<-NULL
for (n in 1:8) {
  genesetlist<-append(genesetlist,paste0(n*5," HOURS"))
}
heatmap_name <- "Stunnenberg"
filenames<- c("Stunnenberg.352K.176I.Gsea.1470152129639",
              "Stunnenberg.352K.106.Gsea.1470152118314",
              "Stunnenberg.176I.106.Gsea.1470152106126",
              "Stunnenberg.352K.352KCQ.Gsea.1470152090848",
              "Stunnenberg.176I.176ICQ.Gsea.1470152072566",
              "Stunnenberg.106.106CQ.Gsea.1470152054255")
filenames<- c("Stunnenberg.control.treated.Gsea.1470162489505","Stunnenberg.treated.control.Gsea.1470162668133")
GSEA_heatmap(GSEAhome,filenames,genesetlist,heatmap_path,heatmap_name)

#location
GSEAhome<-"gsea_home/output/aug02"
heatmap_path<- "GSEA_Results_2"
genesetlist <- c("APICAL","APICOPLAST","CYTOPLASM","ER_GOLGI","EXPORTED","FOOD_VACUOLE","IMC","MITOCHONDRION","NUCLEUS","PPM","PV")
filenames <- c("Location.352K.176I.Gsea.1470167440103",
                 "Location.352K.106.Gsea.1470167431303",
                 "Location.176I.106.Gsea.1470167414220",
                 "Location.352K.352KCQ.Gsea.1470167399633",
                 "Location.176I.176ICQ.Gsea.1470167384731",
                 "Location.106.106CQ.Gsea.1470167367271")
heatmap_name<- "Location"
filenames<-c("location.treated.control.Gsea.1470161922178","location.control.treated.Gsea.1470161881830")
GSEA_heatmap(GSEAhome,filenames,genesetlist,heatmap_path,heatmap_name)

#pathway
filenames <- c("Pathway.352K.106.Gsea.1473010790640",
               "Pathway.352K.176I.Gsea.1473010775491",
               "Pathway.176I.106.Gsea.1473010759152",
               "Pathway.352K.352KCQ.Gsea.1473010741287",
               "Pathway.106.106CQ.Gsea.1473010717336",
               "Pathway.176I.176ICQ.Gsea.1473009981195")
heatmap_name<-"Pathways"
heatmap_path<- "GSEA_Results"
GSEAhome<-"gsea_home/output/sep04"
w<-read.table("PathwayGenesetNames",stringsAsFactors = F)
w<-data.frame(lapply(w, function(v) {
  if (is.character(v)) return(toupper(v))
  else return(v)
}))
genesetlist<-w
GSEA_heatmap(GSEAhome,filenames,genesetlist,heatmap_path,heatmap_name)


