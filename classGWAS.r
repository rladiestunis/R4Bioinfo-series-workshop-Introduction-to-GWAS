
install.packages(c("poolr","qqman","BGLR","rrBLUP","DT", "SNPRelate","dplyr"))
install.packages(c("rnaturalearth",'rnaturalearthdata','rgeos','ggspatial'))
library(statgenGWAS )
library(gdata)
library(rrBLUP)
library(BGLR)
library(DT)
library(SNPRelate)
library(dplyr)
library(qqman)
library(poolr)
library(OpenStreetMap)
library(rjson)
library(rgdal)
library(RgoogleMaps)
library(mapproj)
library(sf)
devtools::install_github("dkahle/ggmap", ref = "tidyup")
# load library
library(OpenStreetMap)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)





###



rm(list = ls())
setwd("C:\\Users\\tatia\\Desktop\\Training\\RiceDiversity.44K.MSU6.Genotypes_PLINK\\RiceDiversity_44K_Genotypes_PLINK\\")
Geno <- read_ped("sativas413.ped")
head(Geno)
#conatins the marker allele data in 0, 1, 2, and 3 format; 2 represents missing data
p = Geno$p;p
n = Geno$n;n
Geno = Geno$x;Geno
# Acession information
FAM <- read.table("sativas413.fam");FAM
# Map information
MAP <- read.table("sativas413.map");MAP
# Recode the data in ped file
Geno[Geno == 2] <- NA  # Converting missing data to NA
Geno[Geno == 0] <- 0  # Converting 0 data to 0
Geno[Geno == 1] <- 1  # Converting 1 to 1
Geno[Geno == 3] <- 2  # Converting 3 to 2
# Convert the marker data into matrix and transponse and check dimensions
Geno <- matrix(Geno, nrow = p, ncol = n, byrow = TRUE)
Geno <- t(Geno)
dim(Geno)


##read phenotype
rice.pheno <- read.table("sativas413.pheno", 
                         header = TRUE, stringsAsFactors = FALSE, sep = "\t")
# See first few columns and rows of the data
rice.pheno[1:5, 1:5]
dim(rice.pheno)
rownames(Geno) <- FAM$V2
table(rownames(Geno) == rice.pheno$NSFTVID)


# Now let us extract the first trait and assign it to object y
y <- matrix(rice.pheno$Flowering.time.at.Arkansas)  # # use the first trait 
rownames(y) <- rice.pheno$NSFTVID
index <- !is.na(y)
y <- y[index, 1, drop = FALSE]  
Geno <- Geno[index, ]  
table(rownames(Geno) == rownames(y))


#Here we will be using for loop to identify the missing ones and convert them with mean values across columns.
for (j in 1:ncol(Geno)) {
  Geno[, j] <- ifelse(is.na(Geno[, j]), mean(Geno[, j], na.rm = TRUE), Geno[,j])
}
#Here we will dropping the alleles with frequency < 5% and saving the object as Geno1
#Second we will match and drop the markers from map info file and save it MAP1
# Filter for minor alleles
p <- colSums(Geno)/(2 * nrow(Geno))
maf <- ifelse(p > 0.5, 1 - p, p)
maf.index <- which(maf < 0.05)
Geno1 <- Geno[, -maf.index]
dim(Geno1)
dim(Geno)


#read map file
MAP <- read.table("sativas413.map")
dim(MAP)

#subset based on retained SNPs
MAP1 <- MAP[-maf.index, ]
dim(MAP1)


##POPULATION STRUCTURE
# Create geno matrix file and assign the row and column names from fam and
# map files
Geno1 <- as.matrix(Geno1)
sample <- row.names(Geno1)
length(sample)

colnames(Geno1) <- MAP1$V2
snp.id <- colnames(Geno1)
length(snp.id)


snpgdsCreateGeno("44k.gds", genmat = Geno1, sample.id = sample, snp.id = snp.id, 
                 snp.chromosome = MAP1$V1, snp.position = MAP1$V4, snpfirstdim = FALSE)
# Now open the 44k.gds file
geno_44k <- snpgdsOpen("44k.gds")
snpgdsSummary("44k.gds")

#PCA analysis
pca <- snpgdsPCA(geno_44k, snp.id = colnames(Geno1))
#plot results of PCA
pca <- data.frame(sample.id = row.names(Geno1),
                  EV1 = pca$eigenvect[, 1],
                  EV2 = pca$eigenvect[, 2])
        
                                                                                                                                                                                    dim(pca$eigenvect)                                                                                              dim(pca$eigenvect)                                                                                              2], EV3 = pca$eigenvect[, 3], EV3 = pca$eigenvect[, 4], stringsAsFactors = FALSE)
# Plot the PCA
plot(pca$EV2, pca$EV1, xlab = "eigenvector 3", ylab = "eigenvector 4")


#add population information to the plot
# Now let us add the population information to the plot. Here we will be
# using the population information from the PCA file available online
pca_1 <- read.csv("sativas413.csv", 
                  header = TRUE, skip = 1, stringsAsFactors = FALSE)  # 431 x 12
pca_2 <- pca_1[match(pca$sample.id, pca_1$NSFTV.ID), ]

# Extract the population information and add the pca output file
pca_population <- cbind(pca_2$Sub.population, pca)
colnames(pca_population)[1] <- "population"
# Plot and add the population names
plot(pca_population$EV1, pca_population$EV2, xlab = "PC1", ylab = "PC2", pch=15, col = c(1:6)[factor(pca_population$population)])
legend(x = "bottomright", legend = levels(factor(pca_population$population)), col = c(1:6), 
       pch = 15, cex = 1.2)

##by country
# Extract the location information and add the pca output file
pca_country <- as.data.frame(cbind(as.numeric(pca$EV1),as.numeric(pca$EV2),as.numeric(pca$EV3),as.numeric(pca$EV4),as.numeric(pca_2$Latitude), as.numeric(pca_2$Longitude)));head(pca_country)
names(pca_country)=c('PC1','PC2','PC3','PC4','LAT','LONG')
pca_country=na.omit(pca_country)
#correlation with geography
cor(pca_country[,1:6],use='pairwise.complete.obs')
# load data
world <- ne_countries(scale = "medium", returnclass = "sf")

# gene world map
# extract locations
world_points<- st_centroid(world)
# extract labels
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))
#define colors for the map
NCOLS=max(unique(round(50*(pca_country$PC3+0.11))));NCOLS
COL=rainbow(NCOLS)

ggplot(data = world) +
  geom_sf() +
  labs( x = "Longitude", y = "Latitude") +
  ggtitle("World map", subtitle = paste0("(", length(unique(world$admin)), " countries)"))+
  annotation_scale(location = "bl", width_hint = 0.5) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme_bw()+theme(panel.grid.major = element_line(color = "gray60", linetype = "dashed", size = 0.25), 
                   panel.background = element_rect(fill = "aliceblue")) +
  geom_text(data= world_points,aes(x=X, y=Y, label=name),
            color = "gray20", fontface = "italic", check_overlap = T, size = 3)+
  geom_point(data = pca_country, aes(x = LONG, y = LAT), size = 4,   shape = 24, 
             fill = COL[round(50*(pca_country$PC3+0.11))])

      
      
        
        
##GWAS
#GWAS analysis will be done in rrBLUP package.
#we will prepare the genotypic file including markers, and map information
#We will run the GWAS analysis in rrBLUP package
#Correct for multiple testing
#And finally we will look for significant markers and draw the Manhattan plot.


# create the geno file for rrBLUP package GWAS analysis
geno_final <- data.frame(marker = MAP1[, 2], chrom = MAP1[, 1], pos = MAP1[, 4], t(Geno1 - 1), check.names = FALSE)  
dim(Geno1)
# create the pheno file
pheno_final <- data.frame(NSFTV_ID = rownames(y), y = y)
# Run the GWAS analysis
myGWAS <- GWAS(pheno_final, geno_final, min.MAF = 0.05, P3D = TRUE, plot=FALSE)
myGWAS1 <- GWAS(pheno_final, geno_final, n.PC=4, min.MAF = 0.05, P3D = TRUE, plot=TRUE)

myGWAS1=myGWAS1[order(myGWAS1$y),]
head(myGWAS1)

##PLOT MANHATTAN PLOT
manhattan(x = myGWAS, chr = "chrom", bp = "pos", p = "y", snp = "marker", col = c("blue4", 
                                                                                  "orange3"), suggestiveline = -log10(1e-04), logp = TRUE)
myGWAS=myGWAS[order(myGWAS$y),]
head(myGWAS)

####GPS
#library(ape)
#library(phytools)
ADM=read.csv("Admixture.csv",row.names = 1)
head(ADM)

# Instruction to name bars:
barNaming <- function(vec) {
  retVec <- vec
  for (k in 2:length(vec)) {
    if (vec[k - 1] == vec[k])
      retVec[k] <- ""
  }
  return(retVec)
}


# Size of the plot - these are margins, also related to the end appearance of the plot
par(mar = c(1, 1, 1, 1))

#This is my selection of colors and display for K=8 
x=barplot(t(as.matrix(ADM[, 7:14])),
        col = rainbow(8),
          horiz = FALSE,
          las = 2,
        cex.axis = 0.1,
        cex.names = 0.65,
)
text(x,0.5,barNaming(ADM$Country),srt=45, col='white')

