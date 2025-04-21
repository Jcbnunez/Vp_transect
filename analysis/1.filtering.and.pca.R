library(tidyverse)
library(magrittr)
library(data.table)
library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(adegenet)
library(reshape2)
library(FactoMineR)
require(gtools)
require(foreach)
library("missMDA")

#library(BEDASSLE)

#####
samps <- fread("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/Vp_transect/metadata/metadata.final.Feb24.2025.txt")

#####
genofile <- seqOpen("/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/Vp_transect/data/VeroDDRAD.gds")
seqResetFilter(genofile)


#### Load the filtering SNP object -- which JCBN created
snps.dt <- data.table(chr=seqGetData(genofile, "chromosome"),
                      pos=seqGetData(genofile, "position"),
                      variant.id=seqGetData(genofile, "variant.id"),
                      nAlleles=seqNumAllele(genofile),
                      missing=seqMissing(genofile))

snpgdsGetGeno(genofile, 
              verbose=TRUE, with.id=TRUE) ->
  GENO.matrix
#GENO.matrix$genotype[1:10,1:10]

###
#dp <- seqGetData(genofile, "annotation/format/DP")
GENO.matrix$genotype %>%
  as.data.frame() ->
  GENO.matrix.dt

colnames(GENO.matrix.dt) <- paste(seqGetData(genofile, "chromosome"), seqGetData(genofile, "position") , seqGetData(genofile, "variant.id"), sep="_")
rownames(GENO.matrix.dt) <- seqGetData(genofile, "sample.id")
dim(GENO.matrix.dt)[2] -> totSNPs
dim(GENO.matrix.dt)[1] -> totInds

NAinds <- apply(GENO.matrix.dt, 1, function(x) sum(is.na(x)))/totSNPs
data.frame(
  ind = names(NAinds),
  perc = NAinds) ->
  NAinds

NAinds %>%
  ggplot(aes(
    x=ind,
    y=perc
  )) +
  geom_bar(stat = "identity")+
  coord_flip()

median(NAinds$perc)
NAinds %>%
  filter(perc < 0.75) %>%
  .$ind -> keep_ind


### recalibrate after ind filter
data.frame(
  snp = colnames(GENO.matrix.dt),
  perc = apply(GENO.matrix.dt, 2, function(x) sum(is.na(x)))) ->
  NAsnps

hist(NAsnps$perc)

NAsnps %>%
  filter(perc < 20) %>%
  .$snp -> keep_snp
length(keep_snp)/totSNPs

#### final filter
###              
GENO.matrix.dt[
  which(rownames(GENO.matrix.dt) %in% keep_ind),
  which(colnames(GENO.matrix.dt) %in% keep_snp)
] -> GENO.matrix.dt.flt
GENO.matrix.dt.flt[1:10,1:10]
dim(GENO.matrix.dt.flt)
####
data.frame(
  snp = colnames(GENO.matrix.dt),
  perc = apply(GENO.matrix.dt, 2, function(x) sum(x/2, na.rm = T))) ->
  NA_AF


NA_AF %>%
  filter(perc/96 > 0.1) %>%
  .$snp -> keep_snp_MAF

keep_snp_final = keep_snp[keep_snp %in% keep_snp_MAF]

GENO.matrix.dt[
  which(rownames(GENO.matrix.dt) %in% keep_ind),
  which(colnames(GENO.matrix.dt) %in% keep_snp_final)
] -> GENO.matrix.dt.flt
GENO.matrix.dt.flt[1:10,1:10]
dim(GENO.matrix.dt.flt)


####
imputePCA(GENO.matrix.dt.flt, scale = FALSE) %>%
  .$completeObs ->
  GENO.matrix.dt.flt.imp
dim(GENO.matrix.dt.flt.imp)

save(GENO.matrix.dt.flt.imp, file = "/Users/jcnunez/Library/CloudStorage/OneDrive-UniversityofVermont/Documents/GitHub/Vp_transect/data/GENO.matrix.dt.flt.imp.Rdata")

###
samps %>%
  mutate(sampleid = paste(Site, ind, sep = "_")) %>%
  filter(sampleid %in% rownames(GENO.matrix.dt.flt)) ->
  samps.flt

####
groups <- separate(data.frame(ind = row.names(GENO.matrix.dt.flt.imp)), ind, remove = F, into = c("pop","sample"), sep = "_")
dapc(GENO.matrix.dt.flt.imp, groups$pop, 
     n.pca = 30, n.da=3 ) -> dapc.cal
optim.a.score(dapc.cal)
dapc(GENO.matrix.dt.flt.imp, groups$pop, 
     n.pca = optim.a.score(dapc.cal)$best, 
     n.da=3 ) -> dapc.cal.optim
scatter(dapc.cal.optim)
#paste(samps.flt$Transect_Type, samps.flt$Transect_step)
dapc.cal.optim$ind.coord %>%
  as.data.frame() %>%
  mutate(sampleid = rownames(.)) %>%
  left_join(samps.flt) %>%
  ggplot(aes(
    x=LD1,
    y=LD2,
    col=Transect_Type,
    label=Transect_step
    #col=Transect_step
  )) + #geom_point(size = 2.3) +
  geom_text() 

dapc.cal.optim$ind.coord %>%
  as.data.frame() %>%
  mutate(sampleid = rownames(.)) %>%
  left_join(samps.flt) %>%
  ggplot(aes(
    x=LD1,
    y=Transect_step,
    color=Transect_Type
  )) + geom_point() +
  geom_smooth(method = "lm")

dapc.cal$ind.coord %>%
  as.data.frame() %>%
  mutate(sampleid = rownames(.)) %>%
  left_join(samps.flt) %>%
  ggplot(aes(
    x=LD2,
    y=Transect_step,
    color=Transect_Type
  )) + geom_point() +
  geom_smooth(method = "lm")

dapc.cal$ind.coord %>%
  as.data.frame() %>%
  mutate(sampleid = rownames(.)) %>%
  left_join(samps.flt) %>%
  ggplot(aes(
    x=LD3,
    y=Transect_step,
    color=Transect_Type
  )) + geom_point() +
  geom_smooth(method = "lm")


####  
#### #### #### #### #### #### 
GENO.matrix.dt.flt.imp %>% 
  PCA(graph = F, ncp = 5) ->
  pca.object.flt
#save(pca.object.05, file = "Vero.pca.object.05.Rdata")

pca.object.flt$ind$coord %>%
  as.data.frame() %>% 
  mutate(sampleid = rownames(.)) %>%
  left_join(samps.flt) ->
  pca.meta.dim.flt
  
#### PLOT PCA
pca.meta.dim.flt %>%
separate(sampleId, remove = F, into = c("pop","sample"), sep = "_") %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    label = sampleId,
    color = pop,
  )) +
  geom_point(size = 2.0) + 
  geom_text(size = 1.7)


pca.meta.dim.flt %>%
  ggplot(aes(
    x=Dim.1,
    y=Dim.2,
    #label = sampleId,
    color = Transect_Type,
  )) +
  geom_point(size = 2.0) 



pca.meta.dim.flt %>%
  filter(sampleid != "OG_F2") %>%
  ggplot(aes(
    x=Transect_step,
    y=Dim.1,
    color = Transect_Type,
  )) +
  geom_point(size = 2.0) +
  geom_smooth(method = "lm")

pca.meta.dim.flt %>%
  filter(sampleid != "OG_F2") %>%
  ggplot(aes(
    x=Transect_step,
    y=Dim.2,
    color = Transect_Type,
  )) +
  geom_point(size = 2.0) +
  geom_smooth(method = "lm")
