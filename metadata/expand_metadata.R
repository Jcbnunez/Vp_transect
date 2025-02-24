library(tidyverse)
library(foreach)
library(data.table)

dat <- fread("https://raw.githubusercontent.com/Jcbnunez/Vp_transect/refs/heads/main/metadata/met_dat_fG.csv")

##my_exp_meta =
##  foreach(i=1:dim(dat)[1],
##          .combine = "rbind")%do%{
##            
##            
##            tmp <- dat[i,]
##            
##            data.frame(tmp, rep=1:4) -> out
##            
##            return(out)
##            
##          }
##

####
dat2 <- fread("https://raw.githubusercontent.com/Jcbnunez/Vp_transect/refs/heads/main/metadata/file.names.txt", header = FALSE)

dat2[1:4,]$V1  = gsub("13-","", dat2[1:4,]$V1)

names(dat2) = "file_name"
dat2 %>%
  separate(file_name, into = c("biotag","techtag"), 
           sep = "_", remove = F ) %>%
  separate(biotag, into = c("Site", "ind"), sep = "-") ->
  dat2_exp
  
full_join(dat, dat2_exp) %>% View


write.table(my_exp_meta, 
            file = "meta.clean.exp.txt", 
            append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE, qmethod = c("escape", "double"),
            fileEncoding = "")