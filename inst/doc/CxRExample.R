## 2019 wgaim 2.0.x script

############ Cascades x RAC875-2 Zinc Experiment
############ Example 1:

data("phenoCxR", package = "wgaim")

## scatter plot of shoot vs zinc concentration

ggplot(phenoCxR, aes(y = znconc, x = shoot)) + geom_point(colour = "blue", shape = 16) +
    xlab("Shoot Length") + ylab("Zinc Concentration") + theme_bw()

## create means of shoot and zinc concentration

zincm <- with(phenoCxR, aggregate(cbind(znconc, shoot), list(Block = Block, Type= Type), mean))
zincd <- cbind(val = c(zincm$znconc, zincm$shoot), trait = rep(c("Zn Conc", "Shoot"), each = 6), rbind(zincm[,1:2]), row.names = NULL)

## bar ploot of means

ggplot(zincd, aes(y = val, x = Type, fill = Block)) + facet_wrap(~ trait, ncol = 2, scales = "free") +
    geom_bar(stat = "identity", position = "dodge", colour = "grey50", width = 0.6) + coord_flip() +
    xlab("") + ylab("Mean trait values")

## base model

sh.fm <- asreml(shoot ~ Type, random = ~ Block + id, data = phenoCxR)
summary(sh.fm)$varcomp
wald(sh.fm)

###### Genetic data

## read in genetic data to an R/qtl "cross" object

data("genoCxR", package = "wgaim")
wgpath <- system.file("extdata", package = "wgaim")
read.csv(paste(wgpath, "/genoCxR.csv", sep = ""), header = FALSE)[1:10,1:10]
genoCxR <- read.cross("csvr", file="genoCxR.csv", genotypes=c("AA","BB"),
               dir = wgpath, na.strings = c("-", "NA"))

summary(genoCxR)
class(genoCxR)
names(genoCxR$geno)
names(genoCxR$geno$"1B")
genoCxR$geno$"1B"$data[1:8,1:8]
genoCxR$geno$"1B"$map

## convert "cross" object to a "interval" object for use in wgaim

genoCxRi <- cross2int(genoCxR, consensus.mark=TRUE, id = "id", impute = "MartinezCurnow")
class(genoCxRi)
names(genoCxRi$geno$"1B")
genoCxRi$geno$"1B"$imputed.data[1:8,1:8]
genoCxRi$geno$"1B"$interval.data[1:6,1:6]

## wgaim QTL analyses

## QTL interval analysis

sh.qtlI <- wgaim(sh.fm, genoCxRi, merge.by = "id", gen.type = "interval",
                 method = "fixed", selection = "interval", na.action = na.method(x = "include"))

summary(sh.qtlI, genoCxRi, LOD = FALSE)

## QTL marker analysis

sh.qtlM <- wgaim(sh.fm, genoCxRi, merge.by = "id", gen.type = "marker",
                 method = "fixed", selection = "interval", na.action = na.method(x = "include"))

summary(sh.qtlM, genoCxRi)

########## end script
##########
