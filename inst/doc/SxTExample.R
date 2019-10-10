## 2019 wgaim 2.0.x script

########### Sunco-Tasman Example
########### Example 2:

## import phenotypic data

wgpath <- system.file("extdata", package = "wgaim")

####### initial model and diagnostics

data("phenoSxT", package = "wgaim")
phenoSxT <- asreml.read.table(paste(wgpath,"\\phenoSxT.csv", sep =""), header = TRUE, sep = ",")

st.fmI <- asreml(myield ~ Type, random = ~ id + Rep + Range:Row + Millday, residual = ~ Millday:ar1(Millord), data = phenoSxT)

phenoSxTd <- cbind.data.frame(phenoSxT, Residuals = resid(st.fmI))

ggplot(phenoSxTd, aes(y = Residuals, x = as.numeric(Millord))) + facet_wrap(~ Millday) +
    geom_hline(yintercept = 0, linetype = 2) + geom_point(shape = 16, colour = "blue") +
    xlab("Mill order") + theme_bw()

field.resid <- coef(st.fmI, pattern = "Range:Row")
rrd <- data.frame(field.resid = field.resid, Range = factor(rep(1:12, each = 31)), Row = rep(1:31, 12))

ggplot(rrd, aes(y = field.resid, x = Row)) + facet_wrap(~ Range) +
    geom_hline(yintercept = 0, linetype = 2) + geom_point(shape = 16, colour = "blue") +
    geom_line(colour = "blue") + xlab("Row") + ylab("Field Residuals") + theme_bw()

summary(st.fmI)$varcomp

####### full model

st.fmF <- asreml(myield ~ Type + lord + lrow, random = ~ id + Rep + Range:Row + Millday,
                 residual = ~ Millday:ar1(Millord), data = phenoSxT,
                 na.action = na.method(x = "include"))
summary(st.fmF)$varcomp

# fit null model (for comparison)

st.fmN <- asreml(myield ~ 1, random = ~ id, data = phenoSxT, na.action = na.method(x = "include"))

####### read in genetic data and check

genoSxT <- read.cross("csvr", file="genoSxT.csv", genotypes=c("AA","BB"),
                    dir = wgpath, na.strings = c("-", "NA"))

nmar(genoSxT)
nt <- ntyped(genoSxT, "ind")
nt[nt < 120]
genoSxT <- subset(genoSxT, ind = 1:180)
genoSxT <- cross2int(genoSxT, consensus.mark = TRUE, id = "id", impute = "MartinezCurnow")
nmar(genoSxT)

# plot the linkage map in various ways

linkMap(genoSxT, names(nmar(genoSxT)), m.cex = 0.5)
linkMap(genoSxT, names(nmar(genoSxT)), m.cex = 0.5,
         chr.dist = list(start = 25, end = 180), marker.names = "dist")

####### wgaim QTL analyses for Full and Null model and diagnostics

st.qtlN <- wgaim(st.fmN, genoSxT, merge.by = "id", gen.type = "interval", method = "fixed",
    selection = "interval", trace = "nullmodel.txt")

st.qtlF <- wgaim(st.fmF, genoSxT, merge.by = "id", gen.type = "interval", method = "fixed",
    selection = "interval", trace = "nullmodel.txt")

# diagnostics

# plot the outlier statistics

outStat(st.qtlF, genoSxT, iter = 1:2, statistic = "outlier")
outStat(st.qtlF, genoSxT, iter = 1:2, statistic = "blups")
outStat(st.qtlF, genoSxT, iter = 1:5, chr = c("2B","4B","6B","7D"))

# trace of forward selection proces

tr(st.qtlF, iter = 1:10, digits = 3)

####### visualisation and summary

# summary and table

summary(st.qtlF, genoSxT, LOD = TRUE)
summary(st.qtlN, genoSxT, LOD = TRUE)
qtlTable(st.qtlF, st.qtlN, intervalObj = genoSxT, labels = c("Full", "Null"), columns = 1:8)

# plot qtl on linkage map in various ways

linkMap(st.qtlF, genoSxT, marker.names = "dist", trait.labels = "Full")
linkMap(st.qtlF, genoSxT, chr = c("1B", "2B", "4D"), marker.names = "dist", trait.labels = "Full")

# with both models

linkMap.default(list(st.qtlF, st.qtlN), genoSxT, marker.names = "dist",
                 trait.labels = c("Full", "Null"), list.cex = list(m.cex = 0.7))
linkMap.default(list(st.qtlF, st.qtlN), genoSxT, marker.names = "dist",
                 trait.labels = c("Full", "Null"))
linkMap.default(list(st.qtlF, st.qtlN), genoSxT, chr = c("1B", "2B", "4D"),
                 marker.names = "dist", trait.labels = c("Full", "Null"))

# customize your qtl plots

linkMap.default(list(st.qtlF, st.qtlN), genoSxT, marker.names = "dist",
    trait.labels = c("Full", "Null"), list.col = list(q.col = c("skyblue3",
    "salmon2"), m.col = "red", t.col = c("skyblue3", "salmon2")))

linkMap.default(list(st.qtlF, st.qtlN), genoSxT, marker.names = "dist",
    trait.labels = c("Full", "Null"), list.col = list(q.col = rep(gray(0.8), 2),
    m.col = "black", t.col = "black"), list.cex = list(t.cex = 0.8), col = "gray")

####### whole genome marker analysis and summary

st.qtlFM <- wgaim(st.fmF, genoSxT, merge.by = "id", gen.type = "marker", method = "fixed",
                 selection = "interval", trace = "fullmodel.txt")

summary(st.qtlFM, genoSxT, LOD = TRUE)

outStat(st.qtlFM, genoSxT, iter = 1:5, statistic = "blups")
outStat(st.qtlFM, genoSxT, iter = 1:5, chr = c("2B","4B","5A","6B","7D"))

linkMap(st.qtlFM, genoSxT, chr.dist = list(start = 15), trait.labels = "Full", cex = 2.0, pch = 16)
linkMap(st.qtlFM, genoSxT, trait.labels = "Full", cex = 2.0, pch = 16)

# Null model

st.qtlNM <- wgaim(st.fmN, genoSxT, merge.by = "id", gen.type = "marker", method = "fixed",
                 selection = "interval", trace = "nullmodel.txt")

# linkage map with QTL

linkMap.default(list(st.qtlFM, st.qtlNM), genoSxT, marker.names = "dist",
    trait.labels = c("Full", "Null"), list.col = list(q.col = c("red",
    "light blue"), m.col = "red", t.col = c("red", "light blue")),
    list.cex = list(t.cex = 0.9, m.cex = 0.7), col = "black", cex = 2, pch = 16)

####### Random formulation

outStat(st.qtlF, genoSxT, iter = c(1,2,5), chr = c("1B","2B","6B"), statistic = "blups")

st.qtlFR <- wgaim(st.fmF, genoSxT, merge.by = "id", gen.type = "interval", method = "random",
                 selection = "interval", trace = "fullmodel.txt")

summary(st.qtlFR, genoSxT)

####### end script
