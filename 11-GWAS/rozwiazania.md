1. Wczytaj i załaduj poniższe pakiety:
   - `rrBLUP`
   - `BGLR`
   - `DT`
   - `SNPRelate`
   - `dplyr`
   - `qqman`
   - `poolr`
   - `qqman`

```{r}
packages <- c("rrBLUP"
   , "BGLR"
   , "DT"
   , "SNPRelate"
   , "dplyr"
   , "qqman"
   , "poolr")

{for (pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    }
  }
}

library(pkg, character.only = TRUE)
```

2. Wczytaj dane genotypowe z pliku `.ped`, informacje o osobnikach z pliku `.fam`, informacje o mapowaniu z pliku `.map`.

```{r}
Geno <- read_ped("sativas413.ped")
# wczytujemy kolumny jako osobne wartości
p = Geno$p
n = Geno$n
Geno = Geno$x
head(Geno)
Geno

FAM <- read.table("sativas413.fam")
head(FAM)

MAP <- read.table("sativas413.map")
head(MAP)
```
- Przekoduj wartości markerów zgodnie z poniższym schematem:
  - **2** → **NA** (brakujące dane).  
  - **0** → **0** (homozygota dla allelu głównego).  
  - **1** → **1** (heterozygota).
  - **3** → **2** (homozygota dla allelu mniejszościowego).  
  
```{r}
Geno[Geno == 2] <- NA
Geno[Geno == 0] <- 0
Geno[Geno == 1] <- 1
Geno[Geno == 3] <- 2
```
Przekonwertuj dane na macierz i dokonaj jej transpozycji.

```{r}
Geno <- matrix(Geno, nrow = p, ncol = n, byrow = TRUE)
Geno <- t(Geno)
```
   
Podaj wymiary macierzy (liczbę osobników i markerów SNP).

```{r}
dim(Geno)
```

3. Wczytaj dane fenotypowe i sprawdź ich zgodność z danymi genotypowymi.
```{r}
rice.pheno <- read.table("RiceDiversity_44K_Phenotypes_34traits_PLINK.txt",
                         header = TRUE, stringsAsFactors = FALSE, sep = "\t")
head(rice.pheno)

# wymiary - pierwsza wartość powinna być taka sama jak w `dim(Geno)`
dim(rice.pheno)

# przypisujemy nazwy wierszy dla macierzy Geno na podstawie drugiej kolumny (V2) z ramki FAM, zawierającej identyfikatory próbek

rownames(Geno) <- FAM$V2

# sprawdzenie zgodności
table(rownames(Geno) == rice.pheno$NSFTVID)
```

Wyodrębnij pierwszą cechę.

```{r}
y <- matrix(rice.pheno$Flowering.time.at.Arkansas)
rownames(y) <- rice.pheno$NSFTVID
index <- !is.na(y)
y <- y[index, 1, drop = FALSE]
Geno <- Geno[index, ]
table(rownames(Geno) == rownames(y))
```

4. Przeprowadź kontrolę jakości (QC) danych markerowych:
Zastąp brakujące dane markerowe średnią wartością dla każdego markera.

```{r}
for (j in 1:ncol(Geno)){
  Geno[, j] <- ifelse(is.na(Geno[, j]), mean(Geno[, j], nar.rm = TRUE), Geno[, j])
}
```

Odfiltruj markery z MAF < 5%.
```{r}
# obliczanie frekwencji allelu mniejszościowego dla każdego SNP
p <- colSums(Geno)/(2 * nrow(Geno))

# definiujemy MAF
maf <- ifelse(p > 0.5, 1-p, p)
maf.index <- which(maf < 0.05)
Geno1 <- Geno[, -maf.index]

# sprawdzamy wymiary nowej macierzy
dim(Geno1)
```

Zaktualizuj plik `.map` i podaj nowe wymiary danych genotypowych oraz informacji o markerach.

```{r}
MAP <- read.table("sativas413.map")
dim(MAP)
MAP1 <- MAP[-maf.index, ]
dim(MAP1)
```

5. Wykonaj analizę PCA:
Utwórz macierz markerów.

```{r}
Geno1 <- as.matrix(Geno1)
sample <- row.names(Geno1)
length(sample)

colnames(Geno1) <- MAP1$V2
snp.id <- colnames(Geno1)
length(snp.id)
```

Utwórz plik GDS.

```{r}
snpgdsCreateGeno("44k.gds", genmat = Geno1, sample.id = sample, snp.id = snp.id, 
                 snp.chromosome = MAP1$V1, snp.position = MAP1$V4, snpfirstdim = FALSE)

geno_44k <- snpgdsOpen("44k.gds")
snpgdsSummary("44k.gds")
```

Przeprowadź PCA.

```{r}
pca <- snpgdsPCA(geno_44k, snp.id = colnames(Geno1))
pca <- data.frame(sample.id = row.names(Geno1), 
                  EV1 = pca$eigenvect[, 1], 
                  EV2 = pca$eigenvect[, 2], 
                  EV3 = pca$eigenvect[, 3], 
                  EV4 = pca$eigenvect[, 4], 
                  stringsAsFactors = FALSE)

plot(pca$EV2, pca$EV1, xlab = "PC2", ylab = "PC1")
```

Wczytaj dodatkowe informacje o próbkach z pliku `gerplasm.csv`.

```{r}
pca_1 <- read.csv("RiceDiversity.44K.germplasm.csv", 
                  header = TRUE, skip = 1, stringsAsFactors = FALSE)
pca_2 <- pca_1[match(pca$sample.id, pca_1$NSFTV.ID), ]

pca_population <- cbind(pca_2$Sub.population, pca)
colnames(pca_population)[1] <- "population"

plot(pca_population$EV1, pca_population$EV2, xlab = "PC1", ylab = "PC2", 
     col = c(1:6)[factor(pca_population$population)])
legend(x = "topright", legend = levels(factor(pca_population$population)), 
       col = c(1:6), pch = 1, cex = 0.6)
```

6. Przygotuj dane do analizy GWAS:
Przygotuj dane genotypowe i fenotypowe do analizy.

```{r}
geno_final <- data.frame(marker = MAP1[, 2], chrom = MAP1[, 1], pos = MAP1[, 4], 
                         t(Geno1 - 1), check.names = FALSE)

pheno_final <- data.frame(NSFTV_ID = rownames(y), y = y)
```

Wykonaj analizę GWAS.
```{r}
GWAS <- GWAS(pheno_final, geno_final, min.MAF = 0.05, P3D = TRUE, plot = FALSE)
```


7. Wyodrębnij istotne markery SNP.

```{r}
GWAS_1 <- GWAS %>% filter(y != "0")
GWAS_1 %>% filter(y < 1e-04)
```

Podaj listę markerów SNP spełniających ustalone kryterium p-wartości.

```{r}
head(GWAS_1)
```

8. Stwórz wykres Manhattan.

```{r}
manhattan(x = GWAS_1, chr = "chrom", bp = "pos", p = "y", snp = "marker", 
          col = c("blue4", "orange3"), suggestiveline = -log10(1e-04), logp = TRUE)
```

9. Zinterpretuj wyniki:
Jakie markery SNP są potencjalnie istotne?
Jakie wnioski można wyciągnąć na podstawie analizy PCA i wyników GWAS?
