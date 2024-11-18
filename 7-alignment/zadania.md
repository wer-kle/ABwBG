# Zadania praktyczne

## Zadanie 1: Przygotowanie danych sekwencyjnych

a) **Pobranie plików FASTQ**

Pobierz plik FASTQ z sekwencjami *E. coli* pochodzący z publicznej bazy danych NCBI SRA, uzyskanych z platformy Illumina.

b) **Analiza jakości odczytów**

Zaimportuj pliki FASTQ do R przy użyciu pakietu **ShortRead** i przeprowadź analizę jakości odczytów.

```R
library(ShortRead)

# Zaimportuj plik FASTQ
fq1 <- readFastq("ścieżka/do/pliku.fastq")

# Przeprowadź analizę jakości
qa_result <- qa(fq1)
report(qa_result, dest="ścieżka/do/QA_Report")
```

## Zadanie 2: Przygotowanie genomu referencyjnego

a) **Pobranie genomu referencyjnego**

Pobierz genom referencyjny *E. coli* z bazy danych NCBI lub Ensembl.

```R
# Przykładowe pobranie genomu z NCBI
download.file("ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz",
              destfile = "/ścieżka/do/ecoli_genome.fna.gz")
```

b) **Import genomu referencyjnego do R**

Zaimportuj sekwencję genomu do R przy użyciu pakietu **Biostrings**.

```R
library(Biostrings)

ref_genome <- readDNAStringSet("ścieżka/do/ecoli_genome.fna.gz")
```

## Zadanie 3: Indeksowanie genomu referencyjnego

a) **Budowanie indeksu genomu**

Użyj funkcji `buildindex()` z pakietu **Rsubread** do zindeksowania genomu referencyjnego.

```R
library(Rsubread)

buildindex(basename = "ecoli_index", reference = "ścieżka/do/ecoli_genome.fna.gz")
```

b) **Analiza wyników indeksowania**

Sprawdź w katalogu roboczym, jakie pliki zostały utworzone. *Dla chętnych:* Opisz, co zawiera każdy z nich.

## Zadanie 4: Mapowanie odczytów do genomu referencyjnego

a) **Wykonanie mapowania**

Wykorzystaj funkcję `align()` z pakietu **Rsubread** do mapowania odczytów do genomu referencyjnego.

```R
align(index = "ecoli_index",
      readfile1 = "ścieżka/do/pliku/fq1.fastq",
      input_format = "FASTQ",
      output_file = "ścieżka/do/aligned_sample.BAM")
```

b) **Wstępna analiza wyników mapowania**

- Oblicz procent poprawnie zmapowanych odczytów.
- Oblicz procent odczytów, które nie zostały zmapowane.
- Zastanów się nad możliwymi przyczynami niezmapowania odczytów.

## Zadanie 5: Analiza wyników mapowania

a) **Import zmapowanych odczytów**

Zaimportuj plik BAM do R przy użyciu pakietu **GenomicAlignments**.

```R
library(GenomicAlignments)

aln <- readGAlignments("ścieżka/do/aligned_sample.BAM")
```

b) **Obliczenie pokrycia genomu**

Oblicz pokrycie genomu i zidentyfikuj regiony o najwyższym i najniższym pokryciu.

```R
coverage_data <- coverage(aln)
```

c) **Wizualizacja pokrycia**

Zwizualizuj pokrycie początkowych pozycji genomu na wykresie przy użyciu pakietu **ggplot2**.

```R
library(ggplot2)

# Konwersja pokrycia do data frame
cov_df <- as.data.frame(coverage_data[[1]])
cov_df$position <- as.numeric(rownames(cov_df))

# Wykres pokrycia
pdf("ścieżka/do/pliku.pdf", width = 8, height = 6)

ggplot(cov_df[1:25000, ], aes(x = position, y = value)) +
  geom_line(color = "blue") +
  labs(title = "Pokrycie genomu E. coli",
       x = "Pozycja w genomie",
       y = "Liczba zmapowanych odczytów")
       
dev.off()
```

## Zadanie 6: Dokumentacja i raportowanie

a) **Sporządzenie raportu**

Przygotuj raport podsumowujący:

- Opis użytych danych i metod.
- Wyniki analizy jakości odczytów.
- Statystyki mapowania.
- Analizę pokrycia genomu.
- Wnioski i obserwacje.

## Zadanie 7 (dla chętnych): Porównanie narzędzi do mapowania

a) **Mapowanie przy użyciu alternatywnego narzędzia**

Wykonaj mapowanie odczytów za pomocą narzędzia **Bowtie2** lub **BWA**, korzystając z interfejsu R lub zewnętrznie.

b) **Porównanie wyników**

Porównaj wyniki z uzyskanymi wcześniej:

- Czy liczba zmapowanych odczytów różni się między narzędziami?
- Jakie są różnice w pokryciu genomu?
- Jakie mogą być przyczyny tych różnic?

*Wskazówka:* Możesz użyć pakietu **Rbowtie2** lub uruchomić narzędzia z linii poleceń i zaimportować wyniki do R.
