# Wprowadzenie do Bioconductora

------------------------------------------------------------------------

## Spis treści

1.  [Czym jest Bioconductor?](#1-czym-jest-bioconductor)
2.  [Instalacja Bioconductora](#2-instalacja-bioconductora)
3.  [Najpopularniejsze pakiety Bioconductora](#3-podstawowe-pakiety-bioconductora)
4.  [Przykładowe zastosowania](#5-przykładowe-zastosowania)
    -   [4.1. Praca z sekwencjami](#41-praca-z-sekwencjami)
    -   [4.2. Analiza danych NGS](#42-analiza-danych-ngs)
    -   [4.3. Wizualizacja danych](#43-wizualizacja-danych)
5.  [Dodatkowe zasoby](#5-dodatkowe-zasoby)

------------------------------------------------------------------------

## 1. Czym jest Bioconductor?

[Bioconductor](https://bioconductor.org) to projekt open-source dostarczający zestaw narzędzi programistycznych w języku R, przeznaczonych do analizy i interpretacji danych genomicznych. Oferuje szeroką gamę pakietów dedykowanych różnym typom danych biologicznych, w tym danych z sekwencjonowania nowej generacji (NGS), mikromacierzy, proteomiki i wielu innych.

**Główne cechy Bioconductora:**

-   **Integracja z R:** Wykorzystuje moc języka R do analizy statystycznej i wizualizacji danych.
-   **Reprodukowalność:** Ułatwia tworzenie skryptów i raportów, co zwiększa reprodukowalność analiz.
-   **Społeczność:** Aktywna społeczność naukowa tworząca i utrzymująca pakiety.
-   **Aktualność:** Regularne aktualizacje pakietów i dokumentacji.

------------------------------------------------------------------------

## 2. Instalacja Bioconductora

Bioconductor jest zainstalowany jako zbiór pakietów R. Aby zainstalować Bioconductor, wykonaj następujące kroki w konsoli R lub RStudio.

### Krok 1: Instalacja pakietu BiocManager

Pakiet `BiocManager` jest zalecanym sposobem instalacji pakietów Bioconductor.

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) # ta linia kodu oszczędza czas i zasoby obliczeniowe
    install.packages("BiocManager")
```

### Krok 2: Instalacja podstawowych pakietów Bioconductora

Po zainstalowaniu `BiocManager`, możesz zainstalować podstawowe pakiety:

``` r
BiocManager::install()
```

### Krok 3: Instalacja konkretnych pakietów

Możesz zainstalować konkretne pakiety, np.:

``` r
BiocManager::install("GenomicFeatures")
BiocManager::install("AnnotationDbi")
BiocManager::install("ShortRead")
BiocManager::install("Biostrings")
```

Możesz też zainstalować kilka pakietów naraz wykorzystując wektor:

```r
BiocManager::install(c("GenomicFeatures", "AnnotationDbi", 
"ShortRead", "Biostrings"))
```

------------------------------------------------------------------------

## 3. Najpopularniejsze pakiety Bioconductora

Bioconductor oferuje wiele pakietów do różnych zastosowań. Oto niektóre z nich:

-   **Biostrings:** Operacje na sekwencjach DNA, RNA i białek.
-   **ShortRead:** Praca z danymi z sekwencjonowania nowej generacji (NGS).
-   **GenomicRanges:** Reprezentacja i manipulacja zakresami genomowymi.
-   **Rsamtools:** Interakcja z plikami BAM, SAM i indeksami.
-   **GenomicFeatures:** Tworzenie i używanie obiektów TxDb zawierających informacje o cechach genomowych.
-   **AnnotationDbi:** Interfejs do baz danych z annotacjami biologicznymi.
-   **DESeq2:** Analiza różnicowej ekspresji genów.
-   **edgeR:** Analiza danych z eksperymentów RNA-Seq.

------------------------------------------------------------------------

## 4. Przykładowe zastosowania

### 4.1. Praca z sekwencjami

#### Wczytywanie sekwencji DNA z pliku FASTA

``` r
library(Biostrings)

# Wczytanie sekwencji z pliku FASTA
dna_seqs <- readDNAStringSet("path/to/example.fasta")

# Wyświetlenie pierwszych kilku sekwencji
dna_seqs[1:5]
```

#### Podstawowe operacje na sekwencjach

``` r
# Odwrócenie i komplementarność sekwencji
rev_comp_seqs <- reverseComplement(dna_seqs)
rev_comp_seqs[1:5]

# Obliczanie zawartości GC
gc_content <- letterFrequency(dna_seqs, letters = c("G", "C"), as.prob = TRUE)
gc_content
```

### 5.2. Analiza danych NGS

#### Wczytywanie plików FASTQ i kontrola jakości

``` r
library(ShortRead)

# Ścieżka do pliku FASTQ
fastq_file <- "/path/tp/file/SRR31136237.fastq"

# Wczytanie odczytów
reads <- readFastq(fastq_file)

# Podstawowe informacje
length(reads)  # liczba odczytów
reads[1]       # pierwszy odczyt

# Kontrola jakości
qa_results <- qa(fastq_file, type = "fastq")
report(qa_results, dest = "/selected/path/QA_report.html")
```

#### Mapowanie odczytów do genomu referencyjnego (przykład z użyciem Rsubread)

``` r
# Instalacja pakietu Rsubread
BiocManager::install("Rsubread")
library(Rsubread)

# Indeksowanie genomu referencyjnego
buildindex(basename = "reference_index", reference = "path/to/reference/genome/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna")

# Mapowanie odczytów
align(index = "reference_index",
      readfile1 = "path/to/SRR31136237.fastq",
      output_file = "path/to/aligned_reads.BAM")
```

### 5.3. Wizualizacja danych

#### Wizualizacja sekwencji

``` r
library(Biostrings)
library(ggplot2)

# Obliczenie częstości nukleotydów
nucleotide_freq <- alphabetFrequency(dna_seqs, as.prob = TRUE)

# Konwersja do ramki danych
nuc_df <- as.data.frame(nucleotide_freq)

# Dodanie identyfikatora sekwencji
nuc_df$Sequence <- rownames(nuc_df)

# Wykres słupkowy częstości nukleotydów dla pierwszej sekwencji
nuc_long <- reshape2::melt(nuc_df[1, c("A", "C", "G", "T")])

ggplot(nuc_long, aes(x = variable, y = value)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  labs(title = "Częstość nukleotydów w pierwszej sekwencji",
       x = "Nukleotyd",
       y = "Częstość")
```

------------------------------------------------------------------------

## 6. Dodatkowe zasoby

-   **Strona główna Bioconductor:** <https://www.bioconductor.org/>
-   **Dokumentacja i vignettes:** Każdy pakiet Bioconductora zawiera szczegółową dokumentację i vignettes z przykładami użycia.
-   **Forum Bioconductor:** <https://support.bioconductor.org/>
-   **Kursy online i tutoriale:**
    -   **Bioconductor Course Material:** <https://bioconductor.org/help/course-materials/>
    -   **Kanał Bioconductora na YouTube:** <https://www.youtube.com/user/bioconductor>

------------------------------------------------------------------------
