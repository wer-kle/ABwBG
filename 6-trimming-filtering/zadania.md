# Przycinanie (trimming) i filtrowanie danych NGS

## Zadanie 0: Wczytanie plików FASTQ

```{R}
library(ShortRead)
fq_reads1 <- readFastq("/Volumes/flashdrive/ecoli_simulated1_with_adapters.fq")
fq_reads2 <- readFastq("/Volumes/flashdrive/ecoli_simulated2_with_adapters.fq")
```

## **Zadanie 1: Przycinanie odczytów na podstawie jakości**

**Cel:** Nauka przycinania odczytów na podstawie wartości jakości baz.

**Polecenia:**

1. Przytnij bazy o niskiej jakości z końców odczytów:

```{r}
 # Przycinanie odczytów forward
   trimmed_reads1 <- trimTailw(fq_reads1, k = 2, a = "B", halfwidth = 1)
   
   # Przycinanie odczytów reverse
   trimmed_reads2 <- trimTailw(fq_reads2, k = 2, a = "B", halfwidth = 1)
```
   

2. Sprawdź, ile odczytów zostało przyciętych:

```{r}
length(fq_reads1)
length(trimmed_reads1)

length(fq_reads2)
length(trimmed_reads2)

sum(width(trimmed_reads1) < width(fq_reads1))
sum(width(trimmed_reads2) < width(fq_reads2))
```

**Uwaga**: komunikat ostrzegawczy 'Warning: longer object length is not a multiple of shorter object length' informuje o tym, że część odczytów została przycięta, a więc jest pożądany.

**Zadanie do wykonania:**

- Przytnij odczyty i zanotuj, jaki procent odczytów został zmodyfikowany.

---

## **Zadanie 2: Filtracja odczytów o niskiej jakości**

**Cel:** Usuwanie odczytów o ogólnie niskiej jakości lub zbyt krótkich.

**Polecenia:**

1. Ustal minimalną akceptowalną długość odczytu, np. 50 bp.

2. Filtrowanie odczytów:

```{R}
   # Filtrowanie odczytów forward
   filtered_reads1 <- trimmed_reads1[width(trimmed_reads1) >= 50]
   
   # Filtrowanie odczytów reverse
   filtered_reads2 <- trimmed_reads2[width(trimmed_reads2) >= 50]
```

3. Sprawdź liczbę odczytów przed i po filtracji:

```{R}
   # Odczyty forward
   length(trimmed_reads1)       # Po przycinaniu
   length(filtered_reads1)      # Po filtracji
   
   # Odczyty reverse
   length(trimmed_reads2)
   length(filtered_reads2)
```

**Zadanie do wykonania:**

- Oblicz, jaki procent odczytów został odrzucony podczas filtracji.

---

## **Zadanie 3: Ponowna kontrola jakości po przycinaniu i filtracji**

**Cel:** Sprawdzenie, jak przycinanie i filtracja wpłynęły na jakość danych.

**Polecenia:**

1. Zapisz przetworzone odczyty do nowych plików FASTQ:

```{R}
   writeFastq(filtered_reads1, "/Volumes/flashdrive/ecoli_simulated1_processed.fq")
   writeFastq(filtered_reads2, "/Volumes/flashdrive/ecoli_simulated2_processed.fq")
```

2. Wygeneruj nowe raporty QC dla przetworzonych danych:

```{R}
qa_results1 <- qa("/Volumes/flashdrive/ecoli_simulated1_with_adapters.fq", type = "fastq")   
qa_results1_processed <- qa( "/Volumes/flashdrive/ecoli_simulated1_processed.fq", type = "fastq")
report(qa_results1, dest = "/Volumes/flashdrive/QA_report_read1")
report(qa_results1_processed, dest = "/Volumes/flashdrive/QA_report_read1_processed")

qa_results2 <- qa("/Volumes/flashdrive/ecoli_simulated2_with_adapters.fq", type = "fastq")   
qa_results2_processed <- qa("/Volumes/flashdrive/ecoli_simulated2_processed.fq", type = "fastq")
report(qa_results2, dest = "/Volumes/flashdrive/QA_report_read2")
report(qa_results2_processed, dest = "/Volumes/flashdrive/QA_report_read2_processed")
```

**Zadanie do wykonania:**

- Opisz zmiany w jakości danych po przycinaniu i filtracji.

---

## **Zadanie 4: Analiza rozkładu długości odczytów**

**Cel:** Zbadanie, jak przycinanie wpłynęło na długość odczytów.

**Polecenia:**

1. Porównaj rozkład długości odczytów przed i po przycinaniu:

```{R}
   # Przed przycinaniem (odczyty forward)
   hist(width(fq_reads1), breaks = 50, main = "Długość odczytów forward przed przycinaniem", xlab = "Długość (bp)")
   
   # Po przycinaniu (odczyty forward)
   hist(width(filtered_reads1), breaks = 50, main = "Długość odczytów forward po przycinaniu", xlab = "Długość (bp)")
   
   # Przed przycinaniem (odczyty reverse)
   hist(width(fq_reads2), breaks = 50, main = "Długość odczytów reverse przed przycinaniem", xlab = "Długość (bp)")
   
   # Po przycinaniu (odczyty forward)
   hist(width(filtered_reads2), breaks = 50, main = "Długość odczytów reverse po przycinaniu", xlab = "Długość (bp)")
   
```

**Zadanie do wykonania:**

- Porównaj histogramy i opisz, jak przycinanie wpłynęło na długość odczytów.

---

## **Zadanie 5: Wykrywanie i usuwanie sekwencji adapterów**

**Cel:** Zidentyfikowanie obecności sekwencji adapterów i usunięcie ich z odczytów.

**Polecenia:**

1. Zdefiniuj sekwencję adaptera (np. dla Illumina):

```{R}
  library(Biostrings)
   adapter_seq <- DNAString("AGATCGGAAGAGC")
```

2. Przytnij adaptery z odczytów:

```{R}

# Przycinanie adapterów z odczytów forward:
trimmed_reads1_adapt <- trimLRPatterns(
  Lpattern = adapter_seq,
  subject = filtered_reads1
)

# Defuniujemy odczyty po przycięciu adapterów:
filtered_reads1 <- trimmed_reads1_adapt

# Przycinanie adapterów z odczytów reverse:
trimmed_reads2_adapt <- trimLRPatterns(
  Lpattern = adapter_seq,
  subject = filtered_reads2
)

# Defuniujemy odczyty po przycięciu adapterów:
filtered_reads2 <- trimmed_reads2_adapt

```

3. Sprawdź efekty przycinania:

```{R}
# Porównaj długości przed i po przycięciu adapterów
length(filtered_reads1)
length(trimmed_reads1)

length(filtered_reads2)
length(trimmed_reads2)

# Sprawdź ile odczytów zostało zmodyfikowanych
   sum(width(filtered_reads1) < width(trimmed_reads1))
   sum(width(filtered_reads2) < width(trimmed_reads2))
```

**Zadanie do wykonania:**

- Przytnij odczyty zawierające adaptery i sprawdź, ile odczytów zostało zmodyfikowanych.

---

## **Zadanie 6: Ponowna kontrola jakości po usunięciu adapterów**

**Cel:** Sprawdzenie, jak usunięcie adapterów wpłynęło na jakość danych.

**Polecenia:**

1. Zapisz odczyty po usunięciu adapterów:

```{R}
   writeFastq(filtered_reads1, "/Volumes/flashdrive/ecoli_simulated1_final.fq")
   writeFastq(filtered_reads2, "/Volumes/flashdrive/ecoli_simulated2_final.fq")
```

2. Wygeneruj ostateczne raporty QC:

```{R}
   qa_results1_final <- qa("/Volumes/flashdrive/ecoli_simulated1_final.fq", type = "fastq")
   report(qa_results1_final, dest = "/Volumes/flashdrive/QA_report_read1_final")
   
   qa_results2_final <- qa("/Volumes/flashdrive/ecoli_simulated2_final.fq", type = "fastq")
   report(qa_results2_final, dest = "/Volumes/flashdrive/QA_report_read2_final")
```

3. Porównaj raporty QC przed i po usunięciu adapterów.

**Zadanie do wykonania:**

- Opisz zmiany w jakości danych po usunięciu adapterów.

---

## **Zadanie 7: Dokumentacja i raportowanie wyników**

**Cel:** Przygotowanie kompletnego raportu z przeprowadzonych operacji przycinania i filtracji.

**Polecenia:**

1. Uzupełnij dokument R Markdown (`QC_and_Trimming_Report.Rmd`) o:

   - Kody R użyte w przycinaniu i filtracji.
   - Wyniki w postaci wykresów i tabel.
   - Interpretację wyników i wnioski.

2. Wygeneruj raport w formacie HTML lub PDF.

**Zadanie do wykonania:**

- Przygotuj przejrzysty i kompletny raport z procesu przycinania i filtracji.

---

## **Podsumowanie**

Przeprowadziliśmy pełny proces kontroli jakości oraz przycinania i filtracji danych NGS z wykorzystaniem symulowanych odczytów **Escherichia coli**. Dzięki temu przećwiczyliśmy w praktyce:

- Wczytywanie i podstawową analizę danych sekwencyjnych w R.
- Generowanie i interpretację raportów QC.
- Przycinanie odczytów na podstawie jakości i usuwanie sekwencji adapterów.
- Dokumentację i raportowanie wyników analiz bioinformatycznych.

---

**Dodatkowe zadania dla chętnych:**

- **Zadanie 8:** Wykorzystanie narzędzia FastQC do analizy jakości danych i porównanie wyników z raportami z Bioconductora.
- **Zadanie 9:** Przeprowadzenie zbiorczej analizy QC dla wielu próbek z użyciem MultiQC.

---
