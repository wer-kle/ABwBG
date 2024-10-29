# Kontrola jakości (QC) danych NGS z wykorzystaniem pakietów Bioconductor

## **Zadanie 1: Instalacja i konfiguracja środowiska**

**Cel:** Upewnienie się, że środowisko R i Bioconductor jest poprawnie zainstalowane i skonfigurowane do analizy danych NGS.

**Polecenia:**

1. Zainstaluj najnowszą wersję R i RStudio na swoim komputerze (pomiń jeżeli pracujesz na komputerze uczelnianym).
2. Zainstaluj pakiet `BiocManager`, jeśli jeszcze go nie masz:

   ```R
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   ```
   
3. Za pomocą `BiocManager` zainstaluj pakiety:

   - `ShortRead`
   - `Rqc`
   - `Biostrings`

   ```R
   BiocManager::install(c("ShortRead", "Rqc", "Biostrings"))
   ```
   
4. Załaduj zainstalowane pakiety i sprawdź ich wersje:

   ```R
   library(ShortRead)
   library(Rqc)
   library(Biostrings)
   
   packageVersion("ShortRead")
   packageVersion("Rqc")
   packageVersion("Biostrings")
   ```
   
**Zadanie do wykonania:**

- Upewnij się, że wszystkie pakiety zostały poprawnie zainstalowane i załadowane bez błędów.

---

### **Zadanie 2: Pobieranie danych FASTQ z bazy SRA**

**Cel:** Nauka pobierania danych sekwencyjnych w formacie FASTQ z bazy NCBI SRA.

**Polecenia:**

1. Przejdź na stronę NCBI SRA: [https://www.ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra)
2. Wyszukaj dane sekwencyjne dla **Escherichia coli** z sekwencjonowania całogenomowego (WGS) na platformie Illumina z dostępnym plikiem w formacie FASTQ.
3. Wybierz dowolny wyświetlony rekord.
4. Pobierz plik FASTQ (np. `SRX26518123.fastq.gz`) bezpośrednio ze strony. *Dla chętnych*: użyj [SRA Toolkit](https://github.com/ncbi/sra-tools).

**Zadanie do wykonania:**

- Pobierz plik FASTQ i umieść go w wybranym folderze na komputerze lokalnym.

---

### **Zadanie 3: Wczytywanie danych FASTQ do R**

**Cel:** Nauka wczytywania pliku FASTQ do R za pomocą pakietu `ShortRead`.

**Polecenia:**

1. Wczytaj plik FASTQ do R:

   ```R
   fq_file <- "data/SRX26518123.fastq.gz"
   fq_reads <- readFastq(fq_file)
   ```
   
2. Sprawdź podstawowe informacje o danych:

   - Liczba odczytów:

     ```R
     length(fq_reads)
     ```
     
   - Podgląd pierwszych kilku odczytów:

     ```R
     fq_reads[1:5]
     ```
     
**Zadanie do wykonania:**

- Upewnij się, że dane zostały poprawnie wczytane.
- Zanotuj liczbę odczytów w pliku.

---

### **Zadanie 4: Generowanie raportu QC za pomocą ShortRead**

**Cel:** Przeprowadzić kontrolę jakości danych za pomocą pakietu `ShortRead` i wygenerować raport QC.

**Polecenia:**

1. Wygeneruj obiekt z wynikami kontroli jakości:

   ```R
   qa_results <- qa(fq_file, type = "fastq")
   ```
   
2. Utwórz raport QC:

   ```R
   report(qa_results, dest = "path/to/your/QA_report.html")
   ```
   
3. Otwórz wygenerowany raport (`index.html`) w przeglądarce.

**Zadanie do wykonania:**

- Przejrzyj raport QC i zanotuj kluczowe obserwacje dotyczące jakości danych.
- Przejrzyj wykresy z folderu /images. Czy sekwencje są dobrej jakości?

---

### **Zadanie 5: Analiza wyników kontroli jakości**

**Cel:** Interpretacja raportu QC i identyfikacja potencjalnych problemów w danych.

**Polecenia:**

1. Przeanalizuj następujące sekcje raportu:

   - **Quality Score Distribution**
   - **Sequence Length Distribution**
   - **Per Base Sequence Content**
   - **Adapter Content**
   
2. Odpowiedz na pytania:

   - Czy jakość baz spada na końcach odczytów?
   - Czy odczyty mają jednolitą długość?
   - Czy występują biasy sekwencyjne?
   - Czy obecne są sekwencje adapterów?
   
**Zadanie do wykonania:**

- Sporządź krótkie podsumowanie jakości danych na podstawie raportu.

---

### **Zadanie 6: Generowanie raportu QC za pomocą Rqc**

**Cel:** Porównanie wyników kontroli jakości uzyskane za pomocą innego pakietu Bioconductora.

**Polecenia:**

1. Wygeneruj obiekt z wynikami QC za pomocą `Rqc`:

   ```R
   rqc_results <- rqc(path = "data", pattern = "SRR400294.fastq.gz", sample = TRUE)
   ```
   
2. Utwórz raport QC:

   ```R
   rqcReport(rqc_results, outdir = "Rqc_report")
   ```
   
3. Otwórz raport (`index.html`) w przeglądarce.

**Zadanie do wykonania:**

- Porównaj raport z `Rqc` z raportem z `ShortRead`.
- Zanotuj dodatkowe informacje lub wizualizacje dostępne w raporcie `Rqc`.

---

### **Zadanie 7: Przycinanie odczytów na podstawie jakości**

**Cel:** Nauka przycinania (trimmingu) odczytów na podstawie wartości jakości baz.

**Polecenia:**

1. Zainstaluj i załaduj pakiet `Biostrings` (jeśli jeszcze nie jest załadowany).

2. Przytnij odczyty o niskiej jakości z końców odczytów:

   ```R
   trimmed_reads <- trimTailw(fq_reads, k = 2, a = "B", successive = TRUE)
   ```
   
   - **k**: liczba kolejnych baz o jakości poniżej progu.
   - **a**: znak jakości odpowiadający progu (np. "B" dla Phred33 ~ Q10).
   
3. Sprawdź, ile odczytów zostało przyciętych:

   ```R
   sum(width(trimmed_reads) < width(fq_reads))
   ```
   
**Zadanie do wykonania:**

- Przytnij odczyty i zanotuj, jaki procent odczytów został zmodyfikowany.

---

### **Zadanie 8: Filtracja odczytów o niskiej jakości**

**Cel:** Usunąć odczyty o ogólnie niskiej jakości lub zbyt krótkie.

**Polecenia:**

1. Ustal minimalną akceptowalną długość odczytu, np. 50 bp.

2. Filtrowanie odczytów:

   ```R
   filtered_reads <- trimmed_reads[width(trimmed_reads) >= 50]
   ```
   
3. Sprawdź liczbę odczytów przed i po filtracji:

   ```R
   length(trimmed_reads)
   length(filtered_reads)
   ```
   
**Zadanie do wykonania:**

- Oblicz, jaki procent odczytów został odrzucony podczas filtracji.

---

### **Zadanie 9: Ponowna kontrola jakości po przycinaniu i filtracji**

**Cel:** Sprawdzić, jak przycinanie i filtracja wpłynęły na jakość danych.

**Polecenia:**

1. Zapisz przetworzone odczyty do nowego pliku FASTQ:

   ```R
   writeFastq(filtered_reads, "data/SRR400294_processed.fastq.gz")
   ```
   
2. Wygeneruj nowy raport QC dla przetworzonych danych:

   ```R
   qa_results_processed <- qa("data/SRR400294_processed.fastq.gz", type = "fastq")
   report(qa_results_processed, dest = "QA_report_processed")
   ```
   
3. Porównaj raporty QC przed i po przetwarzaniu.

**Zadanie do wykonania:**

- Opisz zmiany w jakości danych po przycinaniu i filtracji.

---

### **Zadanie 10: Analiza rozkładu długości odczytów**

**Cel:** Zbadać, jak przycinanie wpłynęło na długość odczytów.

**Polecenia:**

1. Porównaj rozkład długości odczytów przed i po przycinaniu:

   ```R
   # Przed przycinaniem
   hist(width(fq_reads), breaks = 50, main = "Długość odczytów przed przycinaniem", xlab = "Długość (bp)")
   
   # Po przycinaniu
   hist(width(filtered_reads), breaks = 50, main = "Długość odczytów po przycinaniu", xlab = "Długość (bp)")
   ```
   
**Zadanie do wykonania:**

- Porównaj histogramy i opisz, jak przycinanie wpłynęło na długość odczytów.

---

### **Zadanie 11: Analiza zawartości GC**

**Cel:** Sprawdzić zawartość GC w odczytach i zidentyfikować ewentualne biasy.

**Polecenia:**

1. Oblicz zawartość GC dla oryginalnych odczytów:

   ```R
   gc_content <- letterFrequency(sread(fq_reads), letters = "GC", as.prob = TRUE)
   ```
   
2. Wyświetl histogram zawartości GC:

   ```R
   hist(gc_content, breaks = 50, main = "Zawartość GC w oryginalnych odczytach", xlab = "Procent GC")
   ```
   
3. Powtórz analizę dla przetworzonych odczytów.

**Zadanie do wykonania:**

- Porównaj zawartość GC przed i po przetwarzaniu i zanotuj obserwacje.

---

### **Zadanie 12: Wykrywanie i usuwanie sekwencji adapterów**

**Cel:** Zidentyfikować obecność sekwencji adapterów i usunąć je z odczytów.

**Polecenia:**

1. Użyj pakietu `Biostrings` do wyszukania sekwencji adapterów w odczytach.

   - Zdefiniuj sekwencję adaptera (np. dla Illumina):

     ```R
     adapter_seq <- DNAString("AGATCGGAAGAGC")
     ```
     
   - Wyszukaj adaptery w odczytach:

     ```R
     match_positions <- vmatchPattern(adapter_seq, sread(fq_reads))
     ```
     
2. Przytnij odczyty zawierające adaptery:

   ```R
   # Funkcja do przycinania odczytów do pozycji początku adaptera
   trim_adapters <- function(reads, matches) {
     for (i in seq_along(reads)) {
       if (length(matches[[i]]) > 0) {
         end_pos <- start(matches[[i]][1]) - 1
         reads[[i]] <- subseq(reads[[i]], start = 1, end = end_pos)
       }
     }
     return(reads)
   }
   
   fq_reads_trimmed <- trim_adapters(fq_reads, match_positions)
   ```
   
**Zadanie do wykonania:**

- Przytnij odczyty zawierające adaptery i sprawdź, ile odczytów zostało zmodyfikowanych.

---

### **Zadanie 13: Automatyczne przycinanie adapterów za pomocą pakietu**

**Cel:** Wykorzystać dedykowany pakiet do automatycznego usuwania adapterów.

**Polecenia:**

1. Zainstaluj pakiet `Cutadapt` dla R (interfejs do programu Cutadapt) lub użyj funkcji `trimLRPatterns` z pakietu `Biostrings`.

2. Przytnij adaptery z odczytów:

   ```R
   # Przy użyciu funkcji trimLRPatterns
   fq_reads_trimmed <- trimLRPatterns(Lpattern = adapter_seq, subject = sread(fq_reads))
   ```
   
**Zadanie do wykonania:**

- Porównaj wyniki z poprzednim zadaniem i zanotuj różnice.

---

### **Zadanie 14: Generowanie zbiorczego raportu QC dla wielu próbek**

**Cel:** Nauczyć się generować raport QC dla wielu plików FASTQ jednocześnie.

**Polecenia:**

1. Pobierz dodatkowe pliki FASTQ (np. SRR400295, SRR400296) i umieść je w folderze `data`.

2. Wygeneruj raport QC dla wszystkich plików:

   ```R
   fq_files <- list.files(path = "data", pattern = "SRR40029[4-6].fastq.gz", full.names = TRUE)
   qa_results <- qa(fq_files, type = "fastq")
   report(qa_results, dest = "QA_report_multi")
   ```
   
**Zadanie do wykonania:**

- Przejrzyj zbiorczy raport i porównaj jakość danych między próbkami.

---

### **Zadanie 15: Dokumentacja i raportowanie wyników**

**Cel:** Przygotować kompletny raport z przeprowadzonych analiz QC.

**Polecenia:**

1. Utwórz dokument w formacie R Markdown (`QC_analysis.Rmd`), w którym umieścisz:

   - Opis celów i metodologii.
   - Kody R użyte w analizach.
   - Wyniki w postaci wykresów i tabel.
   - Interpretację wyników i wnioski.
   
2. Wygeneruj raport w formacie HTML lub PDF.

**Zadanie do wykonania:**

- Przygotuj przejrzysty i kompletny raport z analiz QC.

---

## **Dodatkowe zadania dla chętnych**

### **Zadanie 16: Wykorzystanie narzędzia FastQC**

**Cel:** Przeprowadzić kontrolę jakości danych za pomocą narzędzia FastQC i porównać wyniki z Bioconductorem.

**Polecenia:**

1. Pobierz i zainstaluj narzędzie FastQC: [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

2. Uruchom FastQC na pliku FASTQ:

   ```bash
   fastqc data/SRR400294.fastq.gz -o FastQC_report
   ```
   
3. Otwórz wygenerowany raport i przeanalizuj wyniki.

**Zadanie do wykonania:**

- Porównaj raport z FastQC z raportami wygenerowanymi za pomocą Bioconductora.
- Zanotuj różnice i podobieństwa w prezentowanych informacjach.

---

### **Zadanie 17: Zbiorcza analiza z użyciem MultiQC**

**Cel:** Wykorzystać narzędzie MultiQC do wygenerowania zbiorczego raportu QC dla wielu próbek.

**Polecenia:**

1. Zainstaluj MultiQC:

   ```bash
   pip install multiqc
   ```
   
2. Wykonaj analizę FastQC dla wszystkich plików FASTQ.

3. Uruchom MultiQC na folderze z wynikami FastQC:

   ```bash
   multiqc FastQC_report -o MultiQC_report
   ```
   
**Zadanie do wykonania:**

- Przejrzyj zbiorczy raport wygenerowany przez MultiQC.
- Porównaj sposób prezentacji danych z raportami z Bioconductora.

---

### **Zadanie 18: Wykorzystanie pakietu qrqc**

**Cel:** Przetestować inny pakiet R do kontroli jakości danych NGS.

**Polecenia:**

1. Zainstaluj pakiet `qrqc` (jeśli dostępny).

   ```R
   install.packages("qrqc")
   library(qrqc)
   ```
   
2. Wczytaj plik FASTQ i wygeneruj raport QC.

**Zadanie do wykonania:**

- Porównaj wyniki z innymi narzędziami i zanotuj swoje obserwacje.

---
