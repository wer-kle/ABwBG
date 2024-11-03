# Kontrola jakości (QC) danych NGS z wykorzystaniem pakietów Bioconductor

## **Zadanie 1: Instalacja i konfiguracja środowiska**

**Cel:** Upewnienie się, że środowisko R i Bioconductor są poprawnie zainstalowane i skonfigurowane do analizy danych NGS.

**Polecenia:**

1. Upewnij się, że korzystasz z najnowszej wersji R i RStudio.
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
   BiocManager::install(c("ShortRead", "Rqc", "Biostrings"), force = TRUE)
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

## **Zadanie 2: Pobieranie danych FASTQ z bazy SRA**

**Cel:** Nauka pobierania danych sekwencyjnych w formacie FASTQ z bazy NCBI SRA.

**Polecenia:**

1. Przejdź na stronę NCBI SRA: [https://www.ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra)
2. Wyszukaj dane sekwencyjne dla **Escherichia coli** z sekwencjonowania całogenomowego (WGS) na platformie Illumina z dostępem do surowych plików FASTQ.
3. Wybierz dowolny wyświetlony rekord, który zawiera surowe dane sekwencyjne.
4. Pobierz plik FASTQ (np. `SRRXXXXXXX.fastq.gz`) bezpośrednio ze strony. *Dla chętnych*: użyj [SRA Toolkit](https://github.com/ncbi/sra-tools).

**Zadanie do wykonania:**

- Pobierz plik FASTQ i umieść go w wybranym folderze na komputerze lokalnym.

---

## **Zadanie 3: Wczytywanie danych FASTQ do R**

**Cel:** Nauka wczytywania pliku FASTQ do R za pomocą pakietu `ShortRead`.

**Polecenia:**

1. Wczytaj plik FASTQ do R:

   ```R
   fq_file <- "ścieżka/do/twojego/pliku.fastq.gz" 
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

## **Zadanie 4: Generowanie raportu QC za pomocą ShortRead**

**Cel:** Przeprowadzenie kontroli jakości danych za pomocą pakietu `ShortRead` i wygenerowanie raportu QC.

**Polecenia:**

1. Wygeneruj obiekt z wynikami kontroli jakości:

   ```R
   qa_results <- qa(fq_file, type = "fastq")
   ```
   
2. Utwórz raport QC:

   ```R
   report(qa_results, dest = "ścieżka/do/folderu/z/raportem")
   ```
   
3. Otwórz wygenerowany raport (`index.html`) w przeglądarce.

**Zadanie do wykonania:**

- Przejrzyj raport QC i zanotuj kluczowe obserwacje dotyczące jakości danych.
- Przejrzyj wykresy z folderu `/images`. Czy sekwencje są dobrej jakości?

---

## **Zadanie 5: Analiza wyników kontroli jakości**

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

## **Zadanie 6: Generowanie raportu QC za pomocą Rqc**

**Cel:** Porównanie wyników kontroli jakości uzyskanych za pomocą innego pakietu Bioconductora.

**Polecenia:**

1. Wygeneruj obiekt z wynikami QC za pomocą `Rqc`:

   ```R
   rqc_results <- rqc(path = "ścieżka/do/folderu/z/danymi", pattern = "nazwa_pliku.fastq.gz", sample = TRUE)
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

## **Zadanie 7: Analiza zawartości GC**

**Cel:** Sprawdzenie zawartości GC w odczytach i zidentyfikowanie ewentualnych biasów.

**Polecenia:**

1. Oblicz zawartość GC dla oryginalnych odczytów:

   ```R
   gc_content <- letterFrequency(sread(fq_reads), letters = "GC", as.prob = TRUE)
   ```
   
2. Wyświetl histogram zawartości GC:

   ```R
   hist(gc_content, breaks = 50, main = "Zawartość GC w oryginalnych odczytach", xlab = "Procent GC")
   ```
   
**Zadanie do wykonania:**

- Przeanalizuj histogram i zanotuj, czy zawartość GC jest zgodna z oczekiwaniami dla **Escherichia coli**.

---

## **Zadanie 8: Generowanie zbiorczego raportu QC dla wielu próbek**

**Cel:** Nauka generowania raportu QC dla wielu plików FASTQ jednocześnie.

**Polecenia:**

1. Pobierz dwa dodatkowe pliki FASTQ (np. `SRRXXXXXXX.fastq.gz`) i umieść je w folderze z danymi.

2. Wygeneruj raport QC dla wszystkich plików:

   ```R
   fq_files <- list.files(path = "ścieżka/do/folderu/z/danymi", pattern = "SRRXXXXXXX.fastq.gz", full.names = TRUE)
   qa_results <- qa(fq_files, type = "fastq")
   report(qa_results, dest = "QA_report_multi")
   ```
   
**Zadanie do wykonania:**

- Przejrzyj zbiorczy raport i porównaj jakość danych między próbkami.

---

## **Zadanie 9: Dokumentacja i raportowanie wyników**

**Cel:** Przygotowanie kompletnego raportu z przeprowadzonych analiz QC.

**Polecenia:**

1. Utwórz dokument w formacie R Markdown (`QC_analysis.Rmd`), w którym umieścisz:

   - Opis celów i metodologii.
   - Kody R użyte w analizach.
   - Wyniki w postaci wykresów i tabel.
   - Interpretację wyników i wnioski.
   
2. Wygeneruj raport w formacie HTML lub PDF.

**Zadanie do wykonania:**

- Przygotuj przejrzysty i kompletny raport z analiz QC.
