# Tutorial korzystania z bazy danych SRA (Sequence Read Archive)

**Spis treści**
1. **Wprowadzenie do SRA**
2. **Rejestracja i logowanie**
3. **Nawigacja po stronie SRA**
4. **Wyszukiwanie danych**
5. **Korzystanie z SRA Toolkit**
   - 5.1 **Instalacja SRA Toolkit**
   - 5.2 **Konfiguracja SRA Toolkit**
   - 5.3 **Podstawowe narzędzia i polecenia**
   - 5.4 **Przykłady użycia**
   - 5.5 **Rozwiązywanie problemów**
6. **Analiza danych**
7. **Porady i najlepsze praktyki**

---

## 1. Wprowadzenie do SRA

**Sequence Read Archive (SRA)** to największa publiczna baza danych surowych odczytów sekwencji DNA i RNA. Zarządzana przez **National Center for Biotechnology Information (NCBI)**, służy jako centralne repozytorium danych z eksperymentów sekwencjonowania wysokoprzepustowego.

## 2. Rejestracja i logowanie

Choć wiele funkcji SRA jest dostępnych bez rejestracji, założenie konta NCBI umożliwia korzystanie z dodatkowych opcji, takich jak zapisywanie wyszukiwań czy tworzenie list ulubionych.

**Kroki:**
1. Przejdź na stronę [NCBI](https://www.ncbi.nlm.nih.gov/).
2. Kliknij **"Sign in to NCBI"** w prawym górnym rogu.
3. Wybierz **"Register for an NCBI account"**.
4. Wypełnij formularz rejestracyjny i potwierdź rejestrację poprzez e-mail.

## 3. Nawigacja po stronie SRA

Strona SRA jest dostępna pod adresem: [https://www.ncbi.nlm.nih.gov/sra](https://www.ncbi.nlm.nih.gov/sra).

**Główne sekcje:**
- **Pasek wyszukiwania**: umożliwia wyszukiwanie danych wg różnych kryteriów.
- **Menu nawigacyjne**: dostęp do różnych zasobów i narzędzi NCBI.
- **Linki pomocnicze**: dokumentacja, pomoc techniczna, FAQ.

## 4. Wyszukiwanie danych

**Proste wyszukiwanie:**
- Wpisz słowa kluczowe w pasku wyszukiwania (np. nazwa organizmu, gen, numer dostępu).
- Użyj filtrów po lewej stronie, aby zawęzić wyniki (np. wg typu eksperymentu, platformy sekwencjonowania).

**Zaawansowane wyszukiwanie:**
- Skorzystaj z **"Advanced Search"** pod paskiem wyszukiwania.
- Używaj operatorów logicznych (AND, OR, NOT) i określ pola wyszukiwania.

**Przykład:**
Aby znaleźć odczyty RNA-seq dla *Homo sapiens*:
```
"Homo sapiens"[Organism] AND "RNA-Seq"[Strategy]
```

## 5. Korzystanie z SRA Toolkit

**SRA Toolkit** to zestaw narzędzi umożliwiających pobieranie i przetwarzanie danych z bazy SRA bezpośrednio z linii komend. Jest niezbędny dla efektywnego pobierania dużych zestawów danych i integracji z pipeline'ami bioinformatycznymi.

### 5.1 Instalacja SRA Toolkit

**Dostępne systemy operacyjne:**
- **Linux**
- **macOS**
- **Windows**

**Instalacja na Linux i macOS:**

1. **Pobierz najnowszą wersję SRA Toolkit:**
   - Przejdź na stronę [SRA Toolkit Release Page](https://github.com/ncbi/sra-tools/wiki/Downloads).
   - Wybierz odpowiednią wersję dla Twojego systemu.

2. **Rozpakuj pobrany plik:**
   ```bash
   tar -xvzf sratoolkit.<wersja>-<platforma>.tar.gz
   ```

3. **Dodaj SRA Toolkit do zmiennej środowiskowej PATH:**
   - Edytuj plik `~/.bashrc` lub `~/.bash_profile` i dodaj:
     ```bash
     export PATH=$PATH:/ścieżka/do/sratoolkit.<wersja>-<platforma>/bin
     ```
   - Zastosuj zmiany:
     ```bash
     source ~/.bashrc
     ```

**Instalacja na Windows:**

1. **Pobierz instalator MSI** z [oficjalnej strony](https://github.com/ncbi/sra-tools/wiki/Downloads).
2. **Uruchom instalator** i postępuj zgodnie z instrukcjami.
3. **Dodaj ścieżkę do SRA Toolkit** w zmiennych środowiskowych systemu.

### 5.2 Konfiguracja SRA Toolkit

Po instalacji należy skonfigurować SRA Toolkit, aby poprawnie działał.

**Kroki konfiguracji:**

1. **Uruchom polecenie konfiguracji:**
   ```bash
   vdb-config --interactive
   ```
2. **Ustaw katalog cache:**
   - W menu konfiguracji wybierz opcję **"Preferences"**.
   - Ustaw lokalizację katalogu cache, jeśli chcesz zmienić domyślne ustawienia.

3. **Akceptuj warunki licencyjne:**
   - Upewnij się, że akceptujesz warunki użytkowania danych SRA.

**Alternatywna konfiguracja bez interakcji:**

Jeśli chcesz skonfigurować SRA Toolkit bez interfejsu interaktywnego:

```bash
vdb-config --accept-gcpw-license
```

### 5.3 Podstawowe narzędzia i polecenia

**Główne narzędzia SRA Toolkit:**

- **prefetch**: Pobiera dane z SRA i zapisuje je w formacie SRA.
- **fastq-dump**: Konwertuje pliki SRA na format FASTQ.
- **fasterq-dump**: Szybsza alternatywa dla fastq-dump, zalecana dla dużych danych.
- **sam-dump**: Konwertuje pliki SRA na format SAM.

### 5.4 Przykłady użycia

**1. Pobieranie danych przy użyciu prefetch:**

Pobierz dane dla konkretnego identyfikatora, np. SRR000000:

```bash
prefetch SRR000000
```

Pliki zostaną zapisane w domyślnym katalogu cache SRA (zazwyczaj w `$HOME/ncbi/public/sra`).

**2. Konwersja pliku SRA na FASTQ przy użyciu fasterq-dump:**

```bash
fasterq-dump SRR000000
```

Plik FASTQ zostanie zapisany w bieżącym katalogu.

**Opcje dodatkowe:**

- **Określenie katalogu wyjściowego:**

  ```bash
  fasterq-dump SRR000000 -O /ścieżka/do/katalogu
  ```

- **Włączenie kompresji gzip:**

  Ponieważ fasterq-dump nie obsługuje bezpośrednio kompresji, można użyć narzędzia **pigz** do równoległej kompresji:

  ```bash
  fasterq-dump SRR000000 -O /ścieżka/do/katalogu
  pigz /ścieżka/do/katalogu/SRR000000.fastq
  ```

**3. Pobieranie i konwersja w jednym kroku:**

Jeśli chcesz pobrać dane i od razu je skonwertować, możesz użyć fastq-dump z opcją pobierania:

```bash
fastq-dump --gzip --skip-technical --readids --read-filter pass --dumpbase --split-files --clip SRR000000
```

**Wyjaśnienie opcji:**

- `--gzip`: kompresuje wyjściowy plik FASTQ.
- `--skip-technical`: pomija odczyty techniczne.
- `--readids`: zachowuje oryginalne identyfikatory odczytów.
- `--read-filter pass`: filtruje odczyty o wysokiej jakości.
- `--dumpbase`: używa oryginalnych baz zamiast jakości w formacie kolorów.
- `--split-files`: dzieli odczyty parowane na oddzielne pliki.
- `--clip`: przycina adaptory i słabe jakościowo bazy.

### 5.5 Rozwiązywanie problemów

**1. Problem z dostępem do danych:**

Jeśli otrzymujesz komunikat o błędzie dotyczący braku dostępu do danych:

- Upewnij się, że zaakceptowałeś warunki licencyjne:

  ```bash
  vdb-config --accept-gcpw-license
  ```

- Sprawdź, czy masz aktualną wersję SRA Toolkit.

**2. Błędy podczas pobierania:**

- Sprawdź połączenie internetowe.
- Upewnij się, że identyfikator SRA jest poprawny.
- Spróbuj użyć opcji zwiększających szczegółowość komunikatów:

  ```bash
  prefetch SRR000000 --verbose
  ```

**3. Niska wydajność konwersji:**

- Użyj **fasterq-dump** zamiast **fastq-dump** dla lepszej wydajności.
- Zwiększ liczbę wątków:

  ```bash
  fasterq-dump SRR000000 -e 8
  ```

**4. Brak miejsca na dysku:**

- Sprawdź dostępne miejsce przed pobieraniem.
- Użyj opcji kompresji i usuwaj pliki tymczasowe.

## 6. Analiza danych

Po pobraniu danych możesz przystąpić do ich analizy za pomocą narzędzi bioinformatycznych:

- **Kontrola jakości**: np. **FastQC**.
- **Przycinanie adapterów**: np. **Trimmomatic**, **Cutadapt**.
- **Mapowanie odczytów**: np. **Bowtie2**, **HISAT2**, **BWA**.
- **Analiza ekspresji genów**: np. **Cufflinks**, **DESeq2**, **EdgeR**.

**Przykład pipeline'u dla danych RNA-Seq:**

1. **Kontrola jakości surowych odczytów:**

   ```bash
   fastqc SRR000000_1.fastq.gz SRR000000_2.fastq.gz
   ```

2. **Przycinanie niskiej jakości baz i adapterów:**

   ```bash
   trimmomatic PE SRR000000_1.fastq.gz SRR000000_2.fastq.gz \
   SRR000000_1_paired.fastq.gz SRR000000_1_unpaired.fastq.gz \
   SRR000000_2_paired.fastq.gz SRR000000_2_unpaired.fastq.gz \
   ILLUMINACLIP:adapters.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
   ```

3. **Mapowanie do genomu referencyjnego:**

   ```bash
   hisat2 -p 8 -x indeks_genomu -1 SRR000000_1_paired.fastq.gz -2 SRR000000_2_paired.fastq.gz -S SRR000000.sam
   ```

4. **Konwersja SAM do BAM i sortowanie:**

   ```bash
   samtools view -bS SRR000000.sam | samtools sort -o SRR000000_sorted.bam
   ```

5. **Analiza ekspresji genów:**

   Użyj narzędzi takich jak **featureCounts** do zliczania odczytów na gen.

## 7. Porady i najlepsze praktyki

- **Sprawdzaj metadane**: zawsze przeglądaj informacje towarzyszące zestawom danych, aby zrozumieć kontekst eksperymentu.
- **Uważaj na wersje referencyjne**: upewnij się, że używasz odpowiednich wersji genomów referencyjnych podczas mapowania.
- **Zarządzaj przestrzenią dyskową**: dane sekwencyjne mogą zajmować dużo miejsca; regularnie usuwaj niepotrzebne pliki.
- **Korzystaj z najnowszych wersji narzędzi**: regularnie aktualizuj SRA Toolkit i inne oprogramowanie.
- **Automatyzuj procesy**: rozważ pisanie skryptów lub używanie managerów workflow, takich jak **Snakemake** czy **Nextflow**.

---

**Przydatne linki:**

- [NCBI SRA Help](https://www.ncbi.nlm.nih.gov/sra/docs/)
- [SRA Toolkit Documentation](https://github.com/ncbi/sra-tools/wiki)
- [Forum NCBI](https://www.ncbi.nlm.nih.gov/support/discussion/)
- [Instrukcje SRA Toolkit](https://www.ncbi.nlm.nih.gov/sra/docs/sra-cloud-tools/)
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)
- [HISAT2](https://daehwankimlab.github.io/hisat2/)

Jeśli masz dodatkowe pytania lub napotkasz problemy, warto skonsultować się z społecznością bioinformatyczną lub skorzystać z pomocy technicznej NCBI.
