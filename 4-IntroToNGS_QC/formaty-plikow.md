# Wprowadzenie do BioConductor i formatów danych WGS

Analiza danych WGS (Whole Genome Sequencing) wymaga zrozumienia różnych formatów plików używanych do przechowywania i przekazywania informacji. Każdy format ma swoją specyficzną strukturę i przechowuje inne typy danych. Poniżej znajduje się szczegółowe omówienie kolumn w najważniejszych formatach danych WGS.

------------------------------------------------------------------------

## FASTQ - surowe odczyty sekwencji z informacją o jakości

### Struktura pliku FASTQ:

-   Każdy odczyt w pliku FASTQ jest reprezentowany przez **cztery linie**:
    1.  **Nagłówek odczytu** (zaczyna się od '\@')
    2.  **Sekwencja nukleotydowa**
    3.  **Separator** (zaczyna się od '+')
    4.  **Jakości sekwencji** (kodowane znaki ASCII)

**Opis poszczególnych linii:**

1.  **Linia 1 - Nagłówek odczytu (@SEQ_ID):**
    -   Zawiera unikalny identyfikator odczytu.
    -   Może zawierać dodatkowe informacje, takie jak pozycja w komorze przepływowej sekwenatora.
2.  **Linia 2 - Sekwencja nukleotydowa:**
    -   Ciąg liter reprezentujących nukleotydy: A, T, C, G, N.
    -   N reprezentuje nieokreślony nukleotyd (im więcej 'N', tym niższa jakość sekwencji).
3.  **Linia 3 - Separator (+):**
    -   Czasami zawiera komentarze lub powtarza nagłówek odczytu.
    -   Oddziela sekwencję od jakości.
4.  **Linia 4 - Jakości sekwencji:**
    -   Ciąg znaków ASCII reprezentujących jakość każdego nukleotydu.
    -   Jakość jest kodowana w skali Phred.

**Przykład odczytu w formacie FASTQ:**

```         
@D00360:58:C6JL9ANXX:2:1101:1515:2127 1:N:0:CGATGT 

GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAAT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**
```

**Wyjaśnienie:**

-   Nagłówek (```@D00360:58:C6JL9ANXX:2:1101:1515:2127 1:N:0:CGATGT```): Unikalny identyfikator odczytu. -- D00360: 'Instrument': nazwa sekwenatora -- 58: 'RunID': identyfikator przebiegu sekwencjonowania -- C6JL9ANXX: 'FlowCellID': idenryfikator komory przepływowej -- 2: 'Lane': numer rzędu komory przepływowej -- 1101: 'Tile': numer płytki -- 1515: 'X': współrzędna X klastra -- 2127: 'Y': współrzędna Y klastra -- 1:N:0:CGATGT: dodatkowe informacje = numer odczytu:czy odczyt przeszedł kontrolę jakości (Y - tak, N - nie):0 - numer indeksu:CGATGT - sekwencja indeksu
-   Sekwencja: Sekwencja DNA odczytu.
-   Separator (+): Oddziela sekwencję od jakości.
-   Jakość: Znaki ASCII reprezentujące jakość każdego nukleotydu.

------------------------------------------------------------------------

## SAM/BAM - zmapowane odczyty**

### SAM (Sequence Alignment/Map):

-   Tekstowy format z danymi o zmapowanych odczytach.
-   Każda linia (poza nagłówkiem) reprezentuje jeden zmapowany odczyt.
-   Kolumny są oddzielone tabulatorami.

**Kolumny w formacie SAM:**

1.  **QNAME (Query Name):**
    -   Nazwa odczytu (identyfikator).
2.  **FLAG:**
    -   Liczba całkowita reprezentująca różne właściwości odczytu (np. czy odczyt jest zmapowany, czy jest pierwszym z pary).
3.  **RNAME (Reference Name):**
    -   Nazwa referencyjnej sekwencji (np. chromosomu).
4.  **POS (Position):**
    -   Pozycja początkowa odczytu na referencji (1-bazowa).
5.  **MAPQ (Mapping Quality):**
    -   Jakość mapowania odczytu.
6.  **CIGAR:**
    -   Ciąg opisujący dopasowanie odczytu do referencji (np. liczbę dopasowań, wstawień).
7.  **RNEXT (Reference Name of the Mate/Next Read):**
    -   Nazwa referencji dla odczytu z pary.
8.  **PNEXT (Position of the Mate/Next Read):**
    -   Pozycja początkowa odczytu z pary.
9.  **TLEN (Template Length):**
    -   Długość fragmentu (od pierwszego odczytu do drugiego).
10. **SEQ (Sequence):**
    -   Sekwencja nukleotydowa odczytu.
11. **QUAL (Quality):**
    -   Jakości nukleotydów w odczycie.

**Opcjonalne pola:**

-   Dodatkowe informacje w formacie <TAG:TYPE:VALUE> (np. NM:i:1 - liczba niedopasowań).

**Przykład linii w pliku SAM:**

```         
r001    99  chr1    7   30  8M2I4M1D3M    =   37  39  TTAGATAAAGAGGATACTG  *   NM:i:1  AS:i:23 XS:i:17
```

**Wyjaśnienie kolumn:**

-   **QNAME (r001):** Identyfikator odczytu.
-   **FLAG (99):** Wskazuje, że odczyt jest zmapowany i jest pierwszym z pary.
-   **RNAME (chr1):** Odczyt zmapowany do chromosomu 1.
-   **POS (7):** Pozycja początkowa na chromosomie.
-   **MAPQ (30):** Jakość mapowania.
-   **CIGAR (8M2I4M1D3M):** Wzorzec dopasowania (8 dopasowań, 2 wstawienia itd.).
-   **RNEXT (=):** Drugi odczyt z pary zmapowany do tej samej referencji.
-   **PNEXT (37):** Pozycja drugiego odczytu z pary.
-   **TLEN (39):** Długość fragmentu.
-   **SEQ:** Sekwencja odczytu.
-   \*\*QUAL (\*):\*\* Brak informacji o jakości (oznaczone jako '\*').

**Format BAM:**

-   Binarny odpowiednik formatu SAM.
-   Zawiera te same informacje, ale jest skompresowany dla efektywnego przechowywania i przetwarzania.

------------------------------------------------------------------------

## **3. Format VCF - warianty genetyczne**

### Struktura pliku VCF:

-   Plik tekstowy z liniami nagłówkowymi i danymi wariantów.
-   Nagłówek zaczyna się od '\##' i definiuje metadane.
-   Linia z nazwami kolumn zaczyna się od '\#'.

### Kolumny w formacie VCF:

1.  **CHROM:**
    -   Chromosom, na którym znajduje się wariant.
2.  **POS (Position):**
    -   Pozycja wariantu na chromosomie (1-bazowa).
3.  **ID:**
    -   Identyfikator wariantu (np. z bazy dbSNP) lub '.' jeśli brak.
4.  **REF (Reference):**
    -   Nukleotyd(y) referencyjne.
5.  **ALT (Alternate):**
    -   Nukleotyd(y) alternatywne (wariant).
6.  **QUAL (Quality):**
    -   Jakość wariantu.
7.  **FILTER:**
    -   Informacje o filtrach (np. 'PASS' jeśli wariant przeszedł wszystkie filtry).
8.  **INFO:**
    -   Dodatkowe informacje o wariancie w formacie klucz=wartość.
9.  **FORMAT (opcjonalne):**
    -   Definicja formatu danych genotypowych.
10. **Dane próbek (opcjonalne):**
    -   Dane genotypowe dla poszczególnych próbek, zgodne z FORMAT.

**Przykład linii w pliku VCF:**

```         
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO                        FORMAT  Sample1
chr1    10177   rs367896724     A       AC      100     PASS    AC=1;AF=0.5;AN=2     GT:DP   0/1:14
```

**Wyjaśnienie kolumn:**

-   **CHROM (chr1):** Chromosom 1.
-   **POS (10177):** Pozycja wariantu.
-   **ID (rs367896724):** Identyfikator wariantu.
-   **REF (A):** Nukleotyd referencyjny.
-   **ALT (AC):** Wariant alternatywny (wstawienie C).
-   **QUAL (100):** Wysoka jakość wariantu.
-   **FILTER (PASS):** Wariant przeszedł wszystkie filtry.
-   **INFO (AC=1;AF=0.5;AN=2):** Dodatkowe informacje:
    -   **AC=1:** Liczba alleli alternatywnych w populacji.
    -   **AF=0.5:** Częstość allelu alternatywnego.
    -   **AN=2:** Liczba wszystkich alleli w populacji.
-   **FORMAT (GT:DP):** Definicja danych genotypowych:
    -   **GT:** Genotyp.
    -   **DP:** Głębokość pokrycia.
-   **Sample1 (0/1:14):** Dane dla próbki:
    -   **GT (0/1):** Heterozygota (jeden allel referencyjny, jeden alternatywny).
    -   **DP (14):** Głębokość pokrycia wynosi 14.

------------------------------------------------------------------------

## **4. Format BED - anotacje genomowe**

**Struktura pliku BED:**

-   Tekstowy format z minimalnie trzema kolumnami.
-   Kolumny są oddzielone tabulatorami.
-   Pozwala na przechowywanie informacji o regionach genomowych.

**Kolumny w formacie BED:**

1.  **chrom:**
    -   Nazwa chromosomu (np. 'chr1').
2.  **chromStart:**
    -   Pozycja początkowa regionu (0-bazowa).
3.  **chromEnd:**
    -   Pozycja końcowa regionu (nie włączając tej pozycji, 0-bazowa).

**Dodatkowe kolumny (opcjonalne):**

4.  **name:**
    -   Nazwa regionu (np. nazwa genu).
5.  **score:**
    -   Wartość numeryczna (0-1000), często używana do reprezentacji poziomu istotności.
6.  **strand:**
    -   Nić ('+' lub '-').
7.  **thickStart, thickEnd:**
    -   Pozycje używane w wizualizacji (np. do zaznaczenia regionu kodującego).
8.  **itemRgb:**
    -   Kolor w formacie RGB (np. '255,0,0' dla czerwonego).
9.  **blockCount, blockSizes, blockStarts:**
    -   Informacje o podregionach (np. eksonach w genie).

**Przykład linii w pliku BED:**

```         
chr7    127471196   127495720   uc010nxr.1  0   +   127471196   127495720   255,0,0
```

**Wyjaśnienie kolumn:**

-   **chrom (chr7):** Chromosom 7.
-   **chromStart (127471196):** Początek regionu.
-   **chromEnd (127495720):** Koniec regionu.
-   **name (uc010nxr.1):** Nazwa regionu.
-   **score (0):** Brak wartości (często '0' gdy niewykorzystane).
-   **strand (+):** Nić dodatnia.
-   **thickStart, thickEnd:** Używane w wizualizacji (tutaj takie same jak chromStart i chromEnd).
-   **itemRgb (255,0,0):** Kolor czerwony.

------------------------------------------------------------------------

##### **5. Format GFF/GTF - szczegółowe anotacje genomowe**

**Struktura pliku GFF/GTF:**

-   Tekstowy format z dziewięcioma kolumnami.
-   Używany do przechowywania szczegółowych informacji o genach i innych elementach genomu.

**Kolumny w formacie GFF/GTF:**

1.  **seqname:**
    -   Nazwa sekwencji (np. chromosomu).
2.  **source:**
    -   Źródło anotacji (np. program lub baza danych).
3.  **feature:**
    -   Typ elementu (np. 'gene', 'exon', 'CDS').
4.  **start:**
    -   Pozycja początkowa elementu (1-bazowa).
5.  **end:**
    -   Pozycja końcowa elementu.
6.  **score:**
    -   Wartość numeryczna (np. poziom pewności predykcji), '.' jeśli brak.
7.  **strand:**
    -   Nić ('+' lub '-').
8.  **frame:**
    -   Ramka odczytu (0, 1, 2), '.' jeśli nie dotyczy.
9.  **attribute:**
    -   Dodatkowe informacje w formacie klucz-wartość (np. ID, Name).

**Przykład linii w pliku GTF:**

```         
chr1    HAVANA  gene    11869   14412   .   +   .   gene_id "ENSG00000223972"; gene_name "DDX11L1";
```

**Wyjaśnienie kolumn:**

-   **seqname (chr1):** Chromosom 1.
-   **source (HAVANA):** Źródło anotacji.
-   **feature (gene):** Typ elementu - gen.
-   **start (11869):** Pozycja początkowa genu.
-   **end (14412):** Pozycja końcowa genu.
-   **score (.):** Brak wartości.
-   **strand (+):** Nić dodatnia.
-   **frame (.):** Nie dotyczy (dla genów).
-   **attribute:** Dodatkowe informacje:
    -   **gene_id "ENSG00000223972";** Identyfikator genu.
    -   **gene_name "DDX11L1";** Nazwa genu.

------------------------------------------------------------------------

#### **Podsumowanie znaczenia kolumn w formatach WGS**

-   **Dokładne zrozumienie zawartości kolumn jest kluczowe** dla poprawnej analizy danych genomicznych.
-   **Formaty różnią się strukturą i przeznaczeniem**:
    -   **FASTQ**: Surowe odczyty z informacją o jakości.
    -   **SAM/BAM**: Zmapowane odczyty z informacjami o dopasowaniu.
    -   **VCF**: Informacje o wariantach genetycznych.
    -   **BED**: Proste anotacje regionów genomowych.
    -   **GFF/GTF**: Szczegółowe anotacje genów i innych elementów.
-   **Znajomość tych formatów umożliwia efektywne wykorzystanie narzędzi bioinformatycznych** i interpretację wyników.

------------------------------------------------------------------------

**Wskazówki dla studentów:**

-   **Zapoznaj się z przykładowymi plikami w każdym formacie** aby lepiej zrozumieć strukturę danych.
-   **Przećwicz identyfikację i interpretację poszczególnych kolumn** na rzeczywistych danych.
-   **Pamiętaj o różnicach w numeracji pozycji** między formatami (np. 0-bazowa vs 1-bazowa).
-   **Zwracaj uwagę na informacje w polach opcjonalnych**, które mogą dostarczać istotnych danych do analizy.

------------------------------------------------------------------------

**Dodatkowe materiały:**

-   **Dokumentacja formatów:**
    -   FASTQ: [FASTQ Format Specification](https://en.wikipedia.org/wiki/FASTQ_format)
    -   SAM/BAM: [SAM Format Specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
    -   VCF: [VCF Format Specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
    -   BED: [BED Format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
    -   GFF/GTF: [GFF Format Specification](https://www.ensembl.org/info/website/upload/gff.html)
-   **Narzędzia do pracy z formatami:**
    -   **samtools**: Praca z plikami SAM/BAM.
    -   **bcftools**: Analiza plików VCF.
    -   **bedtools**: Operacje na danych w formacie BED.

------------------------------------------------------------------------

#### **Analiza WGS krok po kroku**

**1. Od surowych odczytów do zidentyfikowanych wariantów**

-   **Etapy:**
    1.  **Kontrola jakości**: Przycinanie, filtrowanie odczytów (FASTQ).
    2.  **Mapowanie**: Dopasowanie odczytów do genomu referencyjnego (SAM/BAM).
    3.  **Wykrywanie wariantów**: Identyfikacja SNP, indeli (VCF).
    4.  **Anotacja wariantów**: Przypisanie funkcji biologicznych.

**2. Rola każdego formatu pliku w procesie analizy**

-   **FASTQ**: Surowe dane od sekwenatora.
-   **SAM/BAM**: Zmapowane odczyty, baza do wykrywania wariantów.
-   **VCF**: Wynik analizy wariantów, dane do interpretacji biologicznej.
-   **BED/GFF/GTF**: Anotacje genomowe niezbędne do interpretacji funkcji.

------------------------------------------------------------------------

### **2. Instalacja i konfiguracja BioConductor**

#### **Kroki instalacji**

**1. Instalacja pakietu `BiocManager`**

-   **Cel**: Zarządzanie instalacją pakietów BioConductor.

-   **Polecenie:**

    ``` r
    install.packages("BiocManager")
    ```

-   **Ładowanie pakietu:**

    ``` r
    library(BiocManager)
    ```

**2. Instalacja pakietów do analizy WGS**

-   **Polecenie:**

    ``` r
    BiocManager::install(c("Rsamtools", "GenomicRanges", "GenomicAlignments", "VariantAnnotation"))
    ```

-   **Opis pakietów:**

    -   **`Rsamtools`**: Praca z plikami BAM/SAM.
    -   **`GenomicRanges`**: Reprezentacja regionów genomowych.
    -   **`GenomicAlignments`**: Analiza zmapowanych odczytów.
    -   **`VariantAnnotation`**: Anotacja i analiza wariantów.

**3. Aktualizacja pakietów**

-   **Polecenie:**

    ``` r
    BiocManager::install(version = "devel") # ostrożnie - instaluje deweloperską wersję pakietu, która może nie być stabilna
    ```

-   **Cel**: Instalacja wersji rozwojowej BioConductor (opcjonalnie).

**4. Rozwiązywanie problemów**

-   **Typowe błędy:**
    -   **Niekompatybilna wersja R**: Sprawdź wymagania pakietów.
    -   **Brakujące zależności**: Upewnij się, że wszystkie wymagane pakiety są zainstalowane.
-   **Rozwiązania:**
    -   **Aktualizacja R**: Pobierz najnowszą wersję ze strony CRAN.
    -   **Konsultacja dokumentacji**: Sprawdź instrukcje instalacji pakietu.

**5. Weryfikacja instalacji**

-   **Ładowanie pakietów:**

    ``` r
    library(Rsamtools)
    library(GenomicRanges)
    library(GenomicAlignments)
    library(VariantAnnotation)
    ```

-   **Sprawdzenie wersji:**

    ``` r
    sessionInfo()
    ```

----
