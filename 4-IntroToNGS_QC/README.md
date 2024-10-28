# Bioinformatyczna analiza danych genomowych (NGS workflow)

---

## **Wprowadzenie**

Surowe dane z jakiejkolwiek platformy NGS nie są zbyt użyteczne. Aby przekuć je w wartościowe dane biologiczne, musimy wykonać szereg procedur bioinformatycznych.

---

## **Etapy analizy danych NGS**

1. **Kontrola jakości danych surowych (QC)**
2. **Przetwarzanie odczytów (przycinanie i filtracja)**
3. **Mapowanie odczytów do genomu referencyjnego**
4. **Przetwarzanie po mapowaniu (sortowanie, indeksowanie, oznaczanie duplikatów)**
5. **Detekcja wariantów (variant calling)**
6. **Filtracja i annotacja wariantów**
7. **Analiza downstream (interpretacja wyników)**

---

### **1. Kontrola jakości danych surowych (QC)**

#### **Opis i znaczenie**

- **Cel:** Ocena jakości surowych odczytów sekwencji DNA/RNA uzyskanych z sekwencjonowania NGS.
- **Znaczenie:** Pozwala na identyfikację potencjalnych problemów, takich jak:
  - Niska jakość baz na końcach odczytów.
  - Obecność sekwencji adapterów.
  - Zanieczyszczenia biologiczne lub techniczne.
  - Nierównomierna reprezentacja fragmentów DNA (np. GC bias).
- **Efekt:** Decyzja o konieczności przetwarzania danych przed dalszą analizą, co zwiększa wiarygodność i dokładność wyników.

#### **Przykłady oprogramowania**

- **Bioconductor:**
  - [ShortRead](https://github.com/Bioconductor/ShortRead): Wczytywanie i podstawowa analiza jakości odczytów.
  - [Rqc](https://bioconductor.org/packages/release/bioc/html/Rqc.html): Generowanie interaktywnych raportów jakości.
  - [qrqc](https://github.com/vsbuffalo/qrqc): Szybka kontrola jakości z wizualizacjami.

- **Inne popularne narzędzia:**
  - [FastQC](https://github.com/s-andrews/FastQC): Standardowe narzędzie do oceny jakości danych NGS z czytelnymi raportami.
  - [MultiQC](https://github.com/MultiQC/MultiQC): Łączy raporty z wielu narzędzi QC w jeden zbiorczy raport.
  - [FastQ Screen](https://github.com/StevenWingett/FastQ-Screen): Sprawdza, czy odczyty nie pochodzą z zanieczyszczeń innych gatunków.

---

### **2. Przetwarzanie odczytów (przycinanie i filtracja)**

#### **Opis i znaczenie**

- **Cel:** Usunięcie sekwencji adapterów, niskiej jakości baz oraz filtrowanie odczytów o zbyt krótkiej długości.
- **Znaczenie:** Poprawia jakość danych przed mapowaniem, co zwiększa dokładność wyrównania odczytów do genomu referencyjnego i redukuje liczbę błędów w dalszych analizach.
- **Efekt:** Zwiększenie wiarygodności detekcji wariantów i innych analiz downstream.

#### **Przykłady oprogramowania**

- **Bioconductor:**
  - **ShortRead:** Funkcje do przycinania i filtracji odczytów.
  - **Biostrings:** Operacje na sekwencjach DNA, w tym manipulacje i przycinanie.

- **Inne popularne narzędzia:**
  - **Trimmomatic:** Narzędzie do przycinania odczytów NGS z licznymi opcjami.
  - **Cutadapt:** Usuwa sekwencje adapterów i przycina odczyty na podstawie jakości.
  - **FASTX-Toolkit:** Zestaw narzędzi do filtrowania i przetwarzania odczytów FASTQ.

---

### **3. Mapowanie odczytów do genomu referencyjnego (alignment)**

#### **Opis i znaczenie**

- **Cel:** Wyrównanie przetworzonych odczytów do znanej sekwencji genomu referencyjnego w celu określenia ich położenia genomowego.
- **Znaczenie:** Kluczowy krok umożliwiający identyfikację wariantów genetycznych oraz analizę ekspresji genów.
- **Efekt:** Powstaje plik BAM/SAM zawierający informacje o wyrównaniu każdego odczytu do genomu.

#### **Przykłady oprogramowania**

- **Bioconductor:**
  - **Rsubread:** Pakiet do mapowania odczytów DNA i RNA.
  - **Rhisat2:** Interfejs R do narzędzia HISAT2.

- **Inne popularne narzędzia:**
  - **BWA (Burrows-Wheeler Aligner):** Szybki algorytm do mapowania krótkich odczytów.
  - **Bowtie2:** Wydajne mapowanie odczytów do dużych genomów.
  - **HISAT2:** Mapowanie odczytów RNA-seq z uwzględnieniem splicingów.

---

### **4. Przetwarzanie po mapowaniu (sortowanie, indeksowanie, oznaczanie duplikatów)**

#### **Opis i znaczenie**

- **Cel:** Dalsze przetwarzanie wyrównanych odczytów w celu poprawy jakości danych przed detekcją wariantów.
- **Znaczenie:** 
  - **Sortowanie:** Ułatwia dostęp do danych i jest wymagane przez niektóre narzędzia.
  - **Indeksowanie:** Przyspiesza wyszukiwanie w plikach BAM.
  - **Oznaczanie duplikatów PCR:** Usuwa powielone odczyty, które mogą prowadzić do fałszywych pozytywów.
- **Efekt:** Przygotowane i oczyszczone dane do dalszej analizy wariantów.

#### **Przykłady oprogramowania**

- **Bioconductor:**
  - **Rsamtools:** Operacje na plikach BAM/SAM, w tym sortowanie i indeksowanie.
  - **GenomicAlignments:** Analiza wyrównań genomowych.
  - **picardr:** Interfejs R do narzędzi Picard, w tym oznaczania duplikatów.

- **Inne popularne narzędzia:**
  - **SAMtools:** Zestaw narzędzi do manipulacji plikami SAM/BAM.
  - **Picard Tools:** Kolekcja narzędzi Java do przetwarzania plików BAM/SAM.
  - **GATK (Genome Analysis Toolkit):** Umożliwia przetwarzanie i analizę danych NGS.

---

### **5. Detekcja wariantów (variant calling)**

#### **Opis i znaczenie**

- **Cel:** Identyfikacja wariantów genetycznych (np. SNP-ów, indeli) poprzez analizę wyrównanych odczytów względem genomu referencyjnego.
- **Znaczenie:** Umożliwia odkrycie mutacji, które mogą mieć znaczenie biologiczne.
- **Efekt:** Uzyskanie listy potencjalnych wariantów genetycznych do dalszej analizy.

#### **Przykłady oprogramowania**

- **Bioconductor:**
  - **VariantTools:** Detekcja wariantów w środowisku R.
  - **VarScan:** Analiza wariantów na poziomie populacji.

- **Inne popularne narzędzia:**
  - **GATK HaplotypeCaller:** Narzędzie standardowe do detekcji wariantów jednonukleotydowych i indeli.
  - **FreeBayes:** Wykrywanie wariantów u jednego lub wielu osobników.
  - **SAMtools mpileup/bcftools:** Detekcja wariantów na podstawie stosów odczytów.

---

### **6. Filtracja i annotacja wariantów**

#### **Opis i znaczenie**

- **Cel:** Usunięcie fałszywie pozytywnych wariantów oraz przypisanie informacji biologicznej do zidentyfikowanych wariantów.
- **Znaczenie:** 
  - **Filtracja:** Zwiększa wiarygodność wyników poprzez eliminację wariantów o niskiej jakości lub małym znaczeniu.
  - **Annotacja:** Dostarcza kontekstu biologicznego, np. wpływ na geny, konsekwencje funkcjonalne.
- **Efekt:** Lista wysokiej jakości wariantów z informacją o ich potencjalnym wpływie.

#### **Przykłady oprogramowania**

- **Bioconductor:**
  - **VariantAnnotation:** Analiza i annotacja plików VCF w R.
  - **Ensembl VEP (Variant Effect Predictor):** Dostępny również jako pakiet w R do annotacji wariantów.

- **Inne popularne narzędzia:**
  - **ANNOVAR:** Annotacja wariantów z wykorzystaniem różnych baz danych.
  - **SnpEff:** Szybka annotacja i predykcja efektu wariantów.
  - **Ensembl VEP:** Narzędzie online i do pobrania do annotacji wariantów.

---

### **7. Analiza downstream (interpretacja wyników)**

#### **Opis i znaczenie**

- **Cel:** Interpretacja biologiczna i kliniczna zidentyfikowanych wariantów, identyfikacja genów kandydujących, analiza funkcjonalna.
- **Znaczenie:** 
  - Pozwala na zrozumienie, jak wykryte warianty wpływają na funkcje biologiczne.
  - Pomaga w identyfikacji potencjalnych celów terapeutycznych lub markerów diagnostycznych.
- **Efekt:** Wnioski biologiczne lub kliniczne, które mogą prowadzić do publikacji naukowych lub zastosowań medycznych.

#### **Przykłady oprogramowania**

- **Bioconductor:**
  - **clusterProfiler:** Analiza wzbogacenia funkcjonalnego i ścieżek.
  - **org.Hs.eg.db:** Anotacje genów dla Homo sapiens.
  - **GenomicRanges:** Praca z zakresami genomowymi i integracja danych.

- **Inne popularne narzędzia:**
  - **DAVID:** Narzędzie online do analizy funkcjonalnej genów.
  - **Reactome Pathway Database:** Analiza ścieżek sygnałowych.
  - **Cytoscape:** Wizualizacja i analiza sieci biologicznych.

---

## **Podsumowanie**

Analiza danych NGS jest procesem wieloetapowym, który wymaga starannego przetwarzania i interpretacji danych. Każdy etap ma kluczowe znaczenie dla uzyskania wiarygodnych i użytecznych wyników. Wykorzystanie odpowiednich narzędzi bioinformatycznych, zarówno z pakietów Bioconductor, jak i innych popularnych programów, pozwala na efektywną analizę dużych zbiorów danych sekwencyjnych.

---

## **Dodatkowe informacje o oprogramowaniu**

### **Bioconductor**

- **Bioconductor** to projekt open-source, który dostarcza narzędzi do analizy danych biologicznych w języku R.
- **Zalety:**
  - Integracja z R pozwala na wykorzystanie zaawansowanych metod statystycznych i wizualizacji.
  - Umożliwia reprodukowalność analiz poprzez skrypty R.
  - Aktywna społeczność i regularne aktualizacje.

### **Inne popularne narzędzia**

- **FastQC:** Szybkie i łatwe w użyciu narzędzie do kontroli jakości, z interfejsem graficznym.
- **Trimmomatic i Cutadapt:** Popularne narzędzia do przycinania odczytów, często używane w pipeline'ach NGS.
- **BWA i Bowtie2:** Standardy w mapowaniu odczytów do genomu referencyjnego.
- **GATK:** Kompleksowe narzędzie do przetwarzania danych NGS, w tym detekcji wariantów i ich filtracji.
- **ANNOVAR i SnpEff:** Umożliwiają szybką i efektywną annotację wariantów z wykorzystaniem różnych baz danych.

---

## **Praktyczne wskazówki**

- **Planowanie analizy:** Zanim rozpoczniesz analizę, zrozum cele badania i wybierz odpowiednie narzędzia dla każdego etapu.
- **Kontrola jakości:** Nigdy nie pomijaj kontroli jakości danych surowych. Wczesne wykrycie problemów oszczędza czas i zasoby.
- **Dokumentacja:** Prowadź dokładną dokumentację swoich analiz, w tym wersji użytych narzędzi i parametrów.
- **Reprodukowalność:** Używaj skryptów i pipeline'ów, które pozwalają na odtworzenie analizy przez innych badaczy.
- **Aktualizacja oprogramowania:** Regularnie sprawdzaj aktualizacje narzędzi, aby korzystać z najnowszych funkcji i poprawek błędów.

---

## **Dodatkowe materiały i lektury**

- **Bioconductor:** [https://www.bioconductor.org/](https://www.bioconductor.org/)
- **Galaxy Project:** Platforma online do analizy danych NGS bez konieczności instalacji oprogramowania. [https://usegalaxy.org/](https://usegalaxy.org/)
- **GATK Best Practices:** Wytyczne dotyczące analizy danych NGS z użyciem GATK. [https://gatk.broadinstitute.org/](https://gatk.broadinstitute.org/)

---

## **Pytania do samodzielnego przemyślenia**

1. **Dlaczego ważne jest usuwanie duplikatów PCR podczas analizy danych NGS? Jakie mogą być konsekwencje ich pozostawienia?**
2. **Porównaj zalety i wady korzystania z pakietów Bioconductor w R z innymi narzędziami bioinformatycznymi. W jakich sytuacjach jedno podejście może być lepsze od drugiego?**
3. **Jakie czynniki należy wziąć pod uwagę przy wyborze genomu referencyjnego do mapowania odczytów?**
4. **W jaki sposób annotacja wariantów przyczynia się do zrozumienia ich potencjalnego wpływu na organizm?**


