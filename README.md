# Analizy bioinformatyczne w badaniach genomowych

**Forma zajęć:** Ćwiczenia (14 x 3 godziny = 42 godziny)\
**Prowadzący:** Weronika Klecel | [weronika\_klecel@sggw.edu.pl](mailto:weronika_klecel@sggw.edu.pl) | pok. 25\
**Sala:** 14 (do 5 listopada), 1011 (od 12 listopada)\
**Terminy zajęć:** Wtorki w godzinach 10:00-13:00\
**Konsultacje:** Do indywidualnego ustalenia

## Opis przedmiotu

Przedmiot ma na celu zapoznanie studentów z metodami i narzędziami bioinformatycznymi wykorzystywanymi w analizie danych genomowych. Studenci zdobędą praktyczne umiejętności programowania w języku R (ze szczególnym uwzględnieniem projektu Bioconductor), obsługi RStudio oraz nauczą się analizować dane z sekwencjonowania całogenomowego (WGS). Wszystkie analizy będą wykonywane z użyciem pakietów Bioconductora; dla chętnych dostępne będą również narzędzia wymagające znajomości bash.

## Efekty kształcenia

Po ukończeniu przedmiotu student będzie potrafił:

- Korzystać z języka R i środowiska RStudio.
- Przeprowadzać podstawowe i zaawansowane analizy danych WGS.
- Wykonywać kontrolę jakości oraz preprocessing danych sekwencjonowania.
- Wykrywać i anotować warianty genetyczne.
- Przeprowadzać analizy asocjacyjne (GWAS) i interpretować ich wyniki.
- Dokonywać wizualizacji wyników analiz genomowych.
- Interpretować wyniki analiz bioinformatycznych w kontekście badań genomowych.

## Plan ćwiczeń

### **Zajęcia 1-3 (8 października, 15 października, 22 października 2024)**

- **Temat:** Wprowadzenie do R i RStudio.
- **Zakres i cele:**
  - **Zapoznanie się z R i RStudio:** Omówienie środowiska RStudio (konsola, edytor skryptów, przeglądarka obiektów).
    - **Ćwiczenie:** Instalacja R i RStudio na komputerach studentów.
  - **Podstawy języka R:** Wprowadzenie do składni języka, typów danych (wektory, ramki danych, listy), podstawowych operacji (tworzenie obiektów, operacje arytmetyczne).
    - **Ćwiczenie:** Proste skrypty do manipulacji danych i wykonywania obliczeń.

### **Zajęcia 4 (29 października 2024)**

- **Temat:** Wprowadzenie do WGS i Bioconductor.
- **Zakres i cele:**
  - **Sekwencjonowanie całogenomowe (WGS):** Omówienie podstaw sekwencjonowania całogenomowego, formatów danych (FASTQ, BAM, VCF), pipeline'ów analizy.
  - **Bioconductor:** Wprowadzenie do projektu Bioconductor, instalacja pakietów, przegląd dostępnych narzędzi.
    - **Ćwiczenie:** Instalacja i konfiguracja pakietów Bioconductora.

### **Zajęcia 5 (5 listopada 2024)**

- **Temat:** Kontrola jakości (QC).
- **Zakres i cele:**
  - **Kontrola jakości surowych danych sekwencjonowania:** Omówienie metod i narzędzi do oceny jakości danych.
    - **Ćwiczenie:** Użycie pakietów `ShortRead` oraz `qcmetrics` do analizy jakości danych FASTQ.
  - **Wykrywanie i usuwanie błędów sekwencjonowania:** Identyfikacja niskiej jakości baz, adapterów, kontaminacji.
    - **Ćwiczenie:** Praktyczna analiza danych z wykorzystaniem Bioconductora.

### **Zajęcia 6 (12 listopada 2024)**
*Uwaga: od 12 listopada zajęcia w sali 1011 (1. piętro)*

- **Temat:** Przycinanie i filtracja.
- **Zakres i cele:**
  - **Przycinanie danych sekwencjonowania:** Omówienie technik przycinania niskiej jakości sekwencji i adapterów.
    - **Ćwiczenie:** Korzystanie z pakietów `ShortRead` lub narzędzi takich jak `TrimGalore` (dla chętnych z użyciem bash).
  - **Filtracja danych:** Usuwanie kontaminacji, filtrowanie według jakości i długości odczytów.
    - **Ćwiczenie:** Praktyczne zastosowanie filtracji z użyciem Bioconductora.

### **Zajęcia 7 (19 listopada 2024)**

- **Temat:** Mapowanie do genomu referencyjnego (alignment).
- **Zakres i cele:**
  - **Mapowanie odczytów:** Omówienie algorytmów mapowania (np. BWA, Bowtie2), znaczenie indeksowania genomu referencyjnego.
    - **Ćwiczenie:** Użycie pakietu `Rsubread` do mapowania odczytów (dla chętnych narzędzia z linii poleceń w bash).
  - **Ocena jakości mapowania:** Analiza statystyk mapowania, identyfikacja problemów.
    - **Ćwiczenie:** Analiza wyników mapowania z użyciem `Rsamtools`.

### **Zajęcia 8 (26 listopada 2024)**

- **Temat:** Kolokwium z zajęć 1-7.
- **Zakres i cele:**
  - **Sprawdzenie wiedzy i umiejętności:** Kolokwium praktyczne obejmujące materiał z zajęć 1-7.
    - **Ćwiczenie:** Rozwiązanie zadań praktycznych w R i Bioconductor.

### **Zajęcia 9 (3 grudnia 2024)**

- **Prowadzący:** dr Marlena Wojciechowska
- **Temat:** Analiza danych z mikromacierzy w programie GenomeStudio.
- **Zakres i cele:**
  - **Wprowadzenie do mikromacierzy:** Omówienie technologii mikromacierzowej, zastosowań w genomice.
  - **Program GenomeStudio:** Przegląd funkcjonalności, import i analiza danych.
    - **Ćwiczenie:** Praktyczna analiza danych mikromacierzowych w GenomeStudio.

### **Zajęcia 10-14 (10 grudnia 2024 – 21 stycznia 2025)**

- **Temat:** Analiza downstream (GWAS, Variant Calling itp.).
- **Zakres i cele:**
  - **Wykrywanie wariantów (Variant Calling):** Omówienie metod detekcji SNPs i indels.
    - **Ćwiczenie:** Użycie pakietu `VariantTools` lub `GATK` (dla chętnych z użyciem bash) do wykrywania wariantów.
  - **Analiza asocjacji genomowej (GWAS):** Wprowadzenie do GWAS, przygotowanie danych, analiza statystyczna.
    - **Ćwiczenie:** Praca z pakietami `snpStats` lub `GenABEL`.
  - **Anotacja wariantów:** Interpretacja biologiczna wykrytych wariantów.
    - **Ćwiczenie:** Użycie pakietów `VariantAnnotation` i `biomaRt`.
  - **Wizualizacja wyników:** Tworzenie wykresów typu Manhattan, Q-Q.
    - **Ćwiczenie:** Wykorzystanie pakietów `qqman` lub `ggplot2`.

*Terminy zajęć:*
- **Zajęcia 10:** 10 grudnia 2024
- **Zajęcia 11:** 17 grudnia 2024
- *Przerwa świąteczna*
- **Zajęcia 12:** 7 stycznia 2025
- **Zajęcia 13:** 14 stycznia 2025
- **Zajęcia 14:** 21 stycznia 2025

### **Zajęcia 14 (21 stycznia 2025)**

- **Temat:** Kolokwium z zajęć 9-14.
- **Zakres i cele:**
  - **Sprawdzenie wiedzy i umiejętności:** Kolokwium praktyczne obejmujące materiał z zajęć 9-14.
    - **Ćwiczenie:** Rozwiązanie zadań praktycznych z analizy danych mikromacierzowych i downstream.

## Metody dydaktyczne

- Wykłady z prezentacjami multimedialnymi.
- Ćwiczenia praktyczne na komputerach z wykorzystaniem R, RStudio i pakietów Bioconductora.
- Dyskusje i analizy case studies.
- Dla chętnych: dodatkowe zadania z użyciem narzędzi wymagających znajomości bash.

## Wymagania wstępne

- Podstawowa wiedza z zakresu genetyki i biologii molekularnej.
- Chęć nauki programowania w języku R.

## Literatura

- **Podstawowa:**

  - "Bioconductor Case Studies" – Hahne F., et al.
  - "R Programming for Bioinformatics" – Gentleman R.
  - "An Introduction to R" – Douglas A., et al.
  - "Primer to Analysis of Genomic Data Using R" – Gondro C.
  - "Computational Genomics with R" – Akalin A. <https://compgenomr.github.io/book>

- **Uzupełniająca:**

  - Artykuły naukowe udostępniane przez prowadzących.

## Warunki zaliczenia ćwiczeń

- Uczestnictwo w zajęciach (dopuszczalne dwie nieobecności).
- Zamieszczanie rozwiązań zadań z zajęć w repozytorium GitHub.
- Wykonanie dwóch kolokwiów zaliczeniowych podczas zajęć 8. (26 listopada) i 14. (21 stycznia).

## Uwagi

- Wszystkie ćwiczenia będą realizowane na komputerach z systemem Windows w środowisku RStudio. Jeżeli preferujesz środowisko Mac/Linux, możesz pracować na własnym laptopie.
- Terminy, godziny zajęć oraz szczegółowa tematyka mogą ulec zmianie – aktualny harmonogram będzie aktualizowany na bieżąco.
- Dla chętnych przewidziano dodatkowe zadania z użyciem narzędzi wymagających znajomości bash.
