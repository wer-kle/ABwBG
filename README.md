# Analizy bioinformatyczne w badaniach genomowych

**Forma zajęć:** Ćwiczenia (33 godziny)\
**Prowadzący:** Weronika Klecel | [weronika\_klecel@sggw.edu.pl](mailto\:weronika_klecel@sggw.edu.pl) | pok. 25\
**Sala:** 14 (do 5 listopada), 1011 (od 12 listopada) \
**Terminy zajęć:** Wtorki w godzinach 10:00-13:00 \
**Konsultacje:** Do indywidualnego ustalenia 

## Opis przedmiotu

Przedmiot ma na celu zapoznanie studentów z metodami i narzędziami bioinformatycznymi wykorzystywanymi w analizie danych genomowych. Studenci zdobędą praktyczne umiejętności programowania w języku R (ze szczególnym uwzględnieniem projektu BioConductor), obsługi RStudio oraz nauczą się analizować dane z sekwencjonowania całogenomowego (WGS).

## Efekty kształcenia

Po ukończeniu przedmiotu student będzie potrafił:

- Korzystać z języka R i środowiska RStudio&#x20;
- Przeprowadzać podstawowe analizy danych WGS.
- Dokonywać wizualizacji wyników analiz WGS.
- Interpretować wyniki analiz bioinformatycznych w kontekście badań genomowych.

## Plan ćwiczeń

### **Zajęcia 1-3 (8 października – 22 października, wtorki 10:00-13:00)**

- **Temat:** Wprowadzenie do R i RStudio.
- **Zakres i cele:**
  - **Zapoznanie się z R i RStudio**: Omówienie środowiska RStudio (konsola, edytor skryptów, przeglądarka obiektów).
    - **Ćwiczenie**: Instalacja R i RStudio na komputerach uczestników.
  - **Podstawy języka R**: Wprowadzenie do składni języka, typów danych (wektory, ramki danych, listy), podstawowych operacji (tworzenie obiektów, operacje arytmetyczne).
    - **Ćwiczenie**: Proste skrypty do manipulacji danych i wykonywania obliczeń.

### **Zajęcia 4-7 (29 października – 19 listopada, wtorki 10:00-13:00)**
*Uwaga: od 12 listopada zajęcia w sali 1011 (1. piętro)*

- **Temat:** Wprowadzenie do analizy danych WGS.
- **Zakres i cele:**
  - **Preprocessing danych WGS**: Wprowadzenie do formatów danych WGS (FASTQ, BAM, VCF) oraz wstępnego przetwarzania danych (filtrowanie, kontrola jakości).
    - **Ćwiczenie**: Korzystanie z pakietu `ShortRead` do analizy jakości surowych danych sekwencjonowania.
  - **Mapowanie sekwencji i detekcja wariantów**: Omówienie procesu mapowania sekwencji na genom referencyjny, wykrywania wariantów (SNPs, indels) i tworzenia plików VCF.
    - **Ćwiczenie**: Użycie pakietu `Rsamtools` do przetwarzania i analizy danych BAM oraz `VariantAnnotation` do pracy z plikami VCF.
  - **Anotacja i interpretacja wariantów**: Omówienie procesów anotacji wariantów i ich interpretacji w kontekście genomu.
    - **Ćwiczenie**: Korzystanie z pakietu `AnnotationHub` do pobierania danych referencyjnych i `GenomicFeatures` do anotacji genomu.

### **Zajęcia 8-11 (26 listopada – 17 grudnia, wtorki 10:00-14:00)**

- **Temat:** Zaawansowane techniki analizy WGS i case studies.
- **Zakres i cele:**
  - **Analiza strukturalnych wariantów genomowych**: Omówienie rodzajów wariantów strukturalnych (np. duplikacje, delecje, translokacje) oraz narzędzi do ich analizy.
    - **Ćwiczenie**: Praca z pakietem `StructuralVariantAnnotation` do detekcji i analizy wariantów strukturalnych.
  - **Zaawansowane wizualizacje danych WGS**: Tworzenie złożonych wizualizacji wariantów, struktury genomu i wyników analizy.
    - **Ćwiczenie**: Użycie pakietów `ggbio` i `circlize` do tworzenia wykresów genomowych.
  - **Omówienie wybranych badań naukowych**: Przegląd publikacji wykorzystujących analizę WGS, dyskusja nad metodologią, wynikami i wnioskami.
    - **Ćwiczenie**: Analiza przypadków w grupach, prezentacja wybranych badań i dyskusja nad zastosowaniem poznanych technik analizy.

## Metody dydaktyczne

- Wykłady z prezentacjami multimedialnymi.
- Ćwiczenia praktyczne na komputerach z wykorzystaniem R i RStudio.
- Dyskusje i analizy case studies.

## Wymagania wstępne

- Podstawowa wiedza z zakresu genetyki i biologii molekularnej.
- Chęć nauki programowania w języku R.

## Literatura

- **Podstawowa:**

  - "Bioconductor Case Studies" – Hahne F., et al.
  - "R Programming for Bioinformatics" – Gentleman R.
  - "An Introduction to R" - Douglas A., et al.
  - "Primer to Analysis of Genomic Data Using R" - Gondro C.
  - "Computational Genomics with R" - Akalin A. <https://compgenomr.github.io/book>

- **Uzupełniająca:**

  - artykuły naukowe udostępniane przez prowadzących

## Warunki zaliczenia ćwiczeń

- Uczestnictwo w zajęciach (dopuszczalne dwie nieobecności).
- Zamieszczanie rozwiązań zadań z zajęć w repozytorium GitHub.
- Wykonanie dwóch projektów zaliczeniowych podczas zajęć 8. (26 listopada) i 15. (21 stycznia).

## Uwagi

- Wszystkie ćwiczenia będą realizowane na komputerach z systemem Windows w środowisku RStudio. Jeżeli preferujesz srodowisko Mac/Linux, możesz pracowac na wlasnym laptopie.

- Terminy, godziny zajęć oraz szczegółowa tematyka mogą ulec zmianie – aktualny harmonogram będzie aktualizowany na bieżąco.
