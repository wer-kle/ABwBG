# Zadanie 0: Przygotowanie środowiska
- Zainstaluj pakiety `VariantAnnotation`, `GenomicRanges`, `AnnotationHub`.
- Załaduj wymagane pakiety.

# Zadanie 1: Wczytanie i eksploracja danych
- Znajdź ścieżkę do przykładowego pliku VCF z pakietu `VariantAnnotation`.
- Wczytaj dane do obiektu typu `VCF`.
- Wyświetl podstawowe informacje o zawartości pliku:
  - Wyświetl nagłówki i omów strukturę metadanych.
  - Sprawdź, jakie kolumny z informacjami o wariantach (INFO, FORMAT) znajdują się w pliku.
  - Policz łączną liczbę wariantów.

# Zadanie 2: Filtracja i analiza jakości
- Sprawdź kolumnę `QUAL` w obiekcie `vcf`.
- Ustal kryterium jakości i odfiltruj warianty.

# Zadanie 3: Anotacja wariantów
- Sprawdź dostępność kolumn z informacjami w `info(vcf_filtered)`.
- Wykorzystaj pakiety anotacyjne, np. `TxDb`, do wzbogacenia wariantów o informacje genowe.
- Przeprowadź anotację wariantów za pomocą funkcji `locateVariants`.

# Zadanie 4: Przykłady dalszej analizy
- Wybierz warianty z regionów UTR (5'UTR lub 3'UTR) i porównaj ich liczebność.
- Wyodrębnij warianty znajdujące się wyłącznie w regionach międzygenowych.

# Zadanie 5: Podsumowanie
- Przygotuj krótkie podsumowanie wyników uzyskanych w powyższych zadaniach:
  - Liczba wariantów przed i po filtracji.
  - Liczba wariantów w regionach kodujących vs. międzygenowych.
- Omów możliwe ścieżki rozszerzenia analizy.
