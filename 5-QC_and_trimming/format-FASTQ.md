Plik **Fastq** to format tekstowy służący do przechowywania sekwencji nukleotydowych wraz z odpowiadającymi im informacjami o jakości odczytu. Jest powszechnie używany w sekwencjonowaniu nowej generacji (NGS) do przechowywania dużych ilości danych sekwencyjnych.

**Struktura pliku Fastq składa się z powtarzających się rekordów, z których każdy zawiera cztery linie:**

1. **Linia 1**: Nagłówek sekwencji
   - Zaczyna się od znaku `@`.
   - Następnie zawiera **identyfikator sekwencji** i opcjonalnie dodatkowy opis.
   - **Przykład:**
     ```
     @SEQ_ID
     ```

2. **Linia 2**: Sekwencja nukleotydów
   - Zawiera ciąg liter reprezentujących nukleotydy: A, T, C, G i N (gdzie N oznacza nieokreślony nukleotyd).
   - **Przykład:**
     ```
     GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTTA
     ```

3. **Linia 3**: Separator
   - Zaczyna się od znaku `+`.
   - Może być pusta lub powtarzać identyfikator z linii 1 (chociaż nie jest to wymagane).
   - Służy jako separator między sekwencją a danymi jakości.
   - **Przykład:**
     ```
     +
     ```

4. **Linia 4**: Dane jakości odczytu
   - Zawiera symbole reprezentujące **jakość każdego nukleotydu** w sekwencji z linii 2.
   - Długość tej linii musi być taka sama jak długość sekwencji w linii 2.
   - Jakość jest zakodowana za pomocą znaków ASCII, zgodnie ze skalą jakości (np. skala Phred).
   - **Przykład:**
     ```
     !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
     ```

**Pełny przykład rekordu w pliku Fastq:**

```
@SEQ_ID
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTTA
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
```

**Dodatkowe informacje:**

- **Identyfikator sekwencji** w linii 1 często zawiera informacje takie jak numer odczytu, pozycja w chipie sekwencjonującym czy nazwa próbki.
- **Dane jakości** w linii 4 są kluczowe dla oceny wiarygodności każdego nukleotydu w sekwencji. Wyższy znak ASCII odpowiada wyższej jakości odczytu.
- Pliki Fastq mogą być bardzo duże ze względu na ilość danych generowanych podczas sekwencjonowania, dlatego często są kompresowane (np. do formatu `.fastq.gz`).

**Zastosowania:**

- Analiza sekwencji DNA i RNA.
- Filtracja i przycinanie danych na podstawie jakości.
- Mapowanie sekwencji do genomu referencyjnego.
- Wykrywanie wariantów genetycznych.
