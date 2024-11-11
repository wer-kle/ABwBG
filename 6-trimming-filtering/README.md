# Przycinanie (trimming) i filtrowanie danych NGS

## Przycinanie odczytów na podstawie jakości

Przycinanie odczytów na podstawie jakości jest kluczowym etapem w przetwarzaniu danych sekwencjonowania następnej generacji (NGS). Wartości jakości (np. Phred score) przypisywane do każdej bazy informują o wiarygodności odczytu tej bazy przez sekwenator. Bazy na końcach odczytów często charakteryzują się niższą jakością ze względu na zmniejszoną wydajność enzymów lub problemy z detekcją sygnału. Usunięcie tych niskiej jakości baz poprzez przycinanie zwiększa dokładność późniejszych analiz, takich jak mapowanie odczytów do genomu referencyjnego czy składanie de novo.

### **Dlaczego jakość baz jest ważna?**

- **Wysoka jakość baz**: Wskazuje na dużą pewność co do poprawności zidentyfikowanego nukleotydu.
- **Niska jakość baz**: Oznacza większe prawdopodobieństwo błędnej identyfikacji nukleotydu, co może prowadzić do błędów w dalszych analizach, takich jak mapowanie odczytów czy wykrywanie wariantów.

---

### **Parametry funkcji `trimTailw`**

Funkcja `trimTailw` z pakietu `ShortRead` służy do przycinania odczytów na podstawie jakości baz. Główne parametry tej funkcji to:

- **`k`**: Liczba kolejnych baz o jakości poniżej progu, która jest wymagana do rozpoczęcia przycinania.
- **`a`**: Symbol jakości (jeden znak), który określa próg jakości baz. Bazy o jakości niższej niż ten symbol są uznawane za niskiej jakości.
- **`halfwidth`**: Połowa szerokości okna średniej ruchomej (w bazach) używana do wygładzania wartości jakości.

Przyjrzyjmy się każdemu z tych parametrów bardziej szczegółowo.

---

#### **1. Symbol jakości (`a`)**

##### **Co to jest symbol jakości?**

- **Symbol jakości (`a`)** to pojedynczy znak (litera lub symbol), który odpowiada konkretnej wartości jakości w skali Phred.
- W plikach FASTQ wartości jakości są kodowane jako znaki ASCII, gdzie każdy znak reprezentuje wartość jakości dla danej bazy.

##### **Skala Phred**

- **Wartości jakości Phred** (Phred Quality Score) są numerycznymi wartościami, które reprezentują prawdopodobieństwo błędnego odczytu danego nukleotydu.
- Wartość jakości Phred jest obliczana jako:

  \[
  Q = -10 \times \log_{10}(P)
  \]

  gdzie \( P \) to prawdopodobieństwo błędu.

- **Wyższa wartość Phred** oznacza niższe prawdopodobieństwo błędu (lepszą jakość).

##### **Kodowanie wartości jakości**

- **ASCII**: Wartości Phred są mapowane na znaki ASCII w plikach FASTQ.
- Istnieją dwa główne standardy kodowania:
  - **Phred+33**: Używane w nowszych wersjach Illumina. Wartości Phred od 0 do 93 są mapowane na znaki ASCII od 33 do 126.
  - **Phred+64**: Używane w starszych wersjach Illumina.

##### **Przykładowa mapa znaków jakości (Phred+33)**

| Wartość Phred (Q) | Symbol jakości (ASCII) | Kod ASCII |
|-------------------|------------------------|-----------|
| 2                 | `"`                    | 34        |
| 10                | `+`                    | 43        |
| 20                | `5`                    | 53        |
| 30                | `?`                    | 63        |
| 40                | `I`                    | 73        |

##### **Przykład użycia symbolu jakości w `trimTailw`**

- Jeśli chcesz przyciąć bazy o jakości **Phred < 20**, musisz znaleźć odpowiedni symbol jakości.
- Dla **Phred+33**, wartość Phred 20 odpowiada znakowi ASCII 53, czyli symbolowi `5`.

- W funkcji `trimTailw` ustawiasz parametr `a = "5"`.

##### **Jak ustalić symbol jakości dla konkretnej wartości Phred?**

1. **Oblicz kod ASCII**:

   \[
   \text{Kod ASCII} = \text{Wartość Phred} + 33
   \]

   (dla kodowania Phred+33)

2. **Zamień kod ASCII na symbol**:

   - W R możesz użyć funkcji `intToUtf8()`:

     ```R
     symbol <- intToUtf8(Phred_value + 33)
     ```

   - **Przykład**: Dla wartości Phred 20:

     ```R
     symbol <- intToUtf8(20 + 33)  # Wynik: "5"
     ```

3. **Ustaw parametr `a`**:

   - W funkcji `trimTailw`, ustaw `a = symbol`.

---

#### **2. Liczba kolejnych baz o niskiej jakości (`k`)**

- **`k`** to minimalna liczba kolejnych baz o jakości poniżej progu (`a`), która musi zostać znaleziona, aby rozpocząć przycinanie.

- **Dlaczego `k` jest ważne?**

  - Pozwala uniknąć przycinania odczytu na podstawie pojedynczych baz o niskiej jakości, które mogą być przypadkowe.
  - Ustawienie `k` na wartość większą niż 1 (np. 2 lub 3) sprawia, że przycinanie jest bardziej konserwatywne i następuje tylko wtedy, gdy występuje ciąg baz o niskiej jakości.

- **Przykład**:

  - Jeśli ustawisz `k = 2`, to funkcja będzie szukać co najmniej dwóch kolejnych baz o jakości poniżej progu `a`, aby rozpocząć przycinanie.

---

#### **3. Szerokość okna średniej ruchomej (`halfwidth`)**

- **`halfwidth`** to połowa szerokości okna używanego do wygładzania wartości jakości za pomocą średniej ruchomej.

- **Jak działa okno średniej ruchomej?**

  - Dla każdej pozycji w odczycie funkcja oblicza średnią wartość jakości w oknie o szerokości \( 2 \times \text{halfwidth} + 1 \).
  - **halfwidth = 1** oznacza, że okno obejmuje 3 bazy (1 baza przed, bieżąca baza i 1 baza po).

- **Dlaczego używamy średniej ruchomej?**

  - Wartości jakości mogą być nieprawidłowe, gdy obliczamy je na poziomie pojedynczych baz.
  - Wygładzanie za pomocą średniej ruchomej pozwala na bardziej stabilną ocenę jakości w danym regionie odczytu.

- **Wpływ `halfwidth` na przycinanie:**

  - **Mniejsze wartości `halfwidth`** (np. 1) sprawiają, że funkcja jest bardziej czuła na lokalne spadki jakości.
  - **Większe wartości `halfwidth`** (np. 5) powodują, że funkcja jest mniej czuła na krótkotrwałe spadki jakości i koncentruje się na ogólnym trendzie.

- **Przykład:**

  - Ustawienie `halfwidth = 2` oznacza okno o szerokości 5 baz.

---

### **Przykład zastosowania parametrów w funkcji `trimTailw`**

Załóżmy, że chcesz przyciąć odczyty, usuwając końcowe bazy, gdy znajdziesz co najmniej **2 kolejne bazy o jakości Phred poniżej 20**, używając średniej ruchomej z oknem o szerokości 3 baz.

1. **Ustal symbol jakości dla Phred 20:**

   ```R
   symbol <- intToUtf8(20 + 33)  # Wynik: "5"
   ```

2. **Ustaw parametry w funkcji `trimTailw`:**

   ```R
   trimmed_reads <- trimTailw(fq_reads, k = 2, a = "5", halfwidth = 1)
   ```

   - **`k = 2`**: Szukaj 2 kolejnych baz o jakości poniżej "5" (Phred 20).
   - **`a = "5"`**: Baza o jakości Phred 20.
   - **`halfwidth = 1`**: Okno średniej ruchomej obejmujące 3 bazy.

---

### **Wskazówki praktyczne**

#### **1. Ustalanie progu jakości**

- **Typowe wartości progu jakości to Phred 20 lub Phred 30**, co odpowiada prawdopodobieństwu błędu odpowiednio 1% i 0,1%.
- Dostosuj próg jakości w zależności od wymagań analizy i jakości danych.

#### **2. Wybór parametru `k`**

- Jeśli dane są ogólnie wysokiej jakości, **mniejsze wartości `k`** mogą być wystarczające.
- Jeśli chcesz być bardziej konserwatywny i uniknąć nadmiernego przycinania, **większe wartości `k`** mogą być lepsze.

#### **3. Dostosowanie `halfwidth`**

- **Mniejsze wartości `halfwidth`** są przydatne, gdy jakość baz szybko się zmienia.
- **Większe wartości `halfwidth`** są odpowiednie, gdy chcesz wygładzić lokalne fluktuacje i skupić się na ogólnym trendzie.

---

### **Przykładowy kod z objaśnieniami**

```R
library(ShortRead)

# Wczytaj plik FASTQ
fq_reads <- readFastq("ścieżka/do/pliku.fastq")

# Ustal wartość Phred dla progu jakości (np. Phred 20)
phred_threshold <- 20

# Oblicz odpowiadający symbol jakości dla Phred+33
quality_symbol <- intToUtf8(phred_threshold + 33)

# Przytnij odczyty
trimmed_reads <- trimTailw(
  fq_reads,
  k = 2,               # Liczba kolejnych baz o niskiej jakości
  a = quality_symbol,  # Symbol jakości odpowiadający Phred 20
  halfwidth = 1        # Połowa szerokości okna średniej ruchomej
)

# Sprawdź liczbę odczytów przed i po przycinaniu
length(fq_reads)        # Przed przycinaniem
length(trimmed_reads)   # Po przycinaniu

# Opcjonalnie: Zapisz przycięte odczyty do nowego pliku FASTQ
writeFastq(trimmed_reads, "ścieżka/do/przyciętego_pliku.fastq")
```

---

### **Podsumowanie**

- **Symbol jakości (`a`)**: Odpowiada wartości Phred, która jest progiem jakości baz do przycinania. Użyj funkcji `intToUtf8()` do konwersji wartości Phred na symbol jakości.

- **Liczba kolejnych baz o niskiej jakości (`k`)**: Określa, ile baz o jakości poniżej progu musi wystąpić kolejno, aby rozpocząć przycinanie.

- **Szerokość okna średniej ruchomej (`halfwidth`)**: Połowa szerokości okna używanego do wygładzania wartości jakości. Pozwala na stabilniejszą ocenę jakości baz.

- **Dostosowanie parametrów**: Parametry te można dostosować w zależności od jakości danych i specyfiki analizy. Ważne jest znalezienie balansu między usunięciem niskiej jakości baz a zachowaniem jak największej ilości wartościowych danych.

---

## Filtracja odczytów o niskiej jakości

Filtracja odczytów o niskiej jakości polega na usunięciu całych odczytów, które nie spełniają określonych kryteriów jakościowych, takich jak minimalna długość czy średnia wartość jakości. Odczyty te mogą wprowadzać szumy i błędy w analizach bioinformatycznych, prowadząc do fałszywych pozytywów lub zaniżenia pokrycia genomu. Poprzez zastosowanie filtracji, uzyskuje się bardziej wiarygodny zestaw danych do dalszych analiz.

---

## Ponowna kontrola jakości po przycinaniu i filtracji

Ponowna kontrola jakości po przycinaniu i filtracji jest niezbędna do oceny skuteczności tych procesów. Analiza danych po przetwarzaniu pozwala na sprawdzenie, czy jakość odczytów uległa poprawie oraz czy zostały usunięte potencjalne artefakty. Porównanie raportów jakości przed i po tych operacjach umożliwia identyfikację ewentualnych problemów i podjęcie dalszych działań optymalizacyjnych.

---

## Analiza rozkładu długości odczytów

Analiza rozkładu długości odczytów dostarcza informacji o wpływie przycinania na dane sekwencyjne. Przycinanie może prowadzić do skrócenia odczytów, co z kolei może wpływać na efektywność mapowania czy składania genomu. Porównanie histogramów długości przed i po przycinaniu pozwala na ocenę, czy odczyty są wystarczająco długie do planowanych analiz oraz czy nie doszło do utraty istotnej ilości danych.

---

## Wykrywanie i usuwanie sekwencji adapterów

Sekwencje adapterów są sztucznymi fragmentami DNA dodawanymi podczas przygotowywania bibliotek NGS, niezbędnymi do przyłączenia odczytów do flow cell. Czasami fragmenty DNA są krótsze niż długość odczytu, co skutkuje sekwencjonowaniem części adaptera. Obecność sekwencji adapterów w danych może zakłócać analizę, wprowadzając błędne dopasowania czy fałszywe warianty. Wykrywanie i usuwanie tych sekwencji jest zatem kluczowe dla poprawy jakości danych i wiarygodności wyników.

---

## Ponowna kontrola jakości po usunięciu adapterów

Po usunięciu sekwencji adapterów konieczne jest przeprowadzenie ponownej kontroli jakości w celu oceny wpływu tego procesu na dane. Usunięcie adapterów powinno skutkować zmniejszeniem liczby odczytów zawierających niepożądane sekwencje, co przekłada się na lepsze wyniki w analizach downstream. Porównanie raportów jakości przed i po tym kroku pozwala na potwierdzenie skuteczności usuwania adapterów oraz identyfikację ewentualnych problemów, takich jak nadmierne skracanie odczytów.

---

## Dokumentacja i raportowanie wyników

Dokumentacja i raportowanie przeprowadzonych operacji są kluczowe dla transparentności i replikowalności badań bioinformatycznych. Tworzenie raportów zawierających użyte metody, parametry oraz wyniki analiz pozwala na lepsze zrozumienie procesu przetwarzania danych i ułatwia współpracę między naukowcami. Wykorzystanie narzędzi takich jak R Markdown umożliwia integrację kodu, wyników i interpretacji w jednym dokumencie, co sprzyja efektywnemu komunikowaniu wniosków i zapewnia ścieżkę audytu dla przeprowadzonych analiz.

---

## Podsumowanie

Proces przycinania i filtracji danych NGS jest niezbędnym etapem przygotowania danych do dalszych analiz. Poprzez usunięcie niskiej jakości baz i sekwencji adapterów oraz odrzucenie odczytów niespełniających kryteriów jakościowych, zwiększa się dokładność i wiarygodność wyników bioinformatycznych. Regularna kontrola jakości na każdym etapie przetwarzania danych pozwala na optymalizację procesu i identyfikację potencjalnych problemów, co jest kluczowe dla sukcesu projektów sekwencjonowania.

---

**Dla chętnych:**

### Wykorzystanie narzędzia FastQC do analizy jakości danych i porównanie wyników z raportami z Bioconductora

FastQC jest powszechnie używanym narzędziem do szybkiej oceny jakości danych NGS. Dostarcza ono szeregu statystyk i wizualizacji, takich jak rozkład wartości jakości, zawartość GC czy obecność sekwencji adapterów. Porównanie wyników z FastQC z raportami generowanymi przez pakiety Bioconductora pozwala na uzyskanie pełniejszego obrazu jakości danych oraz na weryfikację spójności między różnymi metodami analizy.

### Przeprowadzenie zbiorczej analizy QC dla wielu próbek z użyciem MultiQC

MultiQC to narzędzie umożliwiające integrację wyników kontroli jakości z wielu próbek i różnych narzędzi w jeden, zbiorczy raport. Jest szczególnie przydatne w projektach obejmujących dużą liczbę próbek, gdzie manualne porównywanie indywidualnych raportów byłoby nieefektywne. MultiQC ułatwia identyfikację globalnych trendów i potencjalnych problemów w zbiorze danych, co pozwala na bardziej efektywne planowanie dalszych analiz.

---
