# Kodowanie ASCII

Kodowanie ASCII (*American Standard Code for Information Interchange*) to standard kodowania znaków, który przypisuje numeryczne wartości liczbowe literom, cyfrom, znakom interpunkcyjnym oraz innym symbolom używanym w języku angielskim i innych językach opartych na alfabecie łacińskim.

---

## **Podstawowe informacje**

- **Zakres kodów ASCII:**
  - ASCII wykorzystuje **7-bitowe** liczby całkowite do reprezentowania **128 unikalnych znaków** (od 0 do 127).
- **Znaki kontrolne (0-31):**
  - Pierwsze 32 kody są zarezerwowane dla znaków kontrolnych, takich jak:
    - **LF** (Line Feed) - nowa linia.
    - **CR** (Carriage Return) - powrót karetki.
    - **TAB** - tabulacja.
- **Znaki drukowalne (32-126):**
  - Obejmują:
    - **Litery wielkie:** A-Z (kody 65-90).
    - **Litery małe:** a-z (kody 97-122).
    - **Cyfry:** 0-9 (kody 48-57).
    - **Znaki interpunkcyjne i specjalne:** np. !, @, #, $, %, ^, &, *, (, ), -, _, =, +, {, }, [, ], |, \, :, ;, ", ', <, >, ?, /, ~, `.
- **Kod 127:**
  - **DEL (Delete):** znak usuwania.

---

## **Zastosowanie kodowania ASCII w bioinformatyce**

W kontekście **bioinformatyki**, a szczególnie przy analizie danych **WGS (Whole Genome Sequencing)**, kodowanie ASCII jest istotne przy reprezentacji **jakości odczytów sekwencji** w plikach **FASTQ**.

---

#### **Kodowanie jakości w plikach FASTQ**

1. **Skala Phred:**
   - Jakość każdego nukleotydu jest wyrażana w **skali Phred**.
   - **Wzór:**
     \[
     Q = -10 \times \log_{10}(P)
     \]
     gdzie **Q** to wartość jakości, a **P** to prawdopodobieństwo błędnego odczytu nukleotydu.
   - **Interpretacja:**
     - Wyższa wartość **Q** oznacza **wyższą jakość** i **niższe prawdopodobieństwo błędu**.
     - Na przykład, Q=20 odpowiada prawdopodobieństwu błędu 1 na 100 (P=0.01).

2. **Przekształcenie wartości jakości na znaki ASCII:**
   - Wartość **Q** jest przekształcana na odpowiadający jej **znak ASCII** poprzez dodanie **offsetu**.
   - **Offset zależy od platformy sekwencjonującej:**
     - **Sanger / Illumina 1.8+**: Offset 33 (phred+33).
     - **Stare wersje Illumina (1.3-1.7):** Offset 64 (phred+64).

3. **Proces kodowania:**
   - **Formuła:**
     \[
     \text{Kod ASCII} = Q + \text{Offset}
     \]
   - **Przykład:**
     - Jeśli **Q = 30** i używamy offsetu 33:
       \[
       \text{Kod ASCII} = 30 + 33 = 63
       \]
     - Kod ASCII 63 odpowiada znakowi **'?'**.

4. **Linia jakości w pliku FASTQ:**
   - Każdy znak w linii jakości odpowiada jakości jednego nukleotydu w sekwencji.
   - **Przykładowa linia:**
     ```
     !''*((((***+))%%%++)(%%%%).1***-+*''))**
     ```
     - Znaki w tej linii reprezentują wartości jakości odczytów.

---

#### **Przykład obliczeń wartości jakości**

1. **Znak 'I':**
   - **Kod ASCII 'I'**: 73.
   - **Offset 33:** Przyjmujemy phred+33.
   - **Obliczenie Q:**
     \[
     Q = 73 - 33 = 40
     \]
   - **Prawdopodobieństwo błędu:**
     \[
     P = 10^{-Q/10} = 10^{-4} = 0.0001
     \]

2. **Znak '#':**
   - **Kod ASCII '#':** 35.
   - **Q:**
     \[
     Q = 35 - 33 = 2
     \]
   - **Prawdopodobieństwo błędu:**
     \[
     P = 10^{-0.2} \approx 0.631
     \]
     - Bardzo wysoka szansa błędu.

---

### **Dlaczego używa się kodowania ASCII w plikach FASTQ?**

- **Kompaktowość i efektywność:**
  - Reprezentowanie jakości jako pojedynczych znaków pozwala na **efektywne przechowywanie** dużych ilości danych.
- **Standardowość:**
  - ASCII jest **powszechnie akceptowanym standardem**, co ułatwia wymianę danych między różnymi systemami i oprogramowaniem.
- **Łatwość przetwarzania:**
  - Tekstowy format plików FASTQ z kodowaniem ASCII jest łatwy do **parsowania** i **analizy** za pomocą skryptów i narzędzi bioinformatycznych.

---

### **Znaczenie dla analiz bioinformatycznych**

- **Poprawna interpretacja jakości:**
  - **Znajomość offsetu** (np. phred+33) jest kluczowa dla poprawnego odczytu wartości jakości z pliku.
  - **Błędne założenie** offsetu może prowadzić do **niewłaściwej oceny jakości** odczytów i błędnych wyników analizy.
- **Filtracja odczytów:**
  - Wartości jakości są używane do **filtrowania niskiej jakości odczytów**, co wpływa na dokładność dalszych analiz, takich jak mapowanie czy wykrywanie wariantów.
- **Przetwarzanie danych:**
  - Narzędzia bioinformatyczne często **automatycznie rozpoznają kodowanie**, ale zawsze warto to **zweryfikować**, szczególnie przy pracy z danymi z różnych źródeł.

---

### **Podsumowanie**

- **Kodowanie ASCII** to metoda przypisywania numerycznych wartości znakom, używana w plikach FASTQ do reprezentowania jakości odczytów sekwencji.
- **Kluczowe elementy:**
  - **Offset (phred+33 lub phred+64)** determinuje, jak przekształcać znaki ASCII na wartości jakości.
  - **Wartości jakości** wpływają na **pewność** analizy sekwencji DNA.
- **Praktyczne wskazówki:**
  - Zawsze **sprawdzaj dokumentację** danych, aby znać używane kodowanie.
  - Przy analizie danych, **weryfikuj**, czy narzędzia poprawnie interpretują kodowanie jakości.

---

### **Dodatkowe informacje**

- **Tabela ASCII:**
  - Dostępna online, zawiera mapowanie znaków na wartości numeryczne (np. [Tabela ASCII](https://www.asciitable.com/)).
- **Standardy sekwencjonowania:**
  - Obecnie **Illumina** i inne platformy stosują **kodowanie phred+33**.
- **Konwersja kodowania:**
  - Jeśli masz pliki z różnymi kodowaniami, możesz użyć narzędzi takich jak **`seqtk`** czy **`FastQC`** do konwersji i weryfikacji.
