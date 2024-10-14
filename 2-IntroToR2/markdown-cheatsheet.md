
# **Formatowanie plików Markdown (`.md`) i RMarkdown (`.Rmd`)**

## **1. Struktura Pliku RMarkdown**

Plik `.Rmd` zawiera dwie główne części:
- **YAML Header**: Określa metadane dokumentu.
- **Treść Dokumentu**: Zawiera tekst, kod R i wyniki.

```markdown
---
title: "Tytuł Dokumentu"
author: "Twoje Imię"
date: "2024-04-27"
output: html_document
---

# Wprowadzenie
Twój tekst tutaj.
```

## **2. Nagłówki**

Używaj znaków `#` do tworzenia nagłówków. Liczba `#` oznacza poziom nagłówka.

```markdown
# Nagłówek poziomu 1
## Nagłówek poziomu 2
### Nagłówek poziomu 3
#### Nagłówek poziomu 4
##### Nagłówek poziomu 5
###### Nagłówek poziomu 6
```

## **3. Tekst z Formatowaniem**

- **Pogrubienie**: `**tekst**` lub `__tekst__`
- *Kursywa*: `*tekst*` lub `_tekst_`
- ***Pogrubiona kursywa***: `***tekst***`
- ~~Przekreślenie~~: `~~tekst~~`

```markdown
**Pogrubiony tekst**

*Tekst kursywą*

***Pogrubiony i kursywa***

~~Przekreślony tekst~~
```

## **4. Listy**

### **Lista Nienumerowana**

Używaj `-`, `*` lub `+` przed elementami listy.

```markdown
- Element 1
- Element 2
  - Podelement 2.1
  - Podelement 2.2
```

### **Lista Numerowana**

Używaj numerów zakończonych kropką.

```markdown
1. Pierwszy punkt
2. Drugi punkt
   1. Podpunkt 2.1
   2. Podpunkt 2.2
```

## **5. Linki i Obrazy**

### **Linki**

```markdown
[Tekst linku](https://przyklad.com)
```

### **Obrazy**

```markdown
![Alt tekst obrazu](ścieżka/do/obrazu.jpg)
```

## **6. Cytaty**

Użyj znaku `>` na początku linii.

```markdown
> To jest cytat.
```

## **7. Kod**

### **Kod Inline**

Użyj pojedynczych backticków `` ` ``.

```markdown
To jest `kod inline`.
```

### **Bloki Kodu**

Użyj potrójnych backticków ``` ``` ``` lub wcięć o 4 spacje.

```markdown
```r
# Przykładowy blok kodu R
summary(mtcars)
```
```

## **8. Tabele**

Użyj myślników `-` do oddzielenia nagłówków od treści i pionowych kresek `|` do rozdzielenia kolumn.

```markdown
| Kolumna 1 | Kolumna 2 | Kolumna 3 |
|-----------|-----------|-----------|
| Wartość 1 | Wartość 2 | Wartość 3 |
| Wartość 4 | Wartość 5 | Wartość 6 |
```

## **9. Linie Oddzielające**

Użyj trzech lub więcej gwiazdek `***`, myślników `---` lub podkreślników `___` na osobnej linii.

```markdown
---
```

## **10. Specjalne Elementy RMarkdown**

### **Bloki R**

Wstawianie kodu R z możliwością wyświetlenia wyników.

```markdown
```{r nazwa-bloku, opcje}
# Twój kod R tutaj
summary(mtcars)
```
```

### **Wstawianie Wyników R**

Wstaw wyniki R bez wyświetlania kodu.

```markdown
`r summary(mtcars)`
```

### **Parametry w YAML Header**

Możesz określić różne opcje wyjściowe, np. format PDF, HTML, Word.

```yaml
output:
  html_document:
    toc: true
    number_sections: true
```

## **11. Dodatkowe Skróty i Porady**

- **Lista punktowana z obrazkami lub linkami**:

  ```markdown
  - [x] Zadanie wykonane
  - [ ] Zadanie do zrobienia
  ```

- **Referencje do innych sekcji**:

  ```markdown
  [Przejdź do wprowadzenia](#wprowadzenie)
  ```

- **Komentarze w RMarkdown**:

  ```markdown
  <!-- To jest komentarz, który nie będzie widoczny w wygenerowanym dokumencie -->
  ```

## **12. Przydatne Narzędzia**

- **Podgląd w RStudio**: Kliknij przycisk `Knit`, aby wygenerować dokument.
- **RStudio Addins**: Użyj dodatków takich jak `formatR` do automatycznego formatowania kodu.
- **RMarkdown Cheat Sheets**: Dostępne na stronie [RStudio Cheat Sheets](https://rstudio.com/resources/cheatsheets/).

---

## **Przykładowy Dokument RMarkdown**

```markdown
---
title: "Przykładowy Raport"
author: "Jan Kowalski"
date: "2024-04-27"
output: html_document
---

# Wprowadzenie

To jest przykładowy raport stworzony w RMarkdown. Poniżej znajduje się przykładowy kod R.

## Opis Środowiska

Zainstalowane pakiety:
- `ggplot2`
- `dplyr`

## Przykładowy Kod

```{r}
# Wyświetlenie powitania
print("Witaj świecie!")
```

## Wnioski

Nauczyłem się podstaw formatowania dokumentów w RMarkdown.
```

---

## **Źródła i Dodatkowe Materiały**

- [Dokumentacja Markdown](https://www.markdownguide.org/basic-syntax/)
- [RMarkdown - RStudio](https://rmarkdown.rstudio.com/)
- [Cheat Sheet - RStudio](https://rstudio.com/resources/cheatsheets/)

---

**Powodzenia w tworzeniu czytelnych i profesjonalnych dokumentów z użyciem Markdown i RMarkdown!**
```

---

### **Jak Zapisz Plik `.md`**

1. **Kopiowanie Treści:**
   - Skopiuj cały tekst powyżej (wszystko między liniami `---`).

2. **Tworzenie Pliku:**
   - Otwórz swój edytor tekstu (np. Notepad, VS Code, RStudio).
   - Wklej skopiowaną treść.

3. **Zapisz Plik:**
   - Wybierz opcję **Zapisz jako**.
   - Nazwij plik, na przykład `sciaga_formatowanie.md`.
   - Upewnij się, że rozszerzenie pliku to `.md`.

4. **Przeglądanie:**
   - Możesz otworzyć plik `.md` w przeglądarce lub edytorze wspierającym podgląd Markdown, aby zobaczyć sformatowany dokument.

---

**Dodatkowe Wskazówki:**

- **RMarkdown vs Markdown:**
  - Pliki `.Rmd` zawierają zarówno Markdown, jak i kod R, co umożliwia dynamiczne generowanie treści.
  - Pliki `.md` zawierają tylko Markdown i są używane głównie do dokumentacji, README, itp.

- **Narzędzia Online do Podglądu Markdown:**
  - [Dillinger](https://dillinger.io/)
  - [Markdown Live Preview](https://markdownlivepreview.com/)

- **RStudio:**
  - RStudio oferuje wbudowany podgląd dla plików `.Rmd` i `.md`, co ułatwia pracę z dokumentami.

---
