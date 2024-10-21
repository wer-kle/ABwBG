# Podstawowe funkcje w R z pakietu base

## 1. Podstawowe operacje matematyczne

### Arytmetyka

```r
a <- 5 + 3    # Dodawanie: a = 8
b <- 10 - 2   # Odejmowanie: b = 8
c <- 4 * 2    # Mnożenie: c = 8
d <- 16 / 2   # Dzielenie: d = 8
```

### Potęgowanie i modulo

```r
e <- 2 ^ 3    # Potęgowanie: e = 8 (2 do potęgi 3)
f <- 17 %% 3  # Modulo: f = 2 (reszta z dzielenia 17 przez 3)
g <- 17 %/% 3 # Dzielenie całkowite: g = 5 (17 podzielone przez 3 bez reszty)
```

### Funkcje matematyczne

```r
h <- abs(-5)        # Wartość bezwzględna: h = 5
i <- sqrt(16)       # Pierwiastek kwadratowy: i = 4
j <- exp(1)         # Funkcja wykładnicza: j ≈ 2.71828
k <- log(10)        # Logarytm naturalny: k ≈ 2.302585
l <- log10(1000)    # Logarytm dziesiętny: l = 3
m <- round(3.14159, digits = 2) # Zaokrąglanie: m = 3.14
n <- ceiling(3.14)  # Zaokrąglenie w górę: n = 4
o <- floor(3.14)    # Zaokrąglenie w dół: o = 3
```

### Funkcje trygonometryczne

```r
p <- sin(pi/2)      # Sinus: p = 1
q <- cos(0)         # Cosinus: q = 1
r <- tan(pi/4)      # Tangens: r = 1
```

---

## 2. Funkcje statystyczne
*Analiza danych genomowych często wymaga zastosowania zaawansowanych metod statystycznych do interpretacji wyników eksperymentów wysokoprzepustowych, takich jak sekwencjonowanie DNA czy RNA. Znajomość podstawowych parametrów statystycznych, takich jak średnia, mediana, wariancja czy odchylenie standardowe, jest kluczowa do oceny jakości danych i identyfikacji istotnych biologicznie zmian.*

### Tworzenie wektora danych

```r
dane <- c(4, 8, 15, 16, 23, 42)
```

### Podstawowe statystyki

```r
srednia <- mean(dane)     # Średnia arytmetyczna
mediana <- median(dane)   # Mediana
wariancja <- var(dane)    # Wariancja
odchylenie <- sd(dane)    # Odchylenie standardowe
minimum <- min(dane)      # Wartość minimalna
maksimum <- max(dane)     # Wartość maksymalna
zakres <- range(dane)     # Zakres (min i max)
suma <- sum(dane)         # Suma elementów
iloczyn <- prod(dane)     # Iloczyn elementów
```

### Podsumowanie danych

```r
summary(dane)             # Podsumowanie statystyczne
```

### Tabela częstości

```r
kategorie <- c("A", "B", "A", "C", "B", "B", "A")
tabela <- table(kategorie) # Liczność poszczególnych kategorii
```

### Korelacja i kowariancja

```r
x <- c(1, 2, 3, 4, 5)
y <- c(2, 4, 6, 8, 10)
korelacja <- cor(x, y)    # Korelacja między x i y
kowariancja <- cov(x, y)  # Kowariancja między x i y
```

---

## 3. Operacje na wektorach i indeksowanie
*Dane genomowe są zwykle przechowywane w dużych wektorach, macierzach lub listach. Umiejętność efektywnego indeksowania, filtrowania i manipulacji tymi strukturami danych jest niezbędna do przetwarzania i analizowania dużych zestawów danych genomowych, takich jak ekspresja genów czy warianty genetyczne.*

### Tworzenie wektorów

```r
wektor1 <- c(1, 2, 3, 4, 5)
wektor2 <- seq(from = 1, to = 10, by = 2) # Sekwencja: 1, 3, 5, 7, 9
wektor3 <- rep(x = 0, times = 5)          # Powtarzanie: 0, 0, 0, 0, 0
```

### Indeksowanie

```r
pierwszy <- wektor1[1]       # Pierwszy element: 1
ostatni <- wektor1[length(wektor1)] # Ostatni element: 5
podwektor <- wektor1[2:4]    # Elementy od 2 do 4: 2, 3, 4
```

### Filtrowanie

```r
powyzej_2 <- wektor1[wektor1 > 2] # Elementy większe niż 2: 3, 4, 5
```

### Sortowanie i porządkowanie

```r
malejaco <- sort(wektor1, decreasing = TRUE) # Sortowanie malejąco
indeksy_sort <- order(wektor1) # Indeksy sortujące
rangi <- rank(wektor1)         # Rangi elementów
```

---

## 4. Praca z ramkami danych (data frames)
*Ramki danych są podstawową strukturą do przechowywania danych w R, w tym danych genomowych. BioConductor wykorzystuje je do organizacji danych o genach, próbkach, anotacjach i wynikach analiz. Umiejętność manipulacji ramkami danych, dodawania i usuwania kolumn czy wierszy oraz dostępu do poszczególnych elementów jest kluczowa.*

### Tworzenie ramki danych

```r
dane_ramka <- data.frame(
  imie = c("Anna", "Piotr", "Maria"),
  wiek = c(28, 34, 23),
  miasto = c("Warszawa", "Kraków", "Gdańsk")
)
```

### Dostęp do danych

```r
dane_ramka$imie          # Dostęp do kolumny 'imie'
dane_ramka[1, ]          # Pierwszy wiersz
dane_ramka[, "wiek"]     # Kolumna 'wiek'
```

### Informacje o danych

```r
str(dane_ramka)          # Struktura obiektu
head(dane_ramka)         # Pierwsze wiersze
tail(dane_ramka)         # Ostatnie wiersze
dim(dane_ramka)          # Wymiary (wiersze, kolumny)
nrow(dane_ramka)         # Liczba wierszy
ncol(dane_ramka)         # Liczba kolumn
names(dane_ramka)        # Nazwy kolumn
rownames(dane_ramka)     # Nazwy wierszy
```

### Dodawanie i usuwanie kolumn

```r
dane_ramka$plec <- c("K", "M", "K") # Dodawanie kolumny 'plec'
dane_ramka$miasto <- NULL            # Usuwanie kolumny 'miasto'
```

---

## 5. Funkcje logiczne i warunkowe
*Operacje logiczne i warunkowe są niezbędne do filtrowania danych genomowych np. na podstawie fenotypów.*

### Operatory logiczne

```r
x <- 5
y <- 10

x == y     # Czy x jest równe y? FALSE
x != y     # Czy x nie jest równe y? TRUE
x < y      # Czy x jest mniejsze od y? TRUE
x > y      # Czy x jest większe od y? FALSE
x <= y     # Czy x jest mniejsze lub równe y? TRUE
x >= y     # Czy x jest większe lub równe y? FALSE
```

### Operatory logiczne AND, OR, NOT

```r
(x > 0) & (y > 0)  # Czy x i y są większe od 0? TRUE
(x > 0) | (y < 0)  # Czy x jest większe od 0 lub y mniejsze od 0? TRUE
!(x == y)          # Negacja: Czy x nie jest równe y? TRUE
```

### Instrukcje warunkowe

```r
if (x > y) {
  wynik <- "x jest większe od y"
} else if (x < y) {
  wynik <- "x jest mniejsze od y"
} else {
  wynik <- "x jest równe y"
}
```

### Funkcja ifelse

```r
wektor <- c(-1, 0, 1)
wyniki <- ifelse(wektor > 0, "Dodatni", "Niedodatni")
# Wynik: "Niedodatni", "Niedodatni", "Dodatni"
```

### Funkcje any i all

```r
any(wektor > 0)  # Czy jakikolwiek element jest większy od 0? TRUE
all(wektor > 0)  # Czy wszystkie elementy są większe od 0? FALSE
```

### Funkcja which

```r
indeksy_pozytywne <- which(wektor > 0) # Indeksy elementów większych od 0
```

---

## 6. Pętle i iteracje
*Chociaż w R preferowane są operacje wektorowe, w analizie danych genomowych często konieczne jest zastosowanie pętli do iteracyjnego przetwarzania danych, szczególnie przy niestandardowych analizach lub gdy funkcje wektorowe nie są dostępne. Znajomość pracy na pętlach jest też bardzo przydatna do pracy w Pythonie.*

### Pętla for

```r
for (i in 1:5) {
  print(paste("Iteracja", i))
}
```

### Pętla while

```r
i <- 1
while (i <= 5) {
  print(paste("i =", i))
  i <- i + 1
}
```

### Pętla repeat

```r
i <- 1
repeat {
  print(i)
  if (i >= 5) break
  i <- i + 1
}
```

### Funkcja apply

```r
macierz <- matrix(1:9, nrow = 3)
suma_wierszy <- apply(macierz, 1, sum) # Suma wierszy
suma_kolumn <- apply(macierz, 2, sum)  # Suma kolumn
```

### Funkcje lapply i sapply

```r
lista <- list(a = 1:5, b = 6:10)
suma_lista <- lapply(lista, sum)   # Zwraca listę sum
suma_wektor <- sapply(lista, sum)  # Zwraca wektor sum
```

### Funkcja tapply

```r
grupy <- c("G1", "G2", "G1", "G2", "G1")
wartosci <- c(10, 20, 30, 40, 50)
srednia_grup <- tapply(wartosci, grupy, mean) # Średnie w grupach
```

---

## 7. Funkcje wejścia/wyjścia
*Importowanie i eksportowanie danych jest fundamentalne w analizie genomowej. Dane mogą pochodzić z różnych źródeł i mogą być zapisywane w różnych formatach (np. FASTQ, BAM, VCF). Znajomość funkcji do czytania i zapisywania danych pozwala na efektywne zarządzanie plikami i integrację różnych typów danych omicznych.*

### Czytanie danych z pliku

```r
# dane <- read.csv("dane.csv", header = TRUE, sep = ",") #jeżeli nie mamy zdefiniowanego katalogu roboczego, musimy wpisać pełny adres ścieżki do pliku wejściowego
```

### Zapisywanie danych do pliku

```r
# write.csv(dane_ramka, "wyniki.csv", row.names = FALSE)
```

### Czytanie linii z pliku

```r
# linie <- readLines("tekst.txt")
```

### Zapisywanie linii do pliku

```r
# writeLines(c("Linia 1", "Linia 2"), "nowy_plik.txt")
```

---

## 8. Funkcje związane z łańcuchami znaków
*Dane genomowe zawierają dużo informacji tekstowych, takich jak nazwy genów, identyfikatory sekwencji czy anotacje funkcjonalne, a często same sekwencje są zapisane jako łańcuchy znaków. Manipulacja łańcuchami znaków jest ważna przy przetwarzaniu informacji, takich jak ekstrakcja określonych fragmentów sekwencji czy wyszukiwanie wzorców*

### Łączenie tekstu

```r
tekst <- paste("Witaj", "Świecie", sep = " ")   # "Witaj Świecie"
```

### Podłańcuchy

```r
napis <- "Bioinformatyka"
podnapis <- substr(napis, start = 4, stop = 7)  # "info"
```

### Rozdzielanie tekstu

```r
rozdzielony <- strsplit("A,B,C,D", split = ",")
# Lista z elementami "A", "B", "C", "D"
```

### Liczba znaków

```r
dlugosc <- nchar(napis)  # 14
```

### Zmiana wielkości liter

```r
duze <- toupper(napis)   # "BIOINFORMATYKA"
male <- tolower(napis)   # "bioinformatyka"
```

### Wyszukiwanie wzorca

```r
wektor_txt <- c("kot", "pies", "kotka")
dopasowanie <- grep("kot", wektor_txt)   # Indeksy elementów zawierających "kot"
czy_zawiera <- grepl("kot", wektor_txt)  # TRUE, FALSE dla każdego elementu
```

### Zastępowanie tekstu

```r
tekst_zmieniony <- sub("kot", "pies", "kot jest tu")   # "pies jest tu"
tekst_zmieniony_all <- gsub("kot", "pies", "kot kotek kotka") # "pies piesek pieska"
```

---

## 9. Funkcje statystyczne i probabilistyczne
*Zaawansowane analizy statystyczne są kluczowe w genomice. Obejmują one modelowanie statystyczne, testy istotności, analizy rozkładów czy estymację parametrów. Znajomość tych funkcji pozwala na przeprowadzenie złożonych analiz, takich jak chociażby GWAS.*

### Rozkład normalny

```r
dnorm(0, mean = 0, sd = 1)      # Gęstość w punkcie 0
pnorm(1.96, mean = 0, sd = 1)   # Dystrybuanta w punkcie 1.96
qnorm(0.975, mean = 0, sd = 1)  # Kwantyl 97.5%
losowe <- rnorm(100, mean = 0, sd = 1)  # 100 losowych liczb
```

### Test t-Studenta

```r
grupa1 <- c(5, 6, 7, 8, 9)
grupa2 <- c(7, 8, 9, 10, 11)
test_t <- t.test(grupa1, grupa2)
```

---

## 10. Funkcje systemowe i diagnostyczne

### Lista obiektów w środowisku

```r
ls()
```

### Usuwanie obiektu

```r
rm(a)  # Usuwa obiekt 'a'
```

### Bieżący katalog roboczy

```r
getwd()
```

### Ustawienie katalogu roboczego

```r
setwd("C:/Nowy/Katalog")
```

### Pomoc dla funkcji

```r
help(mean) # lub ?mean
```

### Wyszukiwanie w dokumentacji

```r
help.search("regression")
```

### Przykłady użycia funkcji

```r
example(mean)
```

### Pomiar czasu wykonania kodu

```r
czas <- system.time({
  suma_duza <- sum(1:1000000)
})
```

### Śledzenie błędów

```r
# Jeśli wystąpi błąd, użyj traceback()
```

---

## 11. Tworzenie własnych funkcji
*Pisanie własnych funkcji umożliwia automatyzację powtarzalnych zadań i tworzenie niestandardowych narzędzi analitycznych dostosowanych do specyficznych potrzeb badawczych. W analizie danych genomowych może to obejmować tworzenie funkcji do specyficznego przetwarzania sekwencji czy agregacji wyników.*

### Definiowanie funkcji

```r
dodaj <- function(x, y) {
  wynik <- x + y
  return(wynik)
}
```

### Wywołanie funkcji

```r
suma <- dodaj(3, 5)  # suma = 8
```

---

## 12. Przekształcenia typów danych
*Dane genomowe mogą występować w różnych formatach i typach. Umiejętność konwersji między typami danych, takimi jak liczby, czynniki (faktory) czy łańcuchy znaków, jest niezbędna do integracji danych z różnych źródeł i przygotowania ich do analizy.*

### Konwersja na liczby

```r
tekst <- "123"
liczba <- as.numeric(tekst)  # liczba = 123
```

### Konwersja na tekst

```r
liczba <- 123
tekst <- as.character(liczba)  # tekst = "123"
```

### Konwersja na czynnik

```r
wektor <- c("tak", "nie", "tak")
czynnik <- as.factor(wektor)
```

---

## 13. Funkcje dotyczące dat i czasu

### Bieżąca data i czas

```r
data <- Sys.Date()   # Np. "2023-10-10"
czas <- Sys.time()   # Np. "2023-10-10 14:30:00"
```

### Formatowanie daty

```r
format(data, "%d-%m-%Y")  # Np. "10-10-2023"
```

### Różnica między datami

```r
data1 <- as.Date("2023-01-01")
data2 <- as.Date("2023-12-31")
roznica <- difftime(data2, data1, units = "days")  # Liczba dni między datami
```

---

## 14. Losowanie i permutacje
*Permutacje i losowe próbkowanie są często wykorzystywane w analizach genomowych do oceny istotności statystycznej, np. w testach permutacyjnych czy podczas tworzenia rozkładów null w analizach asocjacyjnych.*

### Losowy wybór elementów

```r
elementy <- c("A", "B", "C", "D", "E")
proba <- sample(elementy, size = 3, replace = FALSE)  # Losowe 3 różne elementy
```

### Ustawienie ziarna dla powtarzalności

```r
set.seed(123)
proba_powtarzalna <- sample(elementy, size = 3)
```

---

## 15. Inne przydatne funkcje

### Generowanie sekwencji liczbowych

```r
sekwencja <- 1:10          # Liczby od 1 do 10
sekwencja2 <- seq(0, 1, by = 0.1) # Od 0 do 1 co 0.1
```

### Losowe permutacje

```r
permutacja <- sample(1:5)  # Losowa permutacja liczb 1 do 5
```

### Zastosowanie funkcji do elementów listy

```r
lista <- list(a = 1:3, b = 4:6)
suma_elementow <- lapply(lista, sum)  # Sumuje elementy każdej podlisty
```

---

## 16. Profilowanie i debugowanie

### Debugowanie funkcji

```r
debug(nazwa_funkcji)
```

### Wyłączanie debugowania

```r
undebug(nazwa_funkcji)
```

### Śledzenie błędów

```r
# Opcja 'error' pozwala na wejście do debuggera po błędzie
options(error = recover)
```

---

**Zachęcam do eksperymentowania z każdą z funkcji!**
