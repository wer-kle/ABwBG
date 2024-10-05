# -------------------------------
# Podstawowe funkcje w R z pakietu base
# -------------------------------

# 1. Podstawowe operacje matematyczne

# Arytmetyka
a <- 5 + 3    # Dodawanie: a = 8
b <- 10 - 2   # Odejmowanie: b = 8
c <- 4 * 2    # Mnożenie: c = 8
d <- 16 / 2   # Dzielenie: d = 8

# Potęgowanie i modulo
e <- 2 ^ 3    # Potęgowanie: e = 8 (2 do potęgi 3)
f <- 17 %% 3  # Modulo: f = 2 (reszta z dzielenia 17 przez 3)
g <- 17 %/% 3 # Dzielenie całkowite: g = 5 (17 podzielone przez 3 bez reszty)

# Funkcje matematyczne
h <- abs(-5)        # Wartość bezwzględna: h = 5
i <- sqrt(16)       # Pierwiastek kwadratowy: i = 4
j <- exp(1)         # Funkcja wykładnicza: j ≈ 2.71828
k <- log(10)        # Logarytm naturalny: k ≈ 2.302585
l <- log10(1000)    # Logarytm dziesiętny: l = 3
m <- round(3.14159, digits = 2) # Zaokrąglanie: m = 3.14
n <- ceiling(3.14)  # Zaokrąglenie w górę: n = 4
o <- floor(3.14)    # Zaokrąglenie w dół: o = 3

# Funkcje trygonometryczne
p <- sin(pi/2)      # Sinus: p = 1
q <- cos(0)         # Cosinus: q = 1
r <- tan(pi/4)      # Tangens: r = 1

# -------------------------------

# 2. Funkcje statystyczne

# Tworzenie wektora danych
dane <- c(4, 8, 15, 16, 23, 42)

# Podstawowe statystyki
srednia <- mean(dane)     # Średnia arytmetyczna
mediana <- median(dane)   # Mediana
wariancja <- var(dane)    # Wariancja
odchylenie <- sd(dane)    # Odchylenie standardowe
minimum <- min(dane)      # Wartość minimalna
maksimum <- max(dane)     # Wartość maksymalna
zakres <- range(dane)     # Zakres (min i max)
suma <- sum(dane)         # Suma elementów
iloczyn <- prod(dane)     # Iloczyn elementów

# Podsumowanie danych
summary(dane)             # Podsumowanie statystyczne

# Tabela częstości
kategorie <- c("A", "B", "A", "C", "B", "B", "A")
tabela <- table(kategorie) # Liczność poszczególnych kategorii

# Korelacja i kowariancja
x <- c(1, 2, 3, 4, 5)
y <- c(2, 4, 6, 8, 10)
korelacja <- cor(x, y)    # Korelacja między x i y
kowariancja <- cov(x, y)  # Kowariancja między x i y

# -------------------------------

# 3. Operacje na wektorach i indeksowanie

# Tworzenie wektorów
wektor1 <- c(1, 2, 3, 4, 5)
wektor2 <- seq(from = 1, to = 10, by = 2) # Sekwencja: 1, 3, 5, 7, 9
wektor3 <- rep(x = 0, times = 5)          # Powtarzanie: 0, 0, 0, 0, 0

# Indeksowanie
pierwszy <- wektor1[1]       # Pierwszy element: 1
ostatni <- wektor1[length(wektor1)] # Ostatni element: 5
podwektor <- wektor1[2:4]    # Elementy od 2 do 4: 2, 3, 4

# Filtrowanie
filtrowany <- wektor1[wektor1 > 2] # Elementy większe niż 2: 3, 4, 5

# Sortowanie i porządkowanie
posortowany <- sort(wektor1, decreasing = TRUE) # Sortowanie malejąco
indeksy_sort <- order(wektor1) # Indeksy sortujące
rangi <- rank(wektor1)         # Rangi elementów

# -------------------------------

# 4. Praca z ramkami danych (data frames)

# Tworzenie ramki danych
dane_ramka <- data.frame(
  imie = c("Anna", "Piotr", "Maria"),
  wiek = c(28, 34, 23),
  miasto = c("Warszawa", "Kraków", "Gdańsk")
)

# Dostęp do danych
dane_ramka$imie          # Dostęp do kolumny 'imie'
dane_ramka[1, ]          # Pierwszy wiersz
dane_ramka[ , "wiek"]    # Kolumna 'wiek'

# Informacje o danych
str(dane_ramka)          # Struktura obiektu
head(dane_ramka)         # Pierwsze wiersze
tail(dane_ramka)         # Ostatnie wiersze
dim(dane_ramka)          # Wymiary (wiersze, kolumny)
nrow(dane_ramka)         # Liczba wierszy
ncol(dane_ramka)         # Liczba kolumn
names(dane_ramka)        # Nazwy kolumn
rownames(dane_ramka)     # Nazwy wierszy

# Dodawanie i usuwanie kolumn
dane_ramka$plec <- c("K", "M", "K") # Dodawanie kolumny 'plec'
dane_ramka$miasto <- NULL            # Usuwanie kolumny 'miasto'

# -------------------------------

# 5. Funkcje logiczne i warunkowe

# Operatory logiczne
x <- 5
y <- 10

x == y     # Czy x jest równe y? FALSE
x != y     # Czy x nie jest równe y? TRUE
x < y      # Czy x jest mniejsze od y? TRUE
x > y      # Czy x jest większe od y? FALSE
x <= y     # Czy x jest mniejsze lub równe y? TRUE
x >= y     # Czy x jest większe lub równe y? FALSE

# Operatory logiczne AND, OR, NOT
(x > 0) & (y > 0)  # Czy x i y są większe od 0? TRUE
(x > 0) | (y < 0)  # Czy x jest większe od 0 lub y mniejsze od 0? TRUE
!(x == y)          # Negacja: Czy x nie jest równe y? TRUE

# Instrukcje warunkowe
if (x > y) {
  wynik <- "x jest większe od y"
} else if (x < y) {
  wynik <- "x jest mniejsze od y"
} else {
  wynik <- "x jest równe y"
}

# Funkcja ifelse
wektor <- c(-1, 0, 1)
wyniki <- ifelse(wektor > 0, "Dodatni", "Niedodatni") # "Niedodatni", "Niedodatni", "Dodatni"

# Funkcje any i all
any(wektor > 0)  # Czy jakikolwiek element jest większy od 0? TRUE
all(wektor > 0)  # Czy wszystkie elementy są większe od 0? FALSE

# Funkcja which
indeksy_pozytywne <- which(wektor > 0) # Indeksy elementów większych od 0

# -------------------------------

# 6. Pętle i iteracje

# Pętla for
for (i in 1:5) {
  print(paste("Iteracja", i))
}

# Pętla while
i <- 1
while (i <= 5) {
  print(paste("i =", i))
  i <- i + 1
}

# Pętla repeat
i <- 1
repeat {
  print(i)
  if (i >= 5) break
  i <- i + 1
}

# Funkcja apply
macierz <- matrix(1:9, nrow = 3)
suma_wierszy <- apply(macierz, 1, sum) # Suma wierszy
suma_kolumn <- apply(macierz, 2, sum)  # Suma kolumn

# Funkcje lapply i sapply
lista <- list(a = 1:5, b = 6:10)
suma_lista <- lapply(lista, sum)   # Zwraca listę sum
suma_wektor <- sapply(lista, sum)  # Zwraca wektor sum

# Funkcja tapply
grupy <- c("G1", "G2", "G1", "G2", "G1")
wartosci <- c(10, 20, 30, 40, 50)
srednia_grup <- tapply(wartosci, grupy, mean) # Średnie w grupach

# -------------------------------

# 7. Funkcje wejścia/wyjścia

# Czytanie danych z pliku (przykład z plikiem 'dane.csv')
# dane <- read.csv("dane.csv", header = TRUE, sep = ",")

# Zapisywanie danych do pliku
# write.csv(dane_ramka, "wyniki.csv", row.names = FALSE)

# Czytanie linii z pliku (przykład z plikiem 'tekst.txt')
# linie <- readLines("tekst.txt")

# Zapisywanie linii do pliku
# writeLines(c("Linia 1", "Linia 2"), "nowy_plik.txt")

# -------------------------------

# 8. Funkcje związane z łańcuchami znaków

# Łączenie tekstu
tekst <- paste("Witaj", "Świecie", sep = " ")   # "Witaj Świecie"

# Podłańcuchy
napis <- "Bioinformatyka"
podnapis <- substr(napis, start = 4, stop = 7)  # "info"

# Rozdzielanie tekstu
rozdzielony <- strsplit("A,B,C,D", split = ",") # Lista z elementami "A", "B", "C", "D"

# Liczba znaków
dlugosc <- nchar(napis)  # 14

# Zmiana wielkości liter
duze <- toupper(napis)   # "BIOINFORMATYKA"
male <- tolower(napis)   # "bioinformatyka"

# Wyszukiwanie wzorca
wektor_txt <- c("kot", "pies", "kotka")
dopasowanie <- grep("kot", wektor_txt)   # Indeksy elementów zawierających "kot"
czy_zawiera <- grepl("kot", wektor_txt)  # TRUE, FALSE dla każdego elementu

# Zastępowanie tekstu
tekst_zmieniony <- sub("kot", "pies", "kot jest tu")   # "pies jest tu"
tekst_zmieniony_all <- gsub("kot", "pies", "kot kotek kotka") # "pies piesek pieska"

# -------------------------------

# 9. Funkcje statystyczne i probabilistyczne

# Rozkład normalny
dnorm(0, mean = 0, sd = 1)      # Gęstość w punkcie 0
pnorm(1.96, mean = 0, sd = 1)   # Dystrybuanta w punkcie 1.96
qnorm(0.975, mean = 0, sd = 1)  # Kwantyl 97.5%
losowe <- rnorm(100, mean = 0, sd = 1)  # 100 losowych liczb

# Test t-Studenta (przykładowo)
grupa1 <- c(5, 6, 7, 8, 9)
grupa2 <- c(7, 8, 9, 10, 11)
test_t <- t.test(grupa1, grupa2)

# -------------------------------

# 10. Funkcje systemowe i diagnostyczne

# Lista obiektów w środowisku
ls()

# Usuwanie obiektu
rm(a)  # Usuwa obiekt 'a'

# Bieżący katalog roboczy
getwd()

# Ustawienie katalogu roboczego
setwd("C:/Nowy/Katalog")

# Pomoc dla funkcji
help(mean) lub ?mean

# Wyszukiwanie w dokumentacji
help.search("regression")

# Przykłady użycia funkcji
example(mean)

# Pomiar czasu wykonania kodu
czas <- system.time({
  suma_duza <- sum(1:1000000)
})

# Śledzenie błędów
# Jeśli wystąpi błąd, użyj traceback()

# -------------------------------

# 11. Tworzenie własnych funkcji

# Definiowanie funkcji
dodaj <- function(x, y) {
  wynik <- x + y
  return(wynik)
}

# Wywołanie funkcji
suma <- dodaj(3, 5)  # suma = 8

# -------------------------------

# 12. Przekształcenia typów danych

# Konwersja na liczby
tekst <- "123"
liczba <- as.numeric(tekst)  # liczba = 123

# Konwersja na tekst
liczba <- 123
tekst <- as.character(liczba)  # tekst = "123"

# Konwersja na czynnik
wektor <- c("tak", "nie", "tak")
czynnik <- as.factor(wektor)

# -------------------------------

# 13. Funkcje dotyczące dat i czasu

# Bieżąca data i czas
data <- Sys.Date()   # Np. "2023-10-10"
czas <- Sys.time()   # Np. "2023-10-10 14:30:00"

# Formatowanie daty
format(data, "%d-%m-%Y")  # Np. "10-10-2023"

# Różnica między datami
data1 <- as.Date("2023-01-01")
data2 <- as.Date("2023-12-31")
roznica <- difftime(data2, data1, units = "days")  # Liczba dni między datami

# -------------------------------

# 14. Losowanie i permutacje

# Losowy wybór elementów
elementy <- c("A", "B", "C", "D", "E")
proba <- sample(elementy, size = 3, replace = FALSE)  # Losowe 3 różne elementy

# Ustawienie ziarna dla powtarzalności
set.seed(123)
proba_powtarzalna <- sample(elementy, size = 3)

# -------------------------------

# 15. Inne przydatne funkcje

# Generowanie sekwencji liczbowych
sekwencja <- 1:10          # Liczby od 1 do 10
sekwencja2 <- seq(0, 1, by = 0.1) # Od 0 do 1 co 0.1

# Losowe permutacje
permutacja <- sample(1:5)  # Losowa permutacja liczb 1 do 5

# Zastosowanie funkcji do elementów listy
lista <- list(a = 1:3, b = 4:6)
suma_elementow <- lapply(lista, sum)  # Sumuje elementy każdej podlisty

# -------------------------------

# 16. Profilowanie i debugowanie

# Debugowanie funkcji
debug(nazwa_funkcji)
# Wyłączanie debugowania
undebug(nazwa_funkcji)

# Śledzenie błędów
# Opcja 'error' pozwala na wejście do debuggera po błędzie
options(error = recover)

# -------------------------------

# Zachęcam do eksperymentowania i modyfikowania kodu!
