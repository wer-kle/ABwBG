## Zadanie 1: Podstawowe operacje matematyczne

1. Utwórz dwie zmienne `x` i `y` o wartościach 12 i 5. Wykonaj na nich operacje dodawania, odejmowania, mnożenia oraz dzielenia.
2. Oblicz resztę z dzielenia `x` przez `y` (operator modulo) oraz dzielenie całkowite `x` przez `y`.
3. Znajdź wartość bezwzględną z liczby `-15`, pierwiastek kwadratowy z `25`, oraz logarytm dziesiętny z `10000`.
4. Zaokrąglij liczbę `3.45678` do dwóch miejsc po przecinku, a następnie zastosuj funkcje `ceiling()` i `floor()`.

## Zadanie 2: Funkcje statystyczne

1. Utwórz wektor `dane` zawierający liczby `3, 7, 12, 14, 21, 33`.
2. Oblicz średnią, medianę, wariancję, odchylenie standardowe, minimum, maksimum oraz sumę elementów wektora.
3. Sprawdź, jakie statystyki podsumowujące można uzyskać przy użyciu funkcji `summary()`.
4. Utwórz wektor `kolory` zawierający elementy `"czerwony", "zielony", "czerwony", "niebieski", "zielony"`. Policz, ile razy występuje każdy kolor.

## Zadanie 3: Operacje na wektorach i indeksowanie

1. Stwórz wektor `liczby` zawierający liczby od `1` do `10`.
2. Z tego wektora wybierz pierwszy element, ostatni element oraz podwektor zawierający elementy od 3 do 7.
3. Wybierz z wektora `liczby` te elementy, które są większe niż `5`.
4. Posortuj wektor `liczby` malejąco, a następnie znajdź indeksy posortowanych elementów oraz ich rangi.

## Zadanie 4: Praca z ramkami danych (data frames)

1. Stwórz ramkę danych `osoby` z kolumnami `imie`, `wiek`, `miasto`, zawierającą dane dla 3 osób.
2. Wyświetl tylko kolumnę `imie`, pierwszy wiersz oraz kolumnę `miasto`.
3. Dodaj do ramki `osoby` kolumnę `plec` i przypisz do niej odpowiednie wartości.
4. Usuń kolumnę `miasto` z ramki danych.

## Zadanie 5: Funkcje logiczne i warunkowe

1. Utwórz dwie zmienne `a` i `b` o wartościach `8` i `15`. Sprawdź relacje między nimi, wykorzystując operatory logiczne (`==`, `!=`, `<`, `>`, `<=`, `>=`).
2. Sprawdź, czy zmienna `a` jest większa od `0` oraz czy `b` jest mniejsza niż `20`, używając operatorów logicznych `&` oraz `|`.
3. Stwórz zmienną `ocena` o wartości `75`. Użyj instrukcji `ifelse()`, aby przypisać jej etykietę `"Pozytywna"` lub `"Negatywna"` w zależności od tego, czy jest większa lub równa `50`.
4. Sprawdź, czy jakikolwiek element wektora `c(2, -5, 7, -3)` jest mniejszy od zera oraz czy wszystkie są większe od zera.

## Zadanie 6: Pętle i iteracje

1. Użyj pętli `for`, aby wypisać kolejne liczby od 1 do 5.
2. Utwórz pętlę `while`, która będzie wypisywać wartości `i`, dopóki `i` jest mniejsze lub równe `3`.
3. Użyj pętli `repeat`, aby wypisać wartości `j`, zaczynając od `1`, aż do momentu osiągnięcia `j = 4`.
4. Utwórz macierz `macierz_3x3` o wymiarach `3x3`, zawierającą liczby od `1` do `9`. Oblicz sumę wierszy i sumę kolumn tej macierzy za pomocą funkcji `apply()`.

## Zadanie 7: Funkcje wejścia/wyjścia

1. Utwórz ramkę danych `produkty` z kolumnami `produkt`, `cena` i `ilosc`. Zapisz tę ramkę do pliku `produkty.csv`.
2. Wczytaj dane z pliku `produkty.csv` do nowej ramki danych.
3. Utwórz plik tekstowy `notatka.txt` zawierający kilka linii tekstu. Wczytaj wszystkie linie z tego pliku za pomocą funkcji `readLines()`.

## Zadanie 8: Funkcje związane z łańcuchami znaków

1. Stwórz łańcuch `"R jest super!"`. Sprawdź jego długość za pomocą funkcji `nchar()`.
2. Wyciągnij podłańcuch ze swojego łańcucha, zawierający tylko słowo `"super"`.
3. Rozdziel łańcuch `"A-B-C-D"` na elementy oddzielone myślnikiem `-`.
4. Zastąp w łańcuchu `"Kocham koty"` wszystkie wystąpienia słowa `"koty"` na `"psy"`.

## Zadanie 9: Funkcje statystyczne i probabilistyczne

1. Wygeneruj 100 losowych liczb o rozkładzie normalnym ze średnią `0` i odchyleniem standardowym `1`. Zapisz je do wektora `losowe`.
2. Narysuj histogram wektora `losowe` za pomocą funkcji `hist()`.
3. Wykonaj test t-Studenta, porównując dwie grupy danych: `c(3, 5, 7, 9)` i `c(2, 4, 6, 8)`.

## Zadanie 10: Tworzenie własnych funkcji

1. Zdefiniuj funkcję `pomnoz()`, która przyjmuje dwa argumenty i zwraca ich iloczyn.
2. Zastosuj funkcję `pomnoz()`, aby pomnożyć liczby `7` i `8`.
3. Zdefiniuj funkcję `srednia_w_grupach()`, która oblicza średnią wartości dla różnych grup z danych podanych jako argumenty.

## Zadanie 11: Przekształcenia typów danych

1. Utwórz wektor `wartosci_tekstowe` zawierający elementy `"10", "20", "30"`. Zamień je na liczby za pomocą funkcji `as.numeric()`.
2. Utwórz wektor `liczby` zawierający wartości `5, 10, 15`. Zamień je na tekst za pomocą funkcji `as.character()`.
3. Utwórz wektor `odpowiedzi` zawierający elementy `"tak", "nie", "tak"`. Zamień je na czynniki (`factor`).

## Zadanie 12: Funkcje dotyczące dat i czasu

1. Sprawdź bieżącą datę i czas za pomocą funkcji `Sys.Date()` i `Sys.time()`.
2. Utwórz dwie daty: `"2024-01-01"` i `"2024-12-31"`. Oblicz liczbę dni między tymi datami.
3. Sformatuj bieżącą datę, aby wyświetlić ją w formacie `"DD-MM-YYYY"`.

## Zadanie 13: Losowanie i permutacje

1. Utwórz wektor `literki` zawierający litery `"A", "B", "C", "D", "E"`. Wylosuj z niego 2 elementy bez powtórzeń.
2. Ustaw ziarno losowania (`set.seed(42)`) i ponownie wylosuj 2 elementy z `literki`.
3. Stwórz losową permutację liczb od `1` do `6`.
