# Zadania treningowe z R base
***Instrukcja:*** **Zapisz rozwiązania zadań z każdej części jako osobny skrypt w formacie .Rmd lub (ewentualnie) .R, np. Czesc1.Rmd, Czesc2.R. Pamiętaj, że im dokładniej opisane zadanie, tym lepiej!**

## Część 1: Podstawowe operacje matematyczne

1. **Zadanie 1:** Oblicz pole trójkąta o podstawie 10 i wysokości 5, używając operatora mnożenia.

2. **Zadanie 2:** Znajdź resztę z dzielenia 256 przez 7 oraz wynik dzielenia całkowitego.

3. **Zadanie 3:** Oblicz wartość wyrażenia $e^{\ln(5)}$, używając funkcji `exp()` i `log()`.

4. **Zadanie 4:** Zaokrąglij liczbę 7,895 do najbliższej liczby całkowitej: w górę i w dół.

5. **Zadanie 5:** Oblicz wartość funkcji tangens dla kąta $45^\circ$ (pamiętaj o konwersji stopni na radiany).

6. **Zadanie 6:** Oblicz wartość wyrażenia $\sqrt{3^2 + 4^2}$, używając funkcji `sqrt()` i operatora potęgowania.

---

## Część 2: Funkcje statystyczne

1. **Zadanie 1:** Dla wektora danych $[2, 4, 6, 8, 10]$ oblicz średnią arytmetyczną, medianę oraz wariancję.

2. **Zadanie 2:** Stwórz wektor z 50 losowych liczb z przedziału od 1 do 100. Znajdź wartość minimalną, maksymalną oraz sumę wszystkich elementów.

3. **Zadanie 3:** Dla wektora danych $[3, 7, 7, 7, 2, 2, 5]$ utwórz tabelę częstości występowania każdej wartości.

4. **Zadanie 4:** Oblicz korelację i kowariancję między danymi $x = [1, 3, 5, 7, 9]$ i $y = [2, 6, 10, 14, 18]$.

5. **Zadanie 5:** Użyj funkcji `summary()` do podsumowania danych z ramki danych zawierającej informacje o wieku i wysokości pięciu osób.

6. **Zadanie 6:** Oblicz iloczyn wszystkich liczb w wektorze $[1, 2, 3, 4, 5]$.

---

## Część 3: Operacje na wektorach i indeksowanie

1. **Zadanie 1:** Stwórz wektor, który będzie zwracał wszystkie wartości od 0 do 1 co 0,1.

2. **Zadanie 2:** Dla wektora $[2, 4, 6, 8, 10]$ wybierz drugi i czwarty element.

3. **Zadanie 3:** Znajdź wszystkie elementy wektora $[5, 3, 8, 1, 9]$ większe niż 4.

4. **Zadanie 4:** Posortuj wektor $[5, 2, 8, 3, 7]$ malejąco i podaj indeksy sortowania.

5. **Zadanie 5:** Oblicz rangi elementów wektora $[50, 20, 30, 40, 10]$.

6. **Zadanie 6:** Stwórz wektor powtarzający sekwencję $[1, 2, 3]$ trzy razy.

---

## Część 4: Praca z ramkami danych (data frames)

1. **Zadanie 1:** Utwórz ramkę danych zawierającą informacje o produktach: nazwa (tekst), cena (liczba), ilość (liczba całkowita).

2. **Zadanie 2:** Wyświetl kolumnę `'nazwa'` z ramki danych produktów.

3. **Zadanie 3:** Dodaj nową kolumnę `'wartość'` obliczaną jako cena pomnożona przez ilość.

4. **Zadanie 4:** Usuń kolumnę `'ilość'` z ramki danych.

5. **Zadanie 5:** Wyświetl strukturę ramki danych za pomocą funkcji `str()`.

6. **Zadanie 6:** Podaj nazwy kolumn i wierszy w ramce danych.

---

## Część 5: Funkcje logiczne i warunkowe

1. **Zadanie 1:** Sprawdź, czy liczba 7 jest większa od 5 i jednocześnie mniejsza od 10.

2. **Zadanie 2:** Dla wektora $[-3, 2, 5]$ użyj funkcji `ifelse()`, aby przypisać "Dodatnia" lub "Ujemna" do każdego elementu.

3. **Zadanie 3:** Sprawdź, czy jakikolwiek element wektora $[FALSE, FALSE, TRUE]$ jest prawdziwy.

4. **Zadanie 4:** Znajdź indeksy liczb parzystych w wektorze $[1, 2, 3, 4, 5, 6]$.

5. **Zadanie 5:** Napisz instrukcję `if...else`, która sprawdzi, czy dana liczba jest podzielna przez 3.

6. **Zadanie 6:** Użyj operatora NOT, aby odwrócić wynik porównania $x == y$.

---

## Część 6: Pętle i iteracje

1. **Zadanie 1:** Napisz pętlę `for`, która wypisze liczby od 1 do 5.

2. **Zadanie 2:** Użyj pętli `while`, aby obliczyć silnię liczby 5.

3. **Zadanie 3:** Stwórz macierz $2 \times 5$ wypełnioną kolumnami i za pomocą funkcji `apply()` oblicz sumę elementów w każdym wierszu.

4. **Zadanie 4:** Użyj funkcji `sapply()` na liście $[1\!:\!3,\; 4\!:\!6,\; 7\!:\!9]$ do obliczenia sumy elementów każdej podlisty.

5. **Zadanie 5:** Napisz pętlę `repeat`, która będzie zwiększać zmienną $x$ o 2, aż $x$ przekroczy 10.

6. **Zadanie 6:** Użyj funkcji `tapply()` do obliczenia średniej długości słów w grupach tekstów.

---

## Część 7: Funkcje wejścia/wyjścia

1. **Zadanie 1:** Załaduj ramkę danych z pliku `"studenci.csv"` i wyświetl jej strukturę.

2. **Zadanie 2:** Zapisz ramkę danych do pliku `"dane_output.csv"` z nazwami wierszy.

3. **Zadanie 3:** Użyj funkcji `readLines()`, aby wczytać pierwsze 10 linii z pliku `"dokument.txt"`.

4. **Zadanie 4:** Napisz wektor $["Linia1", "Linia2", "Linia3"]$ do pliku `"linie.txt"` za pomocą `writeLines()`.

5. **Zadanie 5:** Sprawdź, czy plik `"dane.csv"` istnieje w bieżącym katalogu roboczym.

6. **Zadanie 6:** Ustaw katalog roboczy na `"C:/Użytkownicy/Jan/Dokumenty"`.

---

## Część 8: Funkcje związane z łańcuchami znaków

1. **Zadanie 1:** Użyj funkcji `paste()`, aby utworzyć zdanie "R jest wspaniały!".

2. **Zadanie 2:** Znajdź długość łańcucha znaków "Analiza danych".

3. **Zadanie 3:** Zamień wszystkie litery w tekście "Witaj Świecie" na duże litery.

4. **Zadanie 4:** Sprawdź, które elementy wektora $["jabłko", "banan", "gruszka"]$ zawierają literę "a".

5. **Zadanie 5:** W tekście "To jest stary komputer" zamień "stary" na "nowy", używając funkcji `sub()`.

6. **Zadanie 6:** Podziel zdanie "Uczę się języka R" na słowa.

---

## Część 9: Funkcje statystyczne i probabilistyczne

1. **Zadanie 1:** Wygeneruj 50 losowych liczb z rozkładu normalnego o średniej 10 i odchyleniu standardowym 2.

2. **Zadanie 2:** Oblicz wartość dystrybuanty rozkładu normalnego w punkcie $x = 1$ dla średniej 0 i odchylenia 1.

3. **Zadanie 3:** Znajdź wartość kwantyla $2{,}5\%$ rozkładu normalnego standardowego.

4. **Zadanie 4:** Przeprowadź jednostronny test t-Studenta sprawdzający, czy średnia próby jest większa od 0.

5. **Zadanie 5:** Narysuj histogram dla wygenerowanych danych z rozkładu normalnego.

6. **Zadanie 6:** Oblicz gęstość rozkładu normalnego w punktach od $-3$ do $3$ z krokiem $0{,}5$.

---

## Część 10: Funkcje systemowe i diagnostyczne

1. **Zadanie 1:** Wyświetl bieżący katalog roboczy.

2. **Zadanie 2:** Zmień katalog roboczy na folder "Projekt".

3. **Zadanie 3:** Użyj funkcji `help.search()`, aby znaleźć informacje o funkcji mediany.

4. **Zadanie 4:** Wyświetl przykłady użycia funkcji `summary()`.

5. **Zadanie 5:** Usuń wszystkie obiekty z bieżącego środowiska.

6. **Zadanie 6:** Zmierz czas wykonania funkcji `sort()` dla wektora 1 000 000 losowych liczb.

---

## Część 11: Tworzenie własnych funkcji

1. **Zadanie 1:** Zdefiniuj funkcję `kwadrat`, która podnosi liczbę do kwadratu.

2. **Zadanie 2:** Napisz funkcję `suma_wektora`, która oblicza sumę elementów wektora.

3. **Zadanie 3:** Utwórz funkcję `maksimum`, która zwraca największy z trzech podanych argumentów.

4. **Zadanie 4:** Napisz funkcję `czy_podzielna`, która sprawdza, czy liczba jest podzielna przez 3 i 5.

5. **Zadanie 5:** Stwórz funkcję `przeksztalc`, która zamienia wszystkie litery w tekście na małe.

6. **Zadanie 6:** Zdefiniuj funkcję `pole_trojkata`, która oblicza pole trójkąta na podstawie długości boków, używając wzoru Herona.

---

## Część 12: Przekształcenia typów danych

1. **Zadanie 1:** Przekształć wektor znakowy $["10", "20", "30"]$ na liczbowy i oblicz ich średnią.

2. **Zadanie 2:** Zamień wektor liczbowy $[5, 10, 15]$ na czynnik.

3. **Zadanie 3:** Sprawdź, jaki typ danych ma wynik funkcji `as.character(123)`.

4. **Zadanie 4:** Dla wektora $["Prawda", "Fałsz"]$ spróbuj przekonwertować go na wartości logiczne.

5. **Zadanie 5:** Przekształć macierz $2 \times 2$ w ramkę danych.

6. **Zadanie 6:** Sprawdź, czy wektor $[1, 2, "3"]$ jest liczbowy po konwersji.

---

## Część 13: Funkcje dotyczące dat i czasu

1. **Zadanie 1:** Wyświetl aktualną datę i czas.

2. **Zadanie 2:** Zamień tekst "01/01/2023" na datę w formacie R.

3. **Zadanie 3:** Oblicz liczbę dni między dzisiejszą datą a Twoimi urodzinami.

4. **Zadanie 4:** Sformatuj bieżący czas w formacie "godzina:minuta:sekunda".

5. **Zadanie 5:** Sprawdź, jaki dzień tygodnia przypada na datę "2023-11-11".

6. **Zadanie 6:** Dodaj 30 dni do dzisiejszej daty i wyświetl wynik.

---

## Część 14: Losowanie i permutacje

1. **Zadanie 1:** Wylosuj 6 liczb z przedziału od 1 do 49 (jak w Lotto).

2. **Zadanie 2:** Ustaw ziarno na 123 i wygeneruj losowy wektor 5 liczb całkowitych od 1 do 100.

3. **Zadanie 3:** Wymieszaj kolejność liter w słowie "PROGRAM".

4. **Zadanie 4:** Wylosuj z powtórzeniami 4 elementy z wektora $["X", "Y", "Z"]$.

5. **Zadanie 5:** Porównaj dwa losowania z tym samym ziarnem i upewnij się, że wyniki są identyczne.

6. **Zadanie 6:** Stwórz losową permutację wektora liczbowego od 1 do 20.

---

## Część 15: Inne przydatne funkcje

1. **Zadanie 1:** Wygeneruj sekwencję liczb od 1 do 100 co 5.

2. **Zadanie 2:** Użyj funkcji `sample()`, aby stworzyć losową permutację liczb od 1 do 10.

3. **Zadanie 3:** Zastosuj funkcję `length()` do każdego elementu listy $["kot", "pies", "słoń"]$, używając `lapply()`.

4. **Zadanie 4:** Oblicz średnią wartość dla każdej kolumny w ramce danych za pomocą `sapply()`.

5. **Zadanie 5:** Stwórz listę zawierającą trzy wektory liczbowe i oblicz sumę elementów każdego z nich.

6. **Zadanie 6:** Wygeneruj sekwencję liczb od 10 do 1, używając funkcji `seq()`.

---

## Część 16: Profilowanie i debugowanie

1. **Zadanie 1:** Włącz debugowanie funkcji `suma_wektora` i sprawdź, jak działa z błędnymi danymi.

2. **Zadanie 2:** Użyj funkcji `traceback()` po otrzymaniu błędu w skrypcie.

3. **Zadanie 3:** Ustaw opcję `error` na `recover` i zobacz, co się stanie po wystąpieniu błędu.

4. **Zadanie 4:** Wyłącz debugowanie funkcji `suma_wektora`.

5. **Zadanie 5:** Znajdź i popraw błąd w funkcji, która nie zwraca oczekiwanego wyniku.

6. **Zadanie 6:** Użyj funkcji `options()`, aby zresetować ustawienia błędów do domyślnych.

---
