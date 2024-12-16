1. Wczytaj i załaduj poniższe pakiety:
   - `rrBLUP`
   - `BGLR`
   - `DT`
   - `SNPRelate`
   - `dplyr`
   - `qqman`
   - `poolr`

2. Wczytaj dane genotypowe z pliku `.ped`, informacje o subpopulacjach z pliku `.fam` i informacje o mapowaniu markerów z pliku `.map`.
- Przekoduj wartości markerów zgodnie z poniższym schematem:
  - **2** → **NA** (brakujące dane).  
  - **0** → **0** (homozygota dla allelu głównego).  
  - **1** → **1** (heterozygota).  a
  - **3** → **2** (homozygota dla allelu mniejszościowego).  
   - Przekonwertuj dane na macierz i dokonaj jej transpozycji.
   - Podaj wymiary macierzy (liczbę osobników i markerów SNP).

3. Wczytaj dane fenotypowe.
   - Wyodrębnij pierwszą cechę.
   - Dopasuj dane fenotypowe do danych genotypowych.
   - Usuń brakujące dane z macierzy.

4. Przeprowadź kontrolę jakości (QC) danych markerowych:
   - Zastąp brakujące dane markerowe średnią wartością dla każdego markera.
   - Odfiltruj markery z MAF < 5%. 
   - Zaktualizuj plik `.map` i podaj nowe wymiary danych genotypowych oraz informacji o markerach.

5. Wykonaj analizę PCA:
   - Utwórz macierz markerów.
   - Przeprowadź PCA.
   - Wczytaj dodatkowe informacje o próbkach z pliku `gerplasm.csv` i przeprowadź PCA z podziałem na subpopulacje.
   - Zinterpretuj i narysuj wykres dla pierwszych dwóch składowych głównych.

6. Przygotuj dane do analizy GWAS:
   - Przygotuj dane genotypowe i fenotypowe do analizy.
   - Wykonaj analizę GWAS.

7. Skoryguj wyniki GWAS dla testów wielokrotnych:
   - Oblicz efektywną liczbę testów niezależnych (Meff) metodą Li i Ji.
   - Wyznacz odpowiedni próg istotności.

8. Wyodrębnij istotne markery SNP.
   - Podaj listę markerów SNP spełniających ustalone kryterium p-wartości.

9. Stwórz wykres Manhattan:
   - Wykorzystaj pakiet `qqman`.
   - Uwzględnij wyróżnienie istotnych markerów.

10. Zinterpretuj wyniki:
    - Jakie markery SNP są potencjalnie istotne?
    - Jakie wnioski można wyciągnąć na podstawie analizy PCA i wyników GWAS?
