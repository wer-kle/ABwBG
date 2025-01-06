# Kolumny metadanych GRanges w `rowRanges(vcf)`

**GRanges** to obiekt w R, który służy do przechowywania informacji o regionach genomu (np. wariantów genetycznych). Zawiera m.in. następujące kolumny metadanych:
- **paramRangeID**: Identyfikator zakresu pozycji w genomie, zwykle określający region, z którego pochodzi wariant (np. używany w filtrach lub obliczeniach zakresów).
- **REF**: Referencyjna sekwencja nukleotydów w tym miejscu w genomie (np. „A” lub „AT”).
- **ALT**: Alternatywne allele, które mogą wystąpić zamiast sekwencji referencyjnej (np. „G” lub „-” w przypadku delecji).
- **QUAL**: Jakość wywołania wariantu, wartość liczbowa, która wskazuje na wiarygodność identyfikacji wariantu.
- **FILTER**: Informacja o zastosowanych filtrach (np. `PASS` oznacza, że wariant spełnia wszystkie kryteria; inne wartości mogą wskazywać problemy, jak niski wskaźnik jakości).

---

# 2. Kolumny `info(vcf)`

`info(vcf)` to tabela zawierająca szczegółowe adnotacje dla każdego wariantu. Każda kolumna opisuje specyficzne atrybuty wariantu:
- **LDAF**: Allele frequency (częstość allelu) wyliczona w modelu ukierunkowanym (LD).
- **AVGPOST**: Średnie posterior probability (prawdopodobieństwo a posteriori) dla wariantu.
- **RSQ**: Współczynnik determinacji (często używany do oceny jakości imputacji wariantów).
- **ERATE**: Szacowana szybkość błędu (error rate).
- **THETA**: Parametr modelu demograficznego, opisujący zmienność.
- **CIEND**: Przedział ufności dla końca wariantu strukturalnego (SV).
- **CIPOS**: Przedział ufności dla pozycji początkowej wariantu.
- **END**: Pozycja końcowa wariantu (szczególnie przy wariantach strukturalnych).
- **HOMLEN**: Długość regionu homologicznego wokół miejsca wariantu.
- **HOMSEQ**: Sekwencja homologiczna.
- **SVLEN**: Długość wariantu strukturalnego (np. delecji, insercji).
- **SVTYPE**: Typ wariantu strukturalnego (np. `DEL` dla delecji, `INS` dla insercji).
- **AC**: Liczba alleli alternatywnych w populacji (ang. *Allele Count*).
- **AN**: Całkowita liczba zliczonych alleli (ang. *Allele Number*).
- **AA**: Ancestral allele (allel pierwotny, obecny w przodkach).
- **AF**: Allele frequency (częstość allelu alternatywnego w populacji).

... i inne kolumny, które zazwyczaj są specyficzne dla projektu lub analizy.

---

# 3. Elementy listy `geno(vcf)`

`geno(vcf)` to lista danych genotypowych przypisanych do każdego próbki:
- **GT (Genotype)**: Genotyp, np. `0/1` (heterozygota) lub `1/1` (homozygota alternatywna).
- **DS (Dosage)**: Dawka allelu alternatywnego, np. wartość ciągła wskazująca, ile kopii alternatywnego allelu przypisano danej próbce (używane w analizie asocjacji).
- **GL (Genotype Likelihood)**: Prawdopodobieństwa genotypu, wyrażone jako logarytmiczne wartości (log10), np. `-10, -2, -0.5` dla genotypów `0/0`, `0/1`, `1/1`.
