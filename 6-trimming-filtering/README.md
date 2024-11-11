# Przycinanie (trimming) i filtrowanie danych NGS

## Przycinanie odczytów na podstawie jakości

Przycinanie odczytów na podstawie jakości jest kluczowym etapem w przetwarzaniu danych sekwencjonowania następnej generacji (NGS). Wartości jakości (np. Phred score) przypisywane do każdej bazy informują o wiarygodności odczytu tej bazy przez sekwenator. Bazy na końcach odczytów często charakteryzują się niższą jakością ze względu na zmniejszoną wydajność enzymów lub problemy z detekcją sygnału. Usunięcie tych niskiej jakości baz poprzez przycinanie zwiększa dokładność późniejszych analiz, takich jak mapowanie odczytów do genomu referencyjnego czy składanie de novo.

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
