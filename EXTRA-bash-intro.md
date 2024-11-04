# Wprowadzenie do Bash

Bash (Bourne Again SHell) to jeden z najpopularniejszych interpreterów poleceń używanych w systemach operacyjnych typu Unix. Pozwala na wykonywanie poleceń, skryptów oraz automatyzację zadań. Ten tutorial ma na celu wprowadzenie Cię w podstawy pracy z Bash, nawet jeśli nigdy wcześniej nie korzystałeś z terminala.

## Spis treści

1. **Wprowadzenie do terminala**
2. **Instalacja i konfiguracja**
   - Windows
   - Mac
3. **Podstawowe polecenia**
4. **Nawigacja po systemie plików**
5. **Operacje na plikach i katalogach**
6. **Przekierowania i potoki**
7. **Zmienne i skrypty**
8. **Zadania praktyczne**

---

## 1. Wprowadzenie do terminala

Terminal to interfejs tekstowy pozwalający na komunikację z systemem operacyjnym poprzez wprowadzanie poleceń. Jest niezwykle potężnym narzędziem dla programistów i administratorów systemów.

## 2. Instalacja i konfiguracja

### Windows

W systemie Windows domyślnie nie ma środowiska Bash. Możesz jednak zainstalować je na kilka sposobów:

#### Opcja 1: Git Bash (zalecane)

1. **Pobierz Git dla Windows** ze strony [git-scm.com](https://git-scm.com/download/win).
2. Uruchom instalator i podczas instalacji zaznacz opcję instalacji Git Bash.
3. Po zakończeniu instalacji możesz uruchomić Git Bash z menu Start.

#### Opcja 2: Windows Subsystem for Linux (WSL)

1. **Włącz WSL**:
   - Otwórz PowerShell jako administrator.
   - Wpisz: `wsl --install` i naciśnij Enter.
2. **Zainstaluj dystrybucję Linux** (np. Ubuntu) z Microsoft Store.
3. Uruchom zainstalowaną dystrybucję z menu Start.

### Mac

MacOS jest systemem Unixowym, więc Bash jest dostępny domyślnie.

1. **Otwórz Terminal**:
   - Przejdź do `Aplikacje` > `Narzędzia` > `Terminal`.
   - Alternatywnie użyj Spotlight (Cmd + Spacja) i wpisz "Terminal".

## 3. Podstawowe polecenia

- `pwd` – wyświetla bieżący katalog.
- `ls` – listuje pliki i katalogi.
- `cd` – zmienia katalog.
- `mkdir` – tworzy nowy katalog.
- `touch` – tworzy nowy plik.
- `cp` – kopiuje pliki lub katalogi.
- `mv` – przenosi lub zmienia nazwy plików/katalogów.
- `rm` – usuwa pliki.
- `rmdir` – usuwa katalogi.

## 4. Nawigacja po systemie plików

- **Przejście do katalogu domowego**: `cd ~` lub samo `cd`.
- **Przejście do katalogu nadrzędnego**: `cd ..`.
- **Przejście do konkretnego katalogu**: `cd /ścieżka/do/katalogu`.

## 5. Operacje na plikach i katalogach

- **Tworzenie nowego katalogu**: `mkdir nowy_katalog`.
- **Tworzenie nowego pliku**: `touch nowy_plik.txt`.
- **Kopiowanie pliku**: `cp plik.txt kopia_plik.txt`.
- **Przenoszenie/zmiana nazwy pliku**: `mv plik.txt nowa_nazwa.txt`.
- **Usuwanie pliku**: `rm plik.txt`.
- **Usuwanie katalogu**: `rmdir pusty_katalog` lub `rm -r katalog`.

## 6. Przekierowania i potoki

- **Przekierowanie wyjścia do pliku**: `polecenie > plik.txt`.
- **Dodanie wyjścia do istniejącego pliku**: `polecenie >> plik.txt`.
- **Przekierowanie wejścia z pliku**: `polecenie < plik.txt`.
- **Łączenie poleceń**: `polecenie1 | polecenie2`.

## 7. Zmienne i skrypty

- **Tworzenie zmiennej**: `zmienna="wartość"`.
- **Wyświetlanie zmiennej**: `echo $zmienna`.
- **Tworzenie skryptu**:
  1. Utwórz plik skryptu, np. `touch skrypt.sh`.
  2. Dodaj prawa wykonywania: `chmod +x skrypt.sh`.
  3. Edytuj plik i dodaj polecenia.
  4. Uruchom skrypt: `./skrypt.sh`.

## 8. Zadania praktyczne

### Zadanie 1: Nawigacja

1. **Otwórz terminal**.
2. **Sprawdź swój bieżący katalog**: użyj `pwd`.
3. **Przejdź do katalogu domowego**: `cd ~`.
4. **Utwórz nowy katalog o nazwie `projekt`**: `mkdir projekt`.
5. **Wejdź do katalogu `projekt`**: `cd projekt`.

### Zadanie 2: Operacje na plikach

1. **W katalogu `projekt` utwórz plik `README.md`**: `touch README.md`.
2. **Sprawdź zawartość katalogu**: `ls`.
3. **Zmień nazwę pliku na `INSTRUKCJA.md`**: `mv README.md INSTRUKCJA.md`.
4. **Utwórz kopię pliku**: `cp INSTRUKCJA.md KOPIA_INSTRUKCJA.md`.
5. **Usuń kopię pliku**: `rm KOPIA_INSTRUKCJA.md`.

### Zadanie 3: Przekierowania i potoki

1. **Wyświetl listę plików i zapisz ją do pliku `lista.txt`**: `ls > lista.txt`.
2. **Wyświetl zawartość pliku `lista.txt`**: `cat lista.txt`.
3. **Dodaj aktualną datę do pliku `lista.txt`**: `date >> lista.txt`.
4. **Wyświetl zawartość pliku ponownie**: `cat lista.txt`.

### Zadanie 4: Skrypt Bash

1. **Wróć do katalogu domowego**: `cd ~`.
2. **Utwórz plik skryptu `witaj.sh`**: `touch witaj.sh`.
3. **Dodaj prawa wykonywania**: `chmod +x witaj.sh`.
4. **Edytuj plik i dodaj następujące linie**:
   ```bash
   #!/bin/bash
   echo "Witaj, $(whoami)! Dzisiaj jest $(date)."
   ```
5. **Zapisz plik i zamknij edytor**.
6. **Uruchom skrypt**: `./witaj.sh`.

### Zadanie 5: Zmienne środowiskowe

1. **Wyświetl listę zmiennych środowiskowych**: `printenv`.
2. **Sprawdź wartość zmiennej `HOME`**: `echo $HOME`.
3. **Utwórz własną zmienną**: `moje_imie="Jan"`.
4. **Wyświetl wartość swojej zmiennej**: `echo $moje_imie`.

---
