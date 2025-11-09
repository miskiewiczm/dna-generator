# 🧬 Deterministyczny Generator Sekwencji DNA

Zaawansowany generator sekwencji DNA z kontrolą jakości biochemicznej i obsługą trybu deterministycznego oraz losowego. Generator wykorzystuje algorytm backtrackingu do tworzenia wysokiej jakości sekwencji DNA spełniających zadane kryteria.

**🔗 Integracja z dna_commons:** Ten projekt wykorzystuje wspólną bibliotekę `dna_commons` do walidacji, analizy sekwencji i obliczeń termodynamicznych, zapewniając spójność z innymi projektami DNA.

## ✨ Funkcjonalności

- **Tryb deterministyczny** - powtarzalne wyniki z tym samym seedem
- **Tryb losowy** - różnorodne sekwencje za każdym razem
- **Kontrola jakości** - walidacja zawartości GC, temperatury topnienia, struktur wtórnych
- **Algorytm backtrackingu** - inteligentne generowanie spełniające kryteria biochemiczne
- **Modularna architektura** - łatwa rozszerzalność i konfiguracja
- **Szczegółowa analiza** - metryki jakości, statystyki generowania
- **Interface linii poleceń** - wygodne użycie z terminala
- **Obsługa primer3** - obliczenia termodynamiczne struktur DNA

## 🚀 Szybki start

### Instalacja

```bash
# Wymagania: Python 3.8+
# Zainstaluj zależności termodynamiczne (opcjonalnie)
pip install primer3-py

# Struktura katalogów:
# - deterministic_generator/    # Ten projekt
# - dna_commons/               # Współdzielona biblioteka (wymagana)
# - dna_encoder_refactored/    # Opcjonalnie, dla kompatybilności
```

**Uwaga:** Projekt wymaga biblioteki `dna_commons` w katalogu nadrzędnym dla pełnej funkcjonalności.

### Podstawowe użycie

```python
from deterministic_generator import DNAGenerator, GeneratorConfig, GenerationMode

# Tryb deterministyczny z zalecana sekwencją startową
config = GeneratorConfig(generation_mode=GenerationMode.DETERMINISTIC, seed=12345)
generator = DNAGenerator(config)
result = generator.generate("CCTGTCATCACGCTAGTAAC", 100)

if result.success:
    print(f"Wygenerowana sekwencja: {result.sequence}")
    print(f"Zawartość GC: {result.quality_metrics.gc_content:.1%}")
    print(f"Temperatura topnienia: {result.quality_metrics.melting_temperature:.1f}°C")
else:
    print(f"Błąd: {result.error_message}")
```

### Uruchomienie z linii poleceń

```bash
# Tryb deterministyczny
python -m deterministic_generator --initial CCTGTCATCACGCTAGTAAC --length 100 --mode deterministic

# Tryb losowy z 3 sekwencjami
python -m deterministic_generator --initial CCTGTCATCACGCTAGTAAC --length 100 --mode random --count 3

# Niestandardowe parametry jakości
python -m deterministic_generator --initial CCTGTCATCACGCTAGTAAC --length 100 --min-gc 0.4 --max-gc 0.6

# Zapis do pliku FASTA
python -m deterministic_generator --initial CCTGTCATCACGCTAGTAAC --length 100 --count 5 --format fasta --output sequences.fasta

# Wyłączenie heurystyk backtrackingu
python -m deterministic_generator --initial CCTGTCATCACGCTAGTAAC --length 100 --no-heuristics
```

#### Profile walidacji

Dodano profil `pcr_friendly` zoptymalizowany pod długie sekwencje używane w kontekście PCR. Utrzymuje aktywne wszystkie sprawdzenia, ale stosuje umiarkowanie szerokie zakresy GC/Tm i nieco bardziej liberalne limity struktur wtórnych, by ułatwić generowanie przy zachowaniu właściwości przyjaznych PCR.

Przykład:

```bash
python -m deterministic_generator \
  --initial CCTGTCATCACGCTAGTAAC \
  --length 200 \
  --profile pcr_friendly
```

#### Profil użytkownika (JSON)

Możesz zdefiniować własny profil walidacji w pliku JSON i użyć go przez `--profile user --profile-file`.

Przykładowy plik `my_profile.json`:

```json
{
  "rules": {
    "gc_content": true,
    "melting_temperature": true,
    "hairpin_structures": true,
    "homodimer_structures": true,
    "homopolymer_runs": true,
    "dinucleotide_repeats": true,
    "three_prime_stability": true
  },
  "params": {
    "min_gc": 0.42,
    "max_gc": 0.58,
    "min_tm": 52.0,
    "max_tm": 68.0,
    "max_hairpin_tm": 28.0,
    "max_homodimer_tm": 28.0,
    "max_homopolymer_length": 4,
    "max_dinucleotide_repeats": 3,
    "max_3prime_gc": 3
  }
}
```

Uruchomienie z własnym profilem:

```bash
python -m deterministic_generator \
  --initial CCTGTCATCACGCTAGTAAC \
  --length 100 \
  --profile user \
  --profile-file my_profile.json
```

#### Wizualizacja metryk (plot_validation) z profilem użytkownika

Skrypt wizualizacyjny czyta sekwencję ze standardowego wejścia i rysuje wykresy GC/Tm/hairpin/homodimer/homopolimer w oknach kroczących.

Przykład z własnym profilem:

```bash
echo CCTGTCATCACGCTAGTAACGATTACAGGCT | \
  python -m deterministic_generator.plot_validation \
  --window-size 25 --step 5 \
  --profile user --profile-file my_profile.json
```

Uwaga: aby widzieć Tm/hairpin/homodimer, zainstaluj `primer3-py`. Do renderowania wykresów wymagane jest `matplotlib`.

### Demo

```bash
# Uruchom demonstrację wszystkich funkcjonalności
python -m deterministic_generator.demo
```

## 🏗️ Architektura

### Modułowa struktura (po refaktoryzacji i integracji dna_commons)

System został zrefaktoryzowany dla lepszej spójności i integracji z `dna_commons`:

```
DNAGenerator (API Layer)           # Wysokopoziomowe API
    ├── BacktrackingEngine          # Algorytm backtrackingu
    └── dna_commons                 # Współdzielona biblioteka
        ├── DNAValidator            # Walidacja sekwencji
        ├── Primer3Adapter          # Obliczenia termodynamiczne
        ├── SequenceAnalyzer        # Analiza sekwencji
        ├── DeterministicRandom     # Generator deterministyczny
        └── ValidationRules         # Konfiguracja walidacji
```

#### Kluczowe komponenty:

- **`DNAGenerator`** - Główny interfejs API, zarządza procesem generowania
- **`BacktrackingEngine`** - Implementacja algorytmu backtrackingu, wydzielona dla lepszej modularności
- **`GeneratorConfig`** - Konfiguracja z profilami walidacji i auto-dostosowaniem
- **`dna_commons`** - Współdzielona biblioteka zawierająca:
  - **`DNAValidator`** - Zunifikowane sprawdzenia jakości sekwencji
  - **`Primer3Adapter`** - Adapter dla biblioteki primer3 z automatycznym fallback
  - **`SequenceAnalyzer`** - Narzędzia analizy sekwencji DNA
  - **`DeterministicRandom`** - Generator deterministyczny z seedem
  - **`ValidationRules`** - Konfiguracja reguł walidacji

### Obsługa primer3

System automatycznie wykrywa dostępność biblioteki `primer3-py`:
- **Z primer3**: Pełne obliczenia termodynamiczne (Tm, hairpin, homodimer)
- **Bez primer3**: Uproszczone obliczenia fallback, automatyczne wyłączenie sprawdzeń termodynamicznych

Sprawdź status:
```bash
python -m deterministic_generator --show-primer3-status
```

## 📖 Dokumentacja API

### Klasy główne

#### `DNAGenerator`

Główna klasa do generowania sekwencji DNA.

```python
generator = DNAGenerator(config)

# Generowanie pojedynczej sekwencji
result = generator.generate(initial_sequence, target_length, seed=None)

# Generowanie wielu sekwencji (do testowania determinizmu)
results = generator.generate_multiple(initial_sequence, target_length, count=5, seed=None)
```

#### `GeneratorConfig`

Konfiguracja parametrów generowania.

```python
config = GeneratorConfig(
    generation_mode=GenerationMode.DETERMINISTIC,  # lub GenerationMode.RANDOM
    seed=12345,                    # Seed dla trybu deterministycznego
    min_gc=0.45,                   # Minimalna zawartość GC (0.0-1.0)
    max_gc=0.55,                   # Maksymalna zawartość GC (0.0-1.0)
    min_tm=55.0,                   # Minimalna temperatura topnienia [°C]
    max_tm=65.0,                   # Maksymalna temperatura topnienia [°C]
    max_hairpin_tm=30.0,           # Maksymalna Tm struktur hairpin [°C]
    max_homodimer_tm=30.0,         # Maksymalna Tm homodimerów [°C]
    window_size=20,                # Rozmiar okna analizy jakości
    max_homopolymer_length=4,      # Maksymalna długość homopolimerów
    max_backtrack_attempts=10000   # Maksymalna liczba prób backtrackingu
)
```

#### `GenerationResult`

Wynik generowania sekwencji.

```python
result = generator.generate(...)

if result.success:
    print(f"Sekwencja: {result.sequence}")
    print(f"Długość: {result.actual_length}")
    print(f"Czas: {result.generation_time:.2f}s")
    
    # Metryki jakości
    metrics = result.quality_metrics
    print(f"GC: {metrics.gc_content:.1%}")
    print(f"Tm: {metrics.melting_temperature:.1f}°C")
    print(f"Spełnia kryteria: {metrics.is_valid}")
    
    # Statystyki generowania
    stats = result.generation_stats
    print(f"Próby backtrackingu: {stats['backtrack_count']}")
    print(f"Całkowite próby: {stats['total_attempts']}")
```

### Tryby generowania

#### `GenerationMode.DETERMINISTIC`
- Używa seeda do generowania powtarzalnych wyników
- Ten sam seed + parametry = identyczne sekwencje
- Przydatne do testów i reprodukowalnych badań

#### `GenerationMode.RANDOM`
- Używa prawdziwie losowego generatora
- Każde uruchomienie daje różne wyniki
- Przydatne do generowania zróżnicowanych bibliotek sekwencji

## 🔬 Przykłady użycia

### 1. Porównanie trybów deterministycznego i losowego

```python
from deterministic_generator import DNAGenerator, GeneratorConfig, GenerationMode

# Deterministyczny - zawsze te same wyniki
det_config = GeneratorConfig(generation_mode=GenerationMode.DETERMINISTIC, seed=42)
det_generator = DNAGenerator(det_config)

print("=== TRYB DETERMINISTYCZNY ===")
for i in range(3):
    result = det_generator.generate("ATGCTAGC", 50)
    print(f"Sekwencja {i+1}: {result.sequence}")
# Wszystkie 3 będą identyczne!

# Losowy - różne wyniki za każdym razem
rand_config = GeneratorConfig(generation_mode=GenerationMode.RANDOM)
rand_generator = DNAGenerator(rand_config)

print("\n=== TRYB LOSOWY ===")
for i in range(3):
    result = rand_generator.generate("ATGCTAGC", 50)
    print(f"Sekwencja {i+1}: {result.sequence}")
# Wszystkie 3 będą różne!
```

### 2. Niestandardowe parametry jakości

```python
# Restrykcyjne parametry dla wysokiej jakości
strict_config = GeneratorConfig(
    min_gc=0.48,
    max_gc=0.52,                    # Wąski zakres GC
    max_homopolymer_length=3,       # Krótsze homopolimery
    max_hairpin_tm=25.0,            # Niższa Tm hairpinów
    window_size=15                  # Mniejsze okno analizy
)

generator = DNAGenerator(strict_config)
result = generator.generate("CCTGTCATCACGCTAGTAAC", 80)

if result.success:
    print("Sekwencja wysokiej jakości wygenerowana!")
    print(f"GC: {result.quality_metrics.gc_content:.2%}")
    print(f"Hairpin Tm: {result.quality_metrics.hairpin_tm:.1f}°C")
else:
    print("Nie udało się wygenerować sekwencji o tak wysokich standardach")
```

### 3. Analiza jakości sekwencji

```python
from deterministic_generator import SequenceAnalyzer

# Wygeneruj sekwencję
generator = DNAGenerator()
result = generator.generate("CCTGTCATCACGCTAGTAAC", 100)

if result.success:
    # Podstawowe metryki z generatora
    print("=== METRYKI PODSTAWOWE ===")
    metrics = result.quality_metrics
    print(f"Zawartość GC: {metrics.gc_content:.2%}")
    print(f"Temperatura topnienia: {metrics.melting_temperature:.1f}°C")
    print(f"Hairpin Tm: {metrics.hairpin_tm:.1f}°C")
    print(f"Długie homopolimery: {'Tak' if metrics.has_homopolymers else 'Nie'}")
    
    # Szczegółowa analiza
    print("\n=== ANALIZA SZCZEGÓŁOWA ===")
    analyzer = SequenceAnalyzer()
    analysis = analyzer.analyze_sequence(result.sequence)
    
    # Rozkład nukleotydów
    print("Rozkład nukleotydów:")
    for nuc, pct in analysis['nucleotide_distribution'].items():
        print(f"  {nuc}: {pct:.1f}%")
    
    # Najdłuższy homopolimer
    longest = analysis['longest_homopolymer']
    print(f"Najdłuższy homopolimer: {longest[0]} x{longest[1]} (pozycja {longest[2]})")
    
    # Złożoność sekwencji
    print(f"Złożoność: {analysis['complexity']:.3f}")
```

### 4. Batch generowanie z różnymi seedami

```python
import json

def generate_sequence_library(base_sequence, length, num_sequences=10):
    """Generuje bibliotekę sekwencji z różnymi seedami."""
    
    sequences = []
    
    for i in range(num_sequences):
        config = GeneratorConfig(
            generation_mode=GenerationMode.DETERMINISTIC,
            seed=i + 1000  # Różne seedy
        )
        
        generator = DNAGenerator(config)
        result = generator.generate(base_sequence, length)
        
        if result.success:
            sequences.append({
                'id': f'seq_{i+1:03d}',
                'seed': i + 1000,
                'sequence': result.sequence,
                'gc_content': result.quality_metrics.gc_content,
                'tm': result.quality_metrics.melting_temperature,
                'generation_time': result.generation_time
            })
        
        print(f"Sekwencja {i+1}/{num_sequences}: {'✅' if result.success else '❌'}")
    
    return sequences

# Użycie
library = generate_sequence_library("CCTGTCATCACGCTAGTAAC", 100, 5)

# Zapisz do JSON
with open('sequence_library.json', 'w') as f:
    json.dump(library, f, indent=2)

print(f"Wygenerowano {len(library)} sekwencji!")
```

### 5. Walidacja i porównanie różnych algorytmów

```python
def compare_generation_strategies():
    """Porównuje różne strategie generowania."""
    
    strategies = [
        ("Standardowe", GeneratorConfig()),
        ("Wysokie GC", GeneratorConfig(min_gc=0.55, max_gc=0.65)),
        ("Niskie GC", GeneratorConfig(min_gc=0.35, max_gc=0.45)),
        ("Małe okno", GeneratorConfig(window_size=10)),
        ("Duże okno", GeneratorConfig(window_size=30)),
        ("Restrykcyjne", GeneratorConfig(
            max_homopolymer_length=3,
            max_hairpin_tm=25.0,
            window_size=15
        ))
    ]
    
    base_seq = "CCTGTCATCACGCTAGTAAC"
    target_len = 80
    
    for name, config in strategies:
        config.generation_mode = GenerationMode.DETERMINISTIC
        config.seed = 12345
        config.enable_progress_logging = False
        
        generator = DNAGenerator(config)
        result = generator.generate(base_seq, target_len)
        
        print(f"\n=== {name.upper()} ===")
        if result.success:
            print(f"Sukces w {result.generation_time:.3f}s")
            print(f"GC: {result.quality_metrics.gc_content:.1%}")
            print(f"Tm: {result.quality_metrics.melting_temperature:.1f}°C")
            
            stats = result.generation_stats
            print(f"Próby: {stats['total_attempts']}")
            print(f"Backtrack: {stats['backtrack_count']}")
        else:
            print(f"Niepowodzenie w {result.generation_time:.3f}s")

# Uruchom porównanie
compare_generation_strategies()
```

## 🔧 Parametry konfiguracji

### Parametry jakości DNA

| Parametr                   | Domyślnie | Opis                                                                         |
|----------------------------|-----------|------------------------------------------------------------------------------|
| `min_gc`                   | 0.45      | Minimalna zawartość GC (0.0-1.0)                                             |
| `max_gc`                   | 0.55      | Maksymalna zawartość GC (0.0-1.0)                                            |
| `min_tm`                   | 55.0      | Minimalna temperatura topnienia [°C]                                         |
| `max_tm`                   | 65.0      | Maksymalna temperatura topnienia [°C]                                        |
| `max_hairpin_tm`           | 30.0      | Maksymalna Tm struktur hairpin [°C]                                          |
| `max_homodimer_tm`         | 30.0      | Maksymalna Tm homodimerów [°C]                                               |
| `max_homopolymer_length`   | 4         | Maksymalna dozwolona długość kolejnych identycznych nukleotydów (<= wartość) |
| `max_dinucleotide_repeats` | 2         | Maksymalna liczba kolejnych powtórzeń tego samego dinukleotydu (<= wartość)  |
| `max_3prime_gc`            | 3         | Max G/C w ostatnich 5 nukleotydach                                           |

### Parametry algorytmu

| Parametr                      | Domyślnie     | Opis                                                       |
|-------------------------------|---------------|------------------------------------------------------------|
| `window_size`                 | 20            | Rozmiar okna analizy jakości                               |
| `max_backtrack_attempts`      | 10000         | Limit prób backtrackingu                                   |
| `generation_mode`             | DETERMINISTIC | Tryb generowania                                           |
| `seed`                        | None          | Seed (None = auto w trybie deterministycznym)              |
| `enable_backtrack_heuristics` | True          | Włącza lekkie heurystyki wyboru baz (szybsza konwergencja) |

### Parametry primer3

| Parametr            | Domyślnie | Opis                                   |
|---------------------|-----------|----------------------------------------|
| `primer3_mv_conc`   | 50.0      | Stężenie jonów jednowartościowych [mM] |
| `primer3_dv_conc`   | 4.0       | Stężenie jonów dwuwartościowych [mM]   |
| `primer3_dntp_conc` | 0.5       | Stężenie dNTP [mM]                     |
| `primer3_dna_conc`  | 50.0      | Stężenie DNA [nM]                      |

## 📊 Formaty wyjściowe

### JSON
```bash
# Pretty JSON
python -m deterministic_generator --initial ATGC --length 50 --format json --output results.json

# Surowe metryki w JSON (liczby)
python -m deterministic_generator --initial ATGC --length 50 --format json --json-raw-metrics --output results.json
```

### FASTA  
```bash
python -m deterministic_generator --initial ATGC --length 50 --count 5 --format fasta --output sequences.fasta
```

### Text (domyślny)
```bash
python -m deterministic_generator --initial ATGC --length 50 --verbose

# Tylko sekwencje (każda w osobnej linii)
python -m deterministic_generator --initial ATGC --length 50 --sequences-only
```

> Uwaga: Walidacja podczas generowania dotyczy okien o rozmiarze `window_size` (kryteria egzekwowane „po drodze”).

Metryki całej sekwencji (np. Tm) są informacyjne i mogą wykraczać poza zakresy okienne dla długich sekwencji.

## 🏗️ Architektura

```sh
deterministic_generator/
├── __init__.py          # API pakietu (z adapterami do dna_commons)
├── config.py           # Konfiguracja i parametry
├── generator.py        # Główna klasa DNAGenerator
├── backtracking_engine.py # Algorytm backtrackingu
├── exceptions.py      # Dedykowane wyjątki
├── __main__.py        # Interface linii poleceń
├── plot_validation.py # Wizualizacja metryk walidacji
├── demo.py           # Demonstracja funkcjonalności
└── README.md         # Ten plik

../dna_commons/         # Współdzielona biblioteka
├── validation/        # Moduły walidacji DNA
├── thermodynamics/    # Obliczenia termodynamiczne
├── sequence/         # Analiza sekwencji
├── random/           # Generatory deterministyczne
└── config/           # Profile konfiguracji
```

## 🤝 Porównanie z oryginalnym kodem

| Aspekt           | Oryginalny kod         | Nowa implementacja            |
|------------------|------------------------|-------------------------------|
| **Struktura**    | Monolityczny skrypt    | Modularna architektura        |
| **Konfiguracja** | Zmienne globalne       | Klasa GeneratorConfig         |
| **Tryby**        | Tylko deterministyczny | Deterministyczny + losowy     |
| **API**          | Funkcje                | Klasy z metodami              |
| **Walidacja**    | Inline                 | Dedykowana klasa DNAValidator |
| **Wyniki**       | Print do konsoli       | Strukturowane obiekty         |
| **CLI**          | Brak                   | Pełny interface linii poleceń |
| **Testy**        | Podstawowe             | Demo + przykłady              |
| **Analiza**      | Ograniczona            | Szczegółowa analiza sekwencji |

## 🔬 Algorytm backtrackingu

Generator wykorzystuje zaawansowany algorytm backtrackingu:

1. **Inicjalizacja** - startuje od sekwencji początkowej
2. **Eksploracja** - dodaje kolejne nukleotydy (A, T, G, C)
3. **Walidacja** - sprawdza jakość w oknie analizy
4. **Akceptacja/Odrzucenie** - akceptuje dobre kandydaty
5. **Backtracking** - cofa się gdy nie ma dobrych opcji
6. **Kontynuacja** - powtarza aż do osiągnięcia długości docelowej

### Kryteria walidacji (w kolejności sprawdzania):
- Brak długich homopolimerów
- Brak powtórzeń dinukleotydów  
- Stabilność końca 3' (max 3 G/C w ostatnich 5 nt)
- Zawartość GC w zadanym zakresie
- Temperatura topnienia w zadanym zakresie
- Niska temperatura struktur hairpin
- Niska temperatura homodimerów

## 💡 Wskazówki użytkowania

### Optymalizacja wydajności
- Użyj większego `window_size` dla lepszej jakości, ale wolniejszego generowania
- Zwiększ `max_backtrack_attempts` dla trudniejszych przypadków
- Rozluźnij parametry jakości jeśli generowanie się nie powodzi

### Najlepsze praktyki
- Używaj znanych dobrych sekwencji początkowych (np. z oryginalnego kodu)
- Testuj różne seedy w trybie deterministycznym
- Zapisuj udane konfiguracje do ponownego użycia
- Analizuj niepowodzenia walidacji z `generation_stats`

### Rozwiązywanie problemów
- **"Wyczerpano wszystkie możliwości backtrackingu"** → Zwiększ `max_backtrack_attempts` lub rozluźnij kryteria
- **Niska jakość sekwencji** → Zmniejsz `window_size` lub dostosuj parametry jakości  
- **Wolne generowanie** → Użyj mniejszego `window_size` lub mniej restrykcyjnych parametrów
 - **Brak primer3-py** → Wyłącz sprawdzenia termodynamiczne (np. profil `sequence_only` lub flagi `--no-tm-check`, `--no-hairpin-check`, `--no-homodimer-check`) albo zainstaluj `primer3-py`

## 📌 Uwagi dot. determinizmu
- Determinizm dotyczy wyborów generatora pseudolosowego (seed i logika backtrackingu).
- Wyniki zależne od `primer3` (Tm/hairpin/homodimer) mogą różnić się między wersjami/systemami. Dla pełnej reprodukowalności utrzymuj tę samą wersję `primer3-py` i zależności.

## 📄 Licencja

Ten projekt jest rozszerzeniem oryginalnego kodu deterministycznego generatora DNA.

---

**🧬 Happy DNA generating! 🧬**
