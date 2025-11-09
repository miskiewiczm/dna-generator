# 📋 TODO - Ulepszenia Generatora Sekwencji DNA

## ✅ **Status Implementacji**

### **Problem 2: Elastyczne sprawdzenia - ROZWIĄZANY** ✅
**Data rozwiązania:** 3 września 2025

**Co zostało zaimplementowane:**
- Profile walidacji: `strict`, `relaxed`, `minimal`, `sequence_only`, `thermodynamic_only`, `none`
- Flagi CLI do wyłączania konkretnych sprawdzeń: `--no-gc-check`, `--no-tm-check`, itd.
- Konfiguracja `validation_profile` i `validation_rules` w `GeneratorConfig`
- Modyfikacja metody `_validate_window_with_stats()` do respektowania flag
- Przykłady użycia w demo.py (demo 5 i 6)

**Jak używać:**
```bash
# Z profilem walidacji
python -m deterministic_generator --initial CCTGTCATCACGCTAGTAAC --length 100 --profile relaxed

# Z wyłączeniem konkretnych sprawdzeń
python -m deterministic_generator --initial CCTGTCATCACGCTAGTAAC --length 100 --no-hairpin-check --no-homodimer-check

# Z profilem 'none' (brak walidacji)
python -m deterministic_generator --initial CCTGTCATCACGCTAGTAAC --length 100 --profile none
```

## 🚨 **Zidentyfikowane Problemy**

### **Problem 1: Konflikt parametrów sekwencja początkowa vs docelowa**

**Opis problemu:**
- Sekwencja początkowa może mieć parametry niezgodne z docelowymi kryteriami
- Przykład: sekwencja początkowa z GC=100%, docelowy zakres GC=20-40%
- Algorytm backtrackingu nie może znaleźć rozwiązania w oknie analizy
- Powód: pierwsze N nukleotydów zawsze będzie miało niewłaściwe parametry

**Przykład konfliktowej sytuacji:**
```
Sekwencja początkowa: CCCGGGCCCGGG (GC = 100%)
Docelowy GC: 20-40%
Okno analizy: 20 bp
Wynik: Pozycja 1-20 zawsze będzie miała >50% GC → backtracking niemożliwy
```

### **Problem 2: Brak elastyczności w wyłączaniu sprawdzeń**

**Opis problemu:**
- Wszystkie sprawdzenia walidacji są "na sztywno" w kodzie
- Brak możliwości wyłączenia konkretnych testów (np. hairpin, homodimer)
- Użytkownik nie może dostosować poziomu restrykcyjności
- CLI nie oferuje flag typu `--no-homodimer-check`

## 🎯 **Proponowane Rozwiązania**

### **Rozwiązanie Problemu 1: Adaptacyjne okno walidacji**

#### **Opcja A: Progresywne zaostrzanie kryteriów**
```python
def get_progressive_limits(current_length, target_length, final_min, final_max):
    progress = current_length / target_length
    if progress < 0.3:  # Pierwsze 30% - bardzo szerokie limity
        return 0.1, 0.9
    elif progress < 0.7:  # Środek - postupnie zaostrzamy
        factor = (progress - 0.3) / 0.4
        min_val = 0.1 + factor * (final_min - 0.1)
        max_val = 0.9 - factor * (0.9 - final_max)
        return min_val, max_val
    else:  # Ostatnie 30% - docelowe kryteria
        return final_min, final_max
```

**Implementacja:**
- Dodaj `adaptive_validation: bool = False` do `GeneratorConfig`
- Dodaj `initial_phase_ratio: float = 0.3` (jaki % sekwencji ma luźne kryteria)
- Modyfikuj `_validate_window_with_stats()` do używania progresywnych limitów
- Dodaj parametr CLI `--adaptive-window`

#### **Opcja B: Generator kompatybilnych sekwencji początkowych**
```python
class InitialSequenceGenerator:
    def generate_compatible_initial(self, target_config, min_length=10):
        """Generuje sekwencję początkową kompatybilną z docelowymi parametrami."""
        target_gc = (target_config.min_gc + target_config.max_gc) / 2
        # Oblicz rozkład nukleotydów
        # Generuj losową sekwencję o odpowiednim składzie
        # Waliduj zgodność z wszystkimi kryteriami
```

**Implementacja:**
- Nowa klasa `InitialSequenceGenerator` w `utils.py`
- Metoda `generate_compatible_initial(config, length)`
- Dodaj do CLI `--generate-initial` lub `--auto-initial`

#### **Opcja C: Dwufazowy algorytm**
- **Faza 1:** "Wyrównywanie" - doprowadź sekwencję do zgodności z docelowymi parametrami
- **Faza 2:** "Rozszerzanie" - generuj z pełnymi kryteriami

### **Rozwiązanie Problemu 2: Elastyczne sprawdzenia**

#### **Rozszerzenie GeneratorConfig**
```python
@dataclass
class GeneratorConfig:
    # Obecne parametry...
    
    # Flagi włączania/wyłączania sprawdzeń
    enable_gc_check: bool = True
    enable_tm_check: bool = True  
    enable_hairpin_check: bool = True
    enable_homodimer_check: bool = True
    enable_homopolymer_check: bool = True
    enable_dinucleotide_repeats_check: bool = True
    enable_3prime_stability_check: bool = True
    
    # Lub bardziej elastycznie:
    validation_rules: Dict[str, bool] = field(default_factory=lambda: {
        'gc_content': True,
        'melting_temperature': True,
        'hairpin_structures': True,
        'homodimer_structures': True,
        'homopolymer_runs': True,
        'dinucleotide_repeats': True,
        'three_prime_stability': True
    })
```

#### **Profile walidacji**
```python
VALIDATION_PROFILES = {
    'strict': {     # Wszystkie sprawdzenia (domyślne)
        'gc_content': True,
        'melting_temperature': True,
        'hairpin_structures': True,
        'homodimer_structures': True,
        'homopolymer_runs': True,
        'dinucleotide_repeats': True,
        'three_prime_stability': True
    },
    'relaxed': {    # Tylko podstawowe
        'gc_content': True,
        'homopolymer_runs': True,
        'melting_temperature': True,
        'dinucleotide_repeats': False,
        'three_prime_stability': False,
        'hairpin_structures': False,
        'homodimer_structures': False
    },
    'minimal': {    # Tylko GC i homopolimery
        'gc_content': True,
        'homopolymer_runs': True,
        'melting_temperature': False,
        'dinucleotide_repeats': False,
        'three_prime_stability': False,
        'hairpin_structures': False,
        'homodimer_structures': False
    },
    'sequence_only': {  # Tylko sprawdzenia sekwencyjne
        'homopolymer_runs': True,
        'dinucleotide_repeats': True,
        'three_prime_stability': True,
        'gc_content': False,
        'melting_temperature': False,
        'hairpin_structures': False,
        'homodimer_structures': False
    },
    'none': {       # Brak sprawdzeń (tylko długość)
        'gc_content': False,
        'melting_temperature': False,
        'hairpin_structures': False,
        'homodimer_structures': False,
        'homopolymer_runs': False,
        'dinucleotide_repeats': False,
        'three_prime_stability': False
    }
}
```

#### **Rozszerzenia CLI**
```bash
# Wyłączenie konkretnych sprawdzeń
python -m deterministic_generator \
    --initial ATGC --length 100 \
    --no-hairpin-check \
    --no-homodimer-check \
    --disable-3prime-check

# Profile walidacji
python -m deterministic_generator \
    --initial ATGC --length 100 \
    --profile relaxed

# Granularna kontrola
python -m deterministic_generator \
    --initial ATGC --length 100 \
    --enable-rules gc_content,homopolymers \
    --disable-rules hairpin,homodimer

# Kombininacja z rozwiązaniem problemu 1
python -m deterministic_generator \
    --initial CCCGGGCCCGGG \
    --length 100 \
    --min-gc 0.2 --max-gc 0.4 \
    --adaptive-window \
    --initial-phase-ratio 0.3 \
    --profile relaxed
```

## 📅 **Plan Implementacji**

### **Etap 1: Rozwiązanie Problemu 2 (łatwiejszy, szybki efekt)**

**Zadania:**
- [ ] Dodać flagi sprawdzeń do `GeneratorConfig`
- [ ] Zdefiniować `VALIDATION_PROFILES` w `config.py`
- [ ] Modyfikować `DNAValidator._validate_window_with_stats()` do respektowania flag
- [ ] Rozszerzyć CLI w `__main__.py` o nowe opcje:
  - `--profile {strict,relaxed,minimal,sequence_only,none}`
  - `--no-gc-check`, `--no-tm-check`, `--no-hairpin-check`, etc.
  - `--enable-rules RULES`, `--disable-rules RULES`
- [ ] Dodać testy dla różnych profili
- [ ] Aktualizować dokumentację w README.md

**Oczekiwany rezultat:**
```python
# Przykład użycia
config = GeneratorConfig(
    validation_profile='relaxed',
    enable_hairpin_check=False,
    enable_homodimer_check=False
)
```

### **Etap 2: Rozwiązanie Problemu 1 (trudniejszy, fundamentalny)**

**Zadania - Opcja A (Adaptacyjne okno):**
- [ ] Dodać `adaptive_validation: bool = False` do `GeneratorConfig`
- [ ] Dodać `initial_phase_ratio: float = 0.3` 
- [ ] Dodać `phase_transition_method: str = 'linear'` (linear/exponential/step)
- [ ] Implementować `get_progressive_gc_limits()` i podobne dla innych parametrów
- [ ] Modyfikować `_validate_window_with_stats()` do używania progresywnych limitów
- [ ] Dodać parametr CLI `--adaptive-window`, `--initial-phase-ratio`
- [ ] Dodać logowanie faz walidacji
- [ ] Dodać testy dla konfliktowych sekwencji początkowych

**Zadania - Opcja B (Generator początkowych):**
- [ ] Implementować `InitialSequenceGenerator` w `utils.py`
- [ ] Metoda `generate_compatible_initial(config, min_length, max_attempts)`
- [ ] Algorytm kompozycji nukleotydów dla docelowego GC
- [ ] Walidacja kompatybilności wygenerowanej sekwencji
- [ ] CLI `--generate-initial LENGTH` lub `--auto-initial`
- [ ] Integracja z głównym generatorem

**Oczekiwany rezultat:**
```python
# Adaptacyjne okno
config = GeneratorConfig(
    adaptive_validation=True,
    initial_phase_ratio=0.3,
    min_gc=0.2, max_gc=0.4
)

# Generator początkowych sekwencji  
initial_gen = InitialSequenceGenerator()
compatible_initial = initial_gen.generate_compatible_initial(config, 20)
```

### **Etap 3: Optymalizacje i dodatkowe funkcje**

**Możliwe ulepszenia:**
- [ ] Cache wyników primer3 dla identycznych sekwencji
- [ ] Równoległa walidacja różnych sprawdzeń
- [ ] Inteligentny wybór nukleotydów bazujący na analizie deficytów
- [ ] Tryb "gentle backtracking" z częściową akceptacją suboptimalnych rozwiązań
- [ ] Analiza przyczyn niepowodzeń backtrackingu
- [ ] Eksport profili walidacji do plików YAML/JSON
- [ ] GUI dla konfiguracji parametrów

## 🧪 **Przypadki Testowe do Zaimplementowania**

### **Testy Problemu 1:**
```python
# Test 1: Wysokie GC początkowe → Niskie GC docelowe
initial = "CCCGGGCCCGGGCCCGGG"  # GC = 100%
target_gc = (0.2, 0.4)  # 20-40%

# Test 2: Niskie GC początkowe → Wysokie GC docelowe  
initial = "ATATATATATATATAT"    # GC = 0%
target_gc = (0.6, 0.8)  # 60-80%

# Test 3: Wysokie Tm początkowe → Niskie Tm docelowe
initial = "CGTACGTACGTACGTA"    # Wysokie Tm
target_tm = (40, 50)    # Niskie Tm
```

### **Testy Problemu 2:**
```python
# Test 1: Wyłączenie sprawdzeń hairpin/homodimer
config = GeneratorConfig(
    enable_hairpin_check=False,
    enable_homodimer_check=False
)

# Test 2: Profile walidacji
for profile in ['strict', 'relaxed', 'minimal', 'none']:
    config = GeneratorConfig(validation_profile=profile)
    
# Test 3: Granularna kontrola
config = GeneratorConfig(
    validation_rules={
        'gc_content': True,
        'homopolymer_runs': True,
        'hairpin_structures': False,
        'homodimer_structures': False
    }
)
```

## 🎯 **Rekomendowana Kolejność**

1. **Problem 2** (flagi sprawdzeń) - szybki efekt, duża użyteczność dla użytkowników
2. **Problem 1 - Opcja A** (adaptacyjne okno) - zachowuje architekturę, rozwiązuje większość przypadków
3. **Problem 1 - Opcja B** (generator początkowych) - dla bardziej zaawansowanych przypadków
4. **Optymalizacje** - po przetestowaniu podstawowych rozwiązań

## 📊 **Metryki Sukcesu**

### **Problem 1:**
- [ ] Generator radzi sobie z konfliktowymi parametrami (success rate >80%)
- [ ] Adaptacyjne okno nie pogarsza jakości finalnych sekwencji
- [ ] Czas generowania nie wzrasta >2x dla adaptacyjnego trybu

### **Problem 2:** ✅ **ZREALIZOWANE**
- [x] Wszystkie profile walidacji działają poprawnie
- [x] CLI obsługuje wszystkie nowe opcje
- [x] Dodano przykłady użycia w demo.py
- [x] Backward compatibility z obecnym API zachowana

## 📝 **Notatki Implementacyjne**

- Zachować backward compatibility - domyślne zachowanie bez zmian
- Dodać obszerną dokumentację dla nowych funkcji
- Rozważyć wydajność - adaptacyjne okno może zwiększyć liczbę obliczeń
- Logowanie - dodać szczegółowe logi dla debugowania nowych funkcji
- Testy - pokrycie >90% dla nowych funkcjonalności

---

**Data utworzenia:** 2 września 2025  
**Status:** Analiza zakończona, gotowy do implementacji  
**Priorytet:** Wysoki (Problem 2), Średni (Problem 1)
