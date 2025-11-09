#!/usr/bin/env python3
"""
Benchmark generatora sekwencji DNA dla różnych trybów i ustawień.

Przykłady:
  # Szybki benchmark bez primer3 (profil sekwencyjny)
  python benchmarks/benchmark_generator.py --profile sequence_only \
      --modes deterministic random \
      --heuristics on off \
      --window-sizes 10 20 30 \
      --runs 10 --length 80 --initial CCTGTCATCACGCTAGTAAC

  # Zapis wyników do CSV
  python benchmarks/benchmark_generator.py --profile sequence_only --runs 20 \
      --csv bench.csv

Uwaga: Profile z termodynamiką (np. strict) wymagają primer3-py.
"""

from __future__ import annotations

import argparse
import time
from statistics import mean
from typing import Dict, Any, List

import sys
from pathlib import Path

# Ensure package import path (repo root)
ROOT = Path(__file__).resolve().parents[1]  # package dir
PKG_PARENT = ROOT.parent
if str(PKG_PARENT) not in sys.path:
    sys.path.insert(0, str(PKG_PARENT))

from dna_generator import (
    DNAGenerator,
    GeneratorConfig,
    GenerationMode,
)


def run_single(config: GeneratorConfig, initial: str, length: int) -> Dict[str, Any]:
    gen = DNAGenerator(config)
    t0 = time.perf_counter()
    result = gen.generate(initial, length)
    dt = time.perf_counter() - t0
    stats = result.generation_stats or {}
    return {
        'success': result.success,
        'time_s': dt,
        'attempts': stats.get('total_attempts', 0),
        'backtrack': stats.get('backtrack_count', 0),
        'max_depth': stats.get('max_depth_reached', 0),
    }


def summarize(results: List[Dict[str, Any]]) -> Dict[str, Any]:
    if not results:
        return {}
    success_count = sum(1 for r in results if r['success'])
    times = [r['time_s'] for r in results]
    attempts = [r['attempts'] for r in results]
    backtracks = [r['backtrack'] for r in results]
    return {
        'success_rate': success_count / len(results),
        'time_avg_s': mean(times),
        'time_p95_s': sorted(times)[int(0.95 * (len(times) - 1))],
        'attempts_avg': mean(attempts),
        'backtrack_avg': mean(backtracks),
        'n': len(results),
    }


def build_config(mode: str, profile: str, heuristics: str, window_size: int, seed: int | None) -> GeneratorConfig:
    m = GenerationMode.DETERMINISTIC if mode == 'deterministic' else GenerationMode.RANDOM
    cfg = GeneratorConfig(
        generation_mode=m,
        seed=seed,
        validation_profile=profile,
        window_size=window_size,
        enable_progress_logging=False,
    )
    cfg.enable_backtrack_heuristics = (heuristics == 'on')
    return cfg


def main():
    p = argparse.ArgumentParser(description='Benchmark generatora sekwencji DNA')
    p.add_argument('--initial', default='CCTGTCATCACGCTAGTAAC', help='Sekwencja początkowa')
    p.add_argument('--length', type=int, default=80, help='Docelowa długość sekwencji')
    p.add_argument('--runs', type=int, default=20, help='Liczba uruchomień na kombinację parametrów')
    p.add_argument('--modes', nargs='+', default=['deterministic', 'random'], help='Tryby: deterministic/random')
    p.add_argument('--profile', default='sequence_only', help='Profil walidacji (np. sequence_only/strict/relaxed)')
    p.add_argument('--heuristics', nargs='+', default=['on', 'off'], help='Heurystyki: on/off')
    p.add_argument('--window-sizes', nargs='+', type=int, default=[10, 20, 30], help='Rozmiary okna')
    p.add_argument('--seed', type=int, default=12345, help='Seed dla trybu deterministycznego')
    p.add_argument('--csv', type=str, help='Plik CSV na wyniki')
    p.add_argument('--json', type=str, help='Plik JSON na wyniki')

    args = p.parse_args()

    rows: List[Dict[str, Any]] = []

    for mode in args.modes:
        for heur in args.heuristics:
            for ws in args.window_sizes:
                cfg = build_config(mode, args.profile, heur, ws, args.seed if mode == 'deterministic' else None)
                results = [run_single(cfg, args.initial, args.length) for _ in range(args.runs)]
                summary = summarize(results)
                row = {
                    'mode': mode,
                    'profile': args.profile,
                    'heuristics': heur,
                    'window_size': ws,
                    **summary,
                }
                rows.append(row)

    # Wypisz w formacie tabelarycznym
    header = ['mode', 'profile', 'heuristics', 'window_size', 'success_rate', 'time_avg_s', 'time_p95_s', 'attempts_avg', 'backtrack_avg', 'n']
    print("\t".join(header))
    for r in rows:
        print("\t".join(str(r[h]) for h in header))

    if args.csv:
        import csv
        with open(args.csv, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=header)
            w.writeheader()
            w.writerows(rows)
        print(f"Zapisano CSV: {args.csv}")

    if args.json:
        import json
        with open(args.json, 'w') as f:
            json.dump(rows, f, indent=2)
        print(f"Zapisano JSON: {args.json}")


if __name__ == '__main__':
    main()
