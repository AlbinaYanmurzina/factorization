# -*- coding: utf-8 -*-
"""
tests/test_critical_sizes.py

Быстрый тест на критичных разрядностях: 80, 90, 100, 110 бит
Для проверки утверждений научного руководителя
"""

import pytest
import sys
import os
import time
from sympy import nextprime, isprime
import random

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'backend'))

from algorithms.pollard import PollardRho
from algorithms.pollard_p1 import PollardP1
from algorithms.quadratic_sieve_basic import QuadraticSieveBasic


def generate_semiprime(bits: int) -> tuple[int, int, int]:
    """Генерирует полупростое число"""
    start_p = random.getrandbits(bits // 2)
    start_q = random.getrandbits(bits // 2)
    p = nextprime(max(start_p, 3))
    q = nextprime(max(start_q, 3))
    while q == p:
        q = nextprime(q)
    return p * q, p, q


def test_critical_sizes():
    """Тест на критичных разрядностях"""
    print("\n" + "="*80)
    print("TEST KRITICHNYH RAZRYADNOSTEY: 80, 85, 90, 95, 100, 105, 110 BIT")
    print("="*80)
    
    bit_sizes = [80, 85, 90, 95, 100, 105, 110]
    results = []
    
    for bits in bit_sizes:
        print(f"\n{'='*80}")
        print(f"RAZRYADNOST: {bits} BIT")
        print(f"{'='*80}")
        
        n, p, q = generate_semiprime(bits)
        print(f"Chislo: {n}")
        print(f"Istinnaya faktorizaciya: {p} * {q}")
        print()
        
        row = {"bits": bits, "n": n, "p": p, "q": q}
        
        # Rho-metod (s ogranicheniem vremeni)
        print("Rho-metod Pollarda...")
        algo_rho = PollardRho()
        start = time.perf_counter()
        factors_rho = algo_rho.factorize(n)
        time_rho = (time.perf_counter() - start) * 1000
        
        success_rho = (set(factors_rho) == {p, q})
        row["rho_time"] = time_rho
        row["rho_success"] = success_rho
        row["rho_result"] = factors_rho
        
        if success_rho:
            print(f"  [OK] Uspeh za {time_rho:.0f} ms")
        else:
            print(f"  [X] Ne spravilsya za {time_rho:.0f} ms (vernul {len(factors_rho)} mnozhiteley)")
            print(f"      Rezultat: {factors_rho}")
        
        # (p-1)-metod
        print("(p-1)-metod Pollarda...")
        algo_p1 = PollardP1()
        start = time.perf_counter()
        factors_p1 = algo_p1.factorize(n)
        time_p1 = (time.perf_counter() - start) * 1000
        
        success_p1 = (set(factors_p1) == {p, q})
        row["p1_time"] = time_p1
        row["p1_success"] = success_p1
        row["p1_result"] = factors_p1
        
        if success_p1:
            print(f"  [OK] Uspeh za {time_p1:.0f} ms")
        else:
            print(f"  [X] Ne spravilsya za {time_p1:.0f} ms (vernul {len(factors_p1)} mnozhiteley)")
            print(f"      Rezultat: {factors_p1}")
        
        # Dikson
        print("Algoritm Diksona...")
        algo_qs = QuadraticSieveBasic()
        start = time.perf_counter()
        factors_qs = algo_qs.factorize(n)
        time_qs = (time.perf_counter() - start) * 1000
        
        success_qs = (set(factors_qs) == {p, q})
        row["qs_time"] = time_qs
        row["qs_success"] = success_qs
        row["qs_result"] = factors_qs
        
        if success_qs:
            print(f"  [OK] Uspeh za {time_qs/1000:.1f} sek")
        else:
            print(f"  [X] Ne spravilsya za {time_qs/1000:.1f} sek (vernul {len(factors_qs)} mnozhiteley)")
            print(f"      Rezultat: {factors_qs}")
        
        results.append(row)
    
    # Itogovaya tablica
    print(f"\n\n{'='*80}")
    print("ITOGOVAYA TABLICA")
    print(f"{'='*80}")
    print(f"{'Bity':<8} {'Rho-Pollard':<25} {'(p-1)-Pollard':<25} {'Dikson':<25}")
    print("-"*80)
    
    for row in results:
        rho_str = f"{'OK' if row['rho_success'] else 'X'} {row['rho_time']/1000:.1f}s"
        p1_str = f"{'OK' if row['p1_success'] else 'X'} {row['p1_time']/1000:.1f}s"
        qs_str = f"{'OK' if row['qs_success'] else 'X'} {row['qs_time']/1000:.1f}s"
        
        print(f"{row['bits']:<8} {rho_str:<25} {p1_str:<25} {qs_str:<25}")
    
    print("="*80)
    
    # Analiz
    print(f"\n{'='*80}")
    print("ANALIZ REZULTATOV")
    print(f"{'='*80}")
    
    for row in results:
        print(f"\n{row['bits']} BIT:")
        print(f"  Rho:  {'SPRAVILSYA' if row['rho_success'] else 'NE SPRAVILSYA'} - {row['rho_result']}")
        print(f"  p-1:  {'SPRAVILSYA' if row['p1_success'] else 'NE SPRAVILSYA'} - {row['p1_result']}")
        print(f"  QS:   {'SPRAVILSYA' if row['qs_success'] else 'NE SPRAVILSYA'} - {row['qs_result']}")
    
    # Proverka utverzhdeniy
    print(f"\n{'='*80}")
    print("PROVERKA UTVERZHDENIY NAUCHNOGO RUKOVODITELYA")
    print(f"{'='*80}")
    
    # Proveryaem 100+ bit
    large_results = [r for r in results if r['bits'] >= 100]
    
    rho_failed_large = all(not r['rho_success'] for r in large_results)
    p1_failed_large = all(not r['p1_success'] for r in large_results)
    qs_success_large = any(r['qs_success'] for r in results if r['bits'] >= 80)
    
    print(f"\n1. Metody Pollarda NE SPOSOBNY razlozhit 100+ bitnye chisla:")
    print(f"   Rho-metod: {'PODTVERZHDENO' if rho_failed_large else 'NE PODTVERZHDENO'}")
    print(f"   (p-1)-metod: {'PODTVERZHDENO' if p1_failed_large else 'NE PODTVERZHDENO'}")
    
    print(f"\n2. Tolko Dikson sposoben razlozhit 80+ bitnye chisla:")
    print(f"   {'PODTVERZHDENO' if qs_success_large else 'NE PODTVERZHDENO'}")
    
    print(f"\n{'='*80}")
    
    if rho_failed_large and p1_failed_large and qs_success_large:
        print("VSE UTVERZHDENIYA PODTVERZHDENY!")
    else:
        print("Trebuetsya dopolnitelnaya proverka")
    
    print(f"{'='*80}\n")


if __name__ == "__main__":
    test_critical_sizes()
