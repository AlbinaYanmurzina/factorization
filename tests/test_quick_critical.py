# -*- coding: utf-8 -*-
"""
tests/test_quick_critical.py

Bystryy test na kritichnyh razryadnostyah: 80, 90, 100 bit
Dlya proverki utverzhdeniy nauchnogo rukovoditelya
"""

import sys
import os
import time
from sympy import nextprime
import random

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'backend'))

from algorithms.pollard import PollardRho
from algorithms.pollard_p1 import PollardP1
from algorithms.quadratic_sieve_basic import QuadraticSieveBasic


def generate_semiprime(bits: int) -> tuple[int, int, int]:
    """Generiruet poluprostoe chislo"""
    start_p = random.getrandbits(bits // 2)
    start_q = random.getrandbits(bits // 2)
    p = nextprime(max(start_p, 3))
    q = nextprime(max(start_q, 3))
    while q == p:
        q = nextprime(q)
    return p * q, p, q


def test_quick_critical():
    """Test na kritichnyh razryadnostyah (bystryy)"""
    print("\n" + "="*80)
    print("BYSTRYY TEST KRITICHNYH RAZRYADNOSTEY: 80, 90, 100 BIT")
    print("="*80)
    
    bit_sizes = [80, 90, 100]
    results = []
    
    # Otkryvaem fayl dlya zapisi rezultatov
    with open('test_results_critical.txt', 'w', encoding='utf-8') as f:
        f.write("="*80 + "\n")
        f.write("REZULTATY TESTA KRITICHNYH RAZRYADNOSTEY\n")
        f.write("="*80 + "\n\n")
        
        for bits in bit_sizes:
            print(f"\n{'='*80}")
            print(f"RAZRYADNOST: {bits} BIT")
            print(f"{'='*80}")
            
            f.write(f"\n{'='*80}\n")
            f.write(f"RAZRYADNOST: {bits} BIT\n")
            f.write(f"{'='*80}\n")
            
            n, p, q = generate_semiprime(bits)
            print(f"Chislo: {n}")
            print(f"Istinnaya faktorizaciya: {p} * {q}")
            print()
            
            f.write(f"Chislo: {n}\n")
            f.write(f"Istinnaya faktorizaciya: {p} * {q}\n\n")
            
            row = {"bits": bits, "n": n, "p": p, "q": q}
            
            # Rho-metod
            print("Rho-metod Pollarda...")
            f.write("Rho-metod Pollarda...\n")
            algo_rho = PollardRho()
            start = time.perf_counter()
            factors_rho = algo_rho.factorize(n)
            time_rho = (time.perf_counter() - start) * 1000
            
            success_rho = (set(factors_rho) == {p, q})
            row["rho_time"] = time_rho
            row["rho_success"] = success_rho
            row["rho_result"] = factors_rho
            
            if success_rho:
                msg = f"  [OK] Uspeh za {time_rho:.0f} ms"
                print(msg)
                f.write(msg + "\n")
            else:
                msg = f"  [X] Ne spravilsya za {time_rho:.0f} ms"
                print(msg)
                print(f"      Rezultat: {factors_rho}")
                f.write(msg + "\n")
                f.write(f"      Rezultat: {factors_rho}\n")
            
            # (p-1)-metod
            print("(p-1)-metod Pollarda...")
            f.write("(p-1)-metod Pollarda...\n")
            algo_p1 = PollardP1()
            start = time.perf_counter()
            factors_p1 = algo_p1.factorize(n)
            time_p1 = (time.perf_counter() - start) * 1000
            
            success_p1 = (set(factors_p1) == {p, q})
            row["p1_time"] = time_p1
            row["p1_success"] = success_p1
            row["p1_result"] = factors_p1
            
            if success_p1:
                msg = f"  [OK] Uspeh za {time_p1:.0f} ms"
                print(msg)
                f.write(msg + "\n")
            else:
                msg = f"  [X] Ne spravilsya za {time_p1:.0f} ms"
                print(msg)
                print(f"      Rezultat: {factors_p1}")
                f.write(msg + "\n")
                f.write(f"      Rezultat: {factors_p1}\n")
            
            # Dikson
            print("Algoritm Diksona...")
            f.write("Algoritm Diksona...\n")
            algo_qs = QuadraticSieveBasic()
            start = time.perf_counter()
            factors_qs = algo_qs.factorize(n)
            time_qs = (time.perf_counter() - start) * 1000
            
            success_qs = (set(factors_qs) == {p, q})
            row["qs_time"] = time_qs
            row["qs_success"] = success_qs
            row["qs_result"] = factors_qs
            
            if success_qs:
                msg = f"  [OK] Uspeh za {time_qs/1000:.1f} sek"
                print(msg)
                f.write(msg + "\n")
            else:
                msg = f"  [X] Ne spravilsya za {time_qs/1000:.1f} sek"
                print(msg)
                print(f"      Rezultat: {factors_qs}")
                f.write(msg + "\n")
                f.write(f"      Rezultat: {factors_qs}\n")
            
            results.append(row)
            f.write("\n")
        
        # Itogovaya tablica
        print(f"\n\n{'='*80}")
        print("ITOGOVAYA TABLICA")
        print(f"{'='*80}")
        
        f.write(f"\n{'='*80}\n")
        f.write("ITOGOVAYA TABLICA\n")
        f.write(f"{'='*80}\n")
        
        header = f"{'Bity':<8} {'Rho-Pollard':<30} {'(p-1)-Pollard':<30} {'Dikson':<30}"
        print(header)
        print("-"*80)
        f.write(header + "\n")
        f.write("-"*80 + "\n")
        
        for row in results:
            rho_str = f"{'[OK]' if row['rho_success'] else '[X]'} {row['rho_time']/1000:.2f}s"
            p1_str = f"{'[OK]' if row['p1_success'] else '[X]'} {row['p1_time']/1000:.2f}s"
            qs_str = f"{'[OK]' if row['qs_success'] else '[X]'} {row['qs_time']/1000:.1f}s"
            
            line = f"{row['bits']:<8} {rho_str:<30} {p1_str:<30} {qs_str:<30}"
            print(line)
            f.write(line + "\n")
        
        print("="*80)
        f.write("="*80 + "\n")
        
        # Analiz
        print(f"\n{'='*80}")
        print("ANALIZ REZULTATOV")
        print(f"{'='*80}")
        
        f.write(f"\n{'='*80}\n")
        f.write("ANALIZ REZULTATOV\n")
        f.write(f"{'='*80}\n")
        
        for row in results:
            msg = f"\n{row['bits']} BIT:"
            print(msg)
            f.write(msg + "\n")
            
            msg = f"  Rho:  {'SPRAVILSYA' if row['rho_success'] else 'NE SPRAVILSYA'} - {row['rho_result']}"
            print(msg)
            f.write(msg + "\n")
            
            msg = f"  p-1:  {'SPRAVILSYA' if row['p1_success'] else 'NE SPRAVILSYA'} - {row['p1_result']}"
            print(msg)
            f.write(msg + "\n")
            
            msg = f"  QS:   {'SPRAVILSYA' if row['qs_success'] else 'NE SPRAVILSYA'} - {row['qs_result']}"
            print(msg)
            f.write(msg + "\n")
        
        # Proverka utverzhdeniy
        print(f"\n{'='*80}")
        print("PROVERKA UTVERZHDENIY NAUCHNOGO RUKOVODITELYA")
        print(f"{'='*80}")
        
        f.write(f"\n{'='*80}\n")
        f.write("PROVERKA UTVERZHDENIY NAUCHNOGO RUKOVODITELYA\n")
        f.write(f"{'='*80}\n")
        
        # Proveryaem 100+ bit
        large_results = [r for r in results if r['bits'] >= 100]
        
        rho_failed_large = all(not r['rho_success'] for r in large_results)
        p1_failed_large = all(not r['p1_success'] for r in large_results)
        qs_success_80_plus = any(r['qs_success'] for r in results if r['bits'] >= 80)
        
        msg = f"\n1. Metody Pollarda NE SPOSOBNY razlozhit 100+ bitnye chisla:"
        print(msg)
        f.write(msg + "\n")
        
        msg = f"   Rho-metod: {'PODTVERZHDENO' if rho_failed_large else 'NE PODTVERZHDENO'}"
        print(msg)
        f.write(msg + "\n")
        
        msg = f"   (p-1)-metod: {'PODTVERZHDENO' if p1_failed_large else 'NE PODTVERZHDENO'}"
        print(msg)
        f.write(msg + "\n")
        
        msg = f"\n2. Tolko Dikson sposoben razlozhit 80+ bitnye chisla:"
        print(msg)
        f.write(msg + "\n")
        
        msg = f"   {'PODTVERZHDENO' if qs_success_80_plus else 'NE PODTVERZHDENO'}"
        print(msg)
        f.write(msg + "\n")
        
        print(f"\n{'='*80}")
        f.write(f"\n{'='*80}\n")
        
        if rho_failed_large and p1_failed_large and qs_success_80_plus:
            msg = "VSE UTVERZHDENIYA PODTVERZHDENY!"
            print(msg)
            f.write(msg + "\n")
        else:
            msg = "Trebuetsya dopolnitelnaya proverka"
            print(msg)
            f.write(msg + "\n")
        
        print(f"{'='*80}\n")
        f.write(f"{'='*80}\n")
        
        print("\nRezultaty sohraneny v fayl: test_results_critical.txt")


if __name__ == "__main__":
    test_quick_critical()
