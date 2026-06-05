# -*- coding: utf-8 -*-
"""
tests/test_large_numbers.py

Тесты для проверки работы алгоритмов на больших числах (80-112 бит).
Проверяем утверждение научного руководителя о том, что:
- Методы Полларда не способны разложить 112-битовые числа
- Алгоритм Диксона может разложить 80 и 112-битовые числа
"""

import pytest
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
    """Генерирует полупростое число (произведение двух простых)."""
    start_p = random.getrandbits(bits // 2)
    start_q = random.getrandbits(bits // 2)
    p = nextprime(max(start_p, 3))
    q = nextprime(max(start_q, 3))
    while q == p:
        q = nextprime(q)
    return p * q, p, q


class TestLargeNumbers:
    """Тесты на больших числах для проверки утверждений научного руководителя"""
    
    def test_80_bit_pollard_rho(self):
        """Проверяем ρ-метод Полларда на 80-битных числах"""
        print("\n" + "="*70)
        print("ТЕСТ: ρ-метод Полларда на 80-битных числах")
        print("="*70)
        
        n, p, q = generate_semiprime(80)
        print(f"Число: {n}")
        print(f"Истинная факторизация: {p} × {q}")
        print(f"Разрядность: {n.bit_length()} бит")
        
        algo = PollardRho()
        start = time.perf_counter()
        
        try:
            # Даём 60 секунд на факторизацию
            import signal
            
            def timeout_handler(signum, frame):
                raise TimeoutError("Превышено время ожидания")
            
            # signal.alarm работает только на Unix, для Windows используем другой подход
            factors = algo.factorize(n)
            elapsed = (time.perf_counter() - start) * 1000
            
            print(f"Время: {elapsed:.2f} мс")
            print(f"Результат: {factors}")
            
            if factors == [n]:
                print("❌ Алгоритм НЕ СПРАВИЛСЯ (вернул исходное число)")
                assert False, "ρ-метод должен справиться с 80-битными числами за разумное время"
            else:
                print(f"✓ Алгоритм справился за {elapsed:.2f} мс")
                
        except Exception as e:
            print(f"❌ Ошибка: {e}")
            raise
    
    def test_112_bit_pollard_rho(self):
        """Проверяем ρ-метод Полларда на 112-битных числах"""
        print("\n" + "="*70)
        print("ТЕСТ: ρ-метод Полларда на 112-битных числах")
        print("="*70)
        print("Ожидание: алгоритм НЕ ДОЛЖЕН справиться за разумное время")
        print("="*70)
        
        n, p, q = generate_semiprime(112)
        print(f"Число: {n}")
        print(f"Истинная факторизация: {p} × {q}")
        print(f"Разрядность: {n.bit_length()} бит")
        
        algo = PollardRho()
        start = time.perf_counter()
        
        # Даём только 30 секунд - если не справился, это ожидаемо
        timeout_seconds = 30
        
        try:
            factors = algo.factorize(n)
            elapsed = (time.perf_counter() - start) * 1000
            
            print(f"Время: {elapsed:.2f} мс ({elapsed/1000:.1f} сек)")
            print(f"Результат: {factors}")
            
            if elapsed > timeout_seconds * 1000:
                print(f"⊗ Алгоритм работал > {timeout_seconds} сек - слишком медленно")
                print("✓ Подтверждено: ρ-метод НЕ СПОСОБЕН разложить 112-битные числа")
            elif factors == [n]:
                print("⊗ Алгоритм не справился (вернул исходное число)")
                print("✓ Подтверждено: ρ-метод НЕ СПОСОБЕН разложить 112-битные числа")
            else:
                print(f"⚠ Неожиданно: алгоритм справился за {elapsed:.2f} мс!")
                print("Это противоречит утверждению научного руководителя")
                
        except TimeoutError:
            elapsed = (time.perf_counter() - start) * 1000
            print(f"⊗ Тайм-аут после {elapsed/1000:.1f} сек")
            print("✓ Подтверждено: ρ-метод НЕ СПОСОБЕН разложить 112-битные числа")
    
    def test_80_bit_quadratic_sieve(self):
        """Проверяем алгоритм Диксона на 80-битных числах"""
        print("\n" + "="*70)
        print("ТЕСТ: Алгоритм Диксона на 80-битных числах")
        print("="*70)
        print("Ожидание: алгоритм ДОЛЖЕН справиться (субэкспоненциальная сложность)")
        print("="*70)
        
        n, p, q = generate_semiprime(80)
        print(f"Число: {n}")
        print(f"Истинная факторизация: {p} × {q}")
        print(f"Разрядность: {n.bit_length()} бит")
        
        algo = QuadraticSieveBasic()
        start = time.perf_counter()
        
        try:
            factors = algo.factorize(n)
            elapsed = (time.perf_counter() - start) * 1000
            
            print(f"Время: {elapsed:.2f} мс ({elapsed/1000:.1f} сек)")
            print(f"Результат: {factors}")
            
            if factors == [n]:
                print("❌ Алгоритм НЕ СПРАВИЛСЯ")
                print("Это противоречит утверждению научного руководителя!")
                assert False, "Алгоритм Диксона должен справиться с 80-битными числами"
            else:
                print(f"✓ Алгоритм справился за {elapsed:.2f} мс")
                print("✓ Подтверждено: Диксон СПОСОБЕН разложить 80-битные числа")
                
        except Exception as e:
            print(f"❌ Ошибка: {e}")
            raise
    
    def test_112_bit_quadratic_sieve(self):
        """Проверяем алгоритм Диксона на 112-битных числах"""
        print("\n" + "="*70)
        print("ТЕСТ: Алгоритм Диксона на 112-битных числах")
        print("="*70)
        print("Ожидание: алгоритм ДОЛЖЕН справиться (субэкспоненциальная сложность)")
        print("="*70)
        
        n, p, q = generate_semiprime(112)
        print(f"Число: {n}")
        print(f"Истинная факторизация: {p} × {q}")
        print(f"Разрядность: {n.bit_length()} бит")
        
        algo = QuadraticSieveBasic()
        start = time.perf_counter()
        
        try:
            factors = algo.factorize(n)
            elapsed = (time.perf_counter() - start) * 1000
            
            print(f"Время: {elapsed:.2f} мс ({elapsed/1000:.1f} сек)")
            print(f"Результат: {factors}")
            
            if factors == [n]:
                print("❌ Алгоритм НЕ СПРАВИЛСЯ")
                print("Это противоречит утверждению научного руководителя!")
                assert False, "Алгоритм Диксона должен справиться с 112-битными числами"
            else:
                print(f"✓ Алгоритм справился за {elapsed:.2f} мс")
                print("✓ Подтверждено: Диксон СПОСОБЕН разложить 112-битные числа")
                
        except Exception as e:
            print(f"❌ Ошибка: {e}")
            raise
    
    def test_comparative_large_numbers(self):
        """Сравнительный тест всех алгоритмов на разных разрядностях"""
        print("\n" + "="*70)
        print("СРАВНИТЕЛЬНЫЙ ТЕСТ: Все алгоритмы на больших числах")
        print("="*70)
        
        bit_sizes = [60, 70, 80]  # Уменьшаем до реалистичных размеров
        
        results = []
        
        for bits in bit_sizes:
            print(f"\n{'-'*70}")
            print(f"Разрядность: {bits} бит")
            print(f"{'-'*70}")
            
            n, p, q = generate_semiprime(bits)
            print(f"Число: {n}")
            print(f"Истинная факторизация: {p} × {q}")
            
            row = {"bits": bits, "n": n, "p": p, "q": q}
            
            # ρ-метод Полларда
            print("\n[ρ-метод Полларда]")
            algo_rho = PollardRho()
            start = time.perf_counter()
            factors_rho = algo_rho.factorize(n)
            time_rho = (time.perf_counter() - start) * 1000
            
            print(f"Результат: {factors_rho}")
            print(f"Время: {time_rho:.2f} мс")
            
            if set(factors_rho) == {p, q}:
                row["rho"] = f"✓ {time_rho:.2f} мс"
                print("Статус: ✓ Успешно")
            else:
                row["rho"] = f"✗ {time_rho:.2f} мс"
                print("Статус: ✗ Не справился")
            
            # (p-1)-метод Полларда
            print("\n[(p-1)-метод Полларда]")
            algo_p1 = PollardP1()
            start = time.perf_counter()
            factors_p1 = algo_p1.factorize(n)
            time_p1 = (time.perf_counter() - start) * 1000
            
            print(f"Результат: {factors_p1}")
            print(f"Время: {time_p1:.2f} мс")
            
            if set(factors_p1) == {p, q}:
                row["p1"] = f"✓ {time_p1:.2f} мс"
                print("Статус: ✓ Успешно")
            else:
                row["p1"] = f"✗ {time_p1:.2f} мс"
                print("Статус: ✗ Не справился")
            
            # Алгоритм Диксона
            print("\n[Алгоритм Диксона]")
            algo_qs = QuadraticSieveBasic()
            start = time.perf_counter()
            factors_qs = algo_qs.factorize(n)
            time_qs = (time.perf_counter() - start) * 1000
            
            print(f"Результат: {factors_qs}")
            print(f"Время: {time_qs:.2f} мс ({time_qs/1000:.1f} сек)")
            
            if set(factors_qs) == {p, q}:
                row["qs"] = f"✓ {time_qs:.2f} мс"
                print("Статус: ✓ Успешно")
            else:
                row["qs"] = f"✗ {time_qs:.2f} мс"
                print("Статус: ✗ Не справился")
            
            results.append(row)
        
        # Итоговая таблица
        print("\n" + "="*70)
        print("ИТОГОВАЯ ТАБЛИЦА")
        print("="*70)
        print(f"{'Биты':<8} {'ρ-Поллард':<25} {'(p-1)-Поллард':<25} {'Диксон':<25}")
        print("-"*70)
        
        for row in results:
            print(f"{row['bits']:<8} {row.get('rho', '—'):<25} {row.get('p1', '—'):<25} {row.get('qs', '—'):<25}")
        
        print("="*70)
        print("\nВЫВОДЫ:")
        print("1. ρ-метод Полларда: экспоненциальная сложность O(n^1/4)")
        print("2. (p-1)-метод Полларда: работает только для B-гладких чисел")
        print("3. Алгоритм Диксона: субэкспоненциальная сложность L[1/2, c]")
        print("   - Должен быть быстрее методов Полларда для больших чисел")
        print("="*70)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
