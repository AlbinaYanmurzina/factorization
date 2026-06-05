# -*- coding: utf-8 -*-
"""
tests/test_factorization_correctness.py

Тесты корректности алгоритмов факторизации.
Проверяют, что все алгоритмы возвращают правильное разложение на простые множители.
"""

import pytest
import sys
import os
import time
from sympy import nextprime, isprime
import random

# Добавляем путь к backend для импорта
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'backend'))

from algorithms.pollard import PollardRho
from algorithms.pollard_p1 import PollardP1
from algorithms.quadratic_sieve_basic import QuadraticSieveBasic


def generate_semiprime(bits: int) -> tuple[int, int, int]:
    """
    Генерирует полупростое число (произведение двух простых).
    
    Возвращает: (n, p, q) где n = p * q
    """
    start_p = random.getrandbits(bits // 2)
    start_q = random.getrandbits(bits // 2)
    p = nextprime(max(start_p, 3))
    q = nextprime(max(start_q, 3))
    while q == p:
        q = nextprime(q)
    return p * q, p, q


def verify_factorization(n: int, factors: list[int]) -> bool:
    """
    Проверяет корректность факторизации:
    1. Произведение множителей равно n
    2. Все множители простые
    """
    # Проверка произведения
    product = 1
    for f in factors:
        product *= f
    
    if product != n:
        return False
    
    # Проверка простоты всех множителей
    for f in factors:
        if f > 1 and not isprime(f):
            return False
    
    return True


class TestFactorizationCorrectness:
    """Тесты корректности факторизации"""
    
    @pytest.mark.parametrize("bits", [20, 30, 40, 50])
    def test_pollard_rho_correctness(self, bits):
        """Тест корректности ρ-метода Полларда"""
        n, p, q = generate_semiprime(bits)
        
        algo = PollardRho()
        factors = algo.factorize(n)
        
        assert verify_factorization(n, factors), \
            f"Неверная факторизация для {n}: {factors}"
        
        # Проверяем, что нашли оба простых множителя
        assert set(factors) == {p, q}, \
            f"Ожидали {{{p}, {q}}}, получили {set(factors)}"
    
    @pytest.mark.parametrize("bits", [20, 30, 40, 50])
    def test_pollard_p1_correctness(self, bits):
        """Тест корректности (p-1)-метода Полларда"""
        n, p, q = generate_semiprime(bits)
        
        algo = PollardP1()
        factors = algo.factorize(n)
        
        # (p-1)-метод может не справиться, если (p-1) и (q-1) не гладкие
        # В этом случае он вернёт [n] - это нормально
        if factors == [n]:
            # Алгоритм не справился - это допустимо для (p-1)-метода
            print(f"\n  (p-1)-метод не справился с {n} ({bits} бит) - (p-1) или (q-1) не гладкие")
            return
        
        # Если алгоритм вернул факторизацию, она должна быть корректной
        assert verify_factorization(n, factors), \
            f"Неверная факторизация для {n}: {factors}"
        
        # Проверяем, что нашли правильные множители
        assert set(factors) == {p, q}, \
            f"Ожидали {{{p}, {q}}}, получили {set(factors)}"
    
    @pytest.mark.parametrize("bits", [20, 30, 40])
    def test_quadratic_sieve_correctness(self, bits):
        """Тест корректности алгоритма Диксона"""
        n, p, q = generate_semiprime(bits)
        
        algo = QuadraticSieveBasic()
        factors = algo.factorize(n)
        
        assert verify_factorization(n, factors), \
            f"Неверная факторизация для {n}: {factors}"
        
        # Квадратичное решето должно найти оба множителя
        if len(factors) == 2 and factors != [n]:
            assert set(factors) == {p, q}, \
                f"Ожидали {{{p}, {q}}}, получили {set(factors)}"


class TestFactorizationPerformance:
    """Тесты производительности - проверяем, что время растёт с разрядностью"""
    
    def test_pollard_rho_time_growth(self):
        """Проверяем, что время ρ-метода растёт с разрядностью"""
        times = []
        bit_sizes = [20, 28, 36, 44]
        
        for bits in bit_sizes:
            n, _, _ = generate_semiprime(bits)
            algo = PollardRho()
            
            start = time.perf_counter()
            factors = algo.factorize(n)
            elapsed = (time.perf_counter() - start) * 1000
            
            times.append(elapsed)
            
            # Проверяем корректность
            assert verify_factorization(n, factors), \
                f"Неверная факторизация для {bits} бит: {n} -> {factors}"
        
        print(f"\n ρ-метод Полларда:")
        for bits, t in zip(bit_sizes, times):
            print(f"  {bits} бит: {t:.2f} мс")
        
        # Проверяем общую тенденцию роста (последнее время > первого)
        assert times[-1] > times[0], \
            f"Время должно расти с разрядностью: {times}"
    
    def test_pollard_p1_time_growth(self):
        """Проверяем, что время (p-1)-метода растёт с разрядностью"""
        times = []
        bit_sizes = [20, 28, 36, 44]
        successful_tests = []
        
        for bits in bit_sizes:
            n, _, _ = generate_semiprime(bits)
            algo = PollardP1()
            
            start = time.perf_counter()
            factors = algo.factorize(n)
            elapsed = (time.perf_counter() - start) * 1000
            
            times.append(elapsed)
            
            # (p-1)-метод может не справиться
            if factors != [n]:
                # Проверяем корректность только если алгоритм справился
                assert verify_factorization(n, factors), \
                    f"Неверная факторизация для {bits} бит: {n} -> {factors}"
                successful_tests.append((bits, elapsed))
            else:
                print(f"\n  (p-1)-метод не справился с {bits} бит: {n}")
        
        print(f"\n (p-1)-метод Полларда:")
        for bits, t in zip(bit_sizes, times):
            print(f"  {bits} бит: {t:.2f} мс")
        
        # Проверяем тенденцию роста только для успешных тестов
        if len(successful_tests) >= 2:
            success_times = [t for _, t in successful_tests]
            assert success_times[-1] >= success_times[0] * 0.5, \
                f"Время должно расти с разрядностью (с учетом вариативности): {successful_tests}"
    
    def test_quadratic_sieve_time_growth(self):
        """Проверяем, что время алгоритма Диксона растёт с разрядностью"""
        times = []
        bit_sizes = [20, 28, 36]
        
        for bits in bit_sizes:
            n, _, _ = generate_semiprime(bits)
            algo = QuadraticSieveBasic()
            
            start = time.perf_counter()
            factors = algo.factorize(n)
            elapsed = (time.perf_counter() - start) * 1000
            
            times.append(elapsed)
            
            # Проверяем корректность
            assert verify_factorization(n, factors), \
                f"Неверная факторизация для {bits} бит: {n} -> {factors}"
        
        print(f"\n Алгоритм Диксона:")
        for bits, t in zip(bit_sizes, times):
            print(f"  {bits} бит: {t:.2f} мс")
        
        # Проверяем общую тенденцию роста
        assert times[-1] > times[0], \
            f"Время должно расти с разрядностью: {times}"


class TestComparativePerformance:
    """Сравнительные тесты производительности"""
    
    def test_comparative_benchmark(self):
        """
        Комплексный бенчмарк всех алгоритмов.
        Проверяем логичность результатов.
        """
        bit_sizes = [20, 28, 36, 44, 52]
        results = {
            "pollard_rho": {},
            "pollard_p1": {},
            "qs_basic": {},
        }
        
        print("\n" + "="*70)
        print("СРАВНИТЕЛЬНЫЙ БЕНЧМАРК АЛГОРИТМОВ ФАКТОРИЗАЦИИ")
        print("="*70)
        
        for bits in bit_sizes:
            print(f"\n{'-'*70}")
            print(f"Разрядность: {bits} бит")
            print(f"{'-'*70}")
            
            # Генерируем одно и то же число для всех алгоритмов
            n, p, q = generate_semiprime(bits)
            print(f"Тестовое число: {n}")
            print(f"Истинная факторизация: {p} × {q}")
            print()
            
            # ρ-метод Полларда
            algo_rho = PollardRho()
            start = time.perf_counter()
            factors_rho = algo_rho.factorize(n)
            time_rho = (time.perf_counter() - start) * 1000
            results["pollard_rho"][bits] = time_rho
            
            assert verify_factorization(n, factors_rho), \
                f"ρ-метод: неверная факторизация {n} -> {factors_rho}"
            
            print(f"✓ ρ-метод Полларда:      {time_rho:8.2f} мс  →  {factors_rho}")
            
            # (p-1)-метод Полларда
            algo_p1 = PollardP1()
            start = time.perf_counter()
            factors_p1 = algo_p1.factorize(n)
            time_p1 = (time.perf_counter() - start) * 1000
            results["pollard_p1"][bits] = time_p1
            
            # (p-1)-метод может не справиться
            if factors_p1 == [n]:
                print(f"⊗ (p-1)-метод Полларда:  {time_p1:8.2f} мс  →  не справился (негладкое p-1)")
            else:
                assert verify_factorization(n, factors_p1), \
                    f"(p-1)-метод: неверная факторизация {n} -> {factors_p1}"
                print(f"✓ (p-1)-метод Полларда:  {time_p1:8.2f} мс  →  {factors_p1}")
            
            # Алгоритм Диксона (только для малых разрядностей)
            if bits <= 44:
                algo_qs = QuadraticSieveBasic()
                start = time.perf_counter()
                factors_qs = algo_qs.factorize(n)
                time_qs = (time.perf_counter() - start) * 1000
                results["qs_basic"][bits] = time_qs
                
                # Квадратичное решето может не справиться
                if factors_qs == [n]:
                    print(f"⊗ Алгоритм Диксона:      {time_qs:8.2f} мс  →  не справился")
                else:
                    assert verify_factorization(n, factors_qs), \
                        f"Диксон: неверная факторизация {n} -> {factors_qs}"
                    print(f"✓ Алгоритм Диксона:      {time_qs:8.2f} мс  →  {factors_qs}")
            else:
                results["qs_basic"][bits] = None
                print(f"⊗ Алгоритм Диксона:      пропущен (слишком медленный)")
        
        # Итоговая таблица
        print("\n" + "="*70)
        print("ИТОГОВАЯ ТАБЛИЦА ВРЕМЕНИ ВЫПОЛНЕНИЯ (мс)")
        print("="*70)
        print(f"{'Разрядность':<15} {'ρ-Поллард':<15} {'(p-1)-Поллард':<15} {'Диксон':<15}")
        print("-"*70)
        
        for bits in bit_sizes:
            rho_time = results["pollard_rho"].get(bits, 0)
            p1_time = results["pollard_p1"].get(bits, 0)
            qs_time = results["qs_basic"].get(bits)
            
            qs_str = f"{qs_time:.2f}" if qs_time else "—"
            print(f"{bits:<15} {rho_time:<15.2f} {p1_time:<15.2f} {qs_str:<15}")
        
        print("="*70)
        
        # Проверяем монотонность роста времени
        for algo_name, algo_results in results.items():
            times = [t for t in algo_results.values() if t is not None]
            if len(times) >= 2:
                # Проверяем, что последнее время больше первого
                assert times[-1] > times[0], \
                    f"{algo_name}: время должно расти с разрядностью! {times}"
                
                print(f"\n✓ {algo_name}: время корректно растёт с разрядностью")


class TestEdgeCases:
    """Тесты граничных случаев"""
    
    @pytest.mark.parametrize("algo_class", [PollardRho, PollardP1, QuadraticSieveBasic])
    def test_small_primes(self, algo_class):
        """Тест на малых простых числах"""
        test_cases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
        
        for n in test_cases:
            algo = algo_class()
            factors = algo.factorize(n)
            
            assert factors == [n], \
                f"{algo_class.__name__}: простое число {n} должно вернуть [{n}], получено {factors}"
    
    @pytest.mark.parametrize("algo_class", [PollardRho, PollardP1, QuadraticSieveBasic])
    def test_small_composites(self, algo_class):
        """Тест на малых составных числах"""
        test_cases = [
            (4, [2, 2]),
            (6, [2, 3]),
            (15, [3, 5]),
            (21, [3, 7]),
            (35, [5, 7]),
            (77, [7, 11]),
        ]
        
        for n, expected in test_cases:
            algo = algo_class()
            factors = algo.factorize(n)
            
            assert verify_factorization(n, factors), \
                f"{algo_class.__name__}: неверная факторизация {n} -> {factors}"
            
            assert sorted(factors) == sorted(expected), \
                f"{algo_class.__name__}: для {n} ожидали {expected}, получили {factors}"
    
    @pytest.mark.parametrize("algo_class", [PollardRho, QuadraticSieveBasic])
    def test_powers_of_two(self, algo_class):
        """Тест на степенях двойки"""
        test_cases = [
            (4, [2, 2]),
            (8, [2, 2, 2]),
            (16, [2, 2, 2, 2]),
            (32, [2, 2, 2, 2, 2]),
        ]
        
        for n, expected in test_cases:
            algo = algo_class()
            factors = algo.factorize(n)
            
            assert verify_factorization(n, factors), \
                f"{algo_class.__name__}: неверная факторизация {n} -> {factors}"
            
            assert sorted(factors) == sorted(expected), \
                f"{algo_class.__name__}: для {n} ожидали {expected}, получили {factors}"
    
    def test_pollard_p1_powers_of_two(self):
        """Отдельный тест для (p-1)-метода на степенях двойки"""
        # (p-1)-метод имеет специальную обработку степеней двойки
        test_cases = [
            (4, [2, 2]),
            (8, [2, 2, 2]),
            (16, [2, 2, 2, 2]),
        ]
        
        for n, expected in test_cases:
            algo = PollardP1()
            factors = algo.factorize(n)
            
            assert verify_factorization(n, factors), \
                f"PollardP1: неверная факторизация {n} -> {factors}"
            
            assert sorted(factors) == sorted(expected), \
                f"PollardP1: для {n} ожидали {expected}, получили {factors}"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
