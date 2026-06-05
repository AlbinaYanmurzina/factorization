# -*- coding: utf-8 -*-
"""
tests/test_very_large_numbers.py

Тесты для очень больших чисел (80-112 бит) с шагом 5 бит.
Проверяем способность алгоритма Диксона факторизовать большие числа.
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
    """Генерирует полупростое число (произведение двух простых)."""
    start_p = random.getrandbits(bits // 2)
    start_q = random.getrandbits(bits // 2)
    p = nextprime(max(start_p, 3))
    q = nextprime(max(start_q, 3))
    while q == p:
        q = nextprime(q)
    return p * q, p, q


def verify_factorization(n: int, factors: list[int], p: int, q: int) -> bool:
    """Проверяет корректность факторизации"""
    if not factors:
        return False
    
    # Проверка произведения
    product = 1
    for f in factors:
        product *= f
    
    if product != n:
        return False
    
    # Проверка что нашли правильные множители
    return set(factors) == {p, q}


class TestVeryLargeNumbers:
    """Тесты на очень больших числах (80-112 бит)"""
    
    def test_large_numbers_step_5(self):
        """Тест всех алгоритмов на числах 80-112 бит с шагом 5"""
        print("\n" + "="*80)
        print("ТЕСТ: Факторизация больших чисел (80-112 бит, шаг 5)")
        print("="*80)
        print("\nЦель: проверить утверждение научного руководителя:")
        print("  - Методы Полларда НЕ СПОСОБНЫ разложить 112-битные числа")
        print("  - Только алгоритм Диксона способен разложить 80 и 112-битные числа")
        print("="*80)
        
        bit_sizes = [80, 85, 90, 95, 100, 105, 110]  # Полный диапазон
        results = []
        
        for bits in bit_sizes:
            print(f"\n{'='*80}")
            print(f"РАЗРЯДНОСТЬ: {bits} бит")
            print(f"{'='*80}")
            
            n, p, q = generate_semiprime(bits)
            print(f"Число n = {n}")
            print(f"Истинная факторизация: {p} * {q}")
            print(f"Фактическая разрядность: {n.bit_length()} бит")
            print()
            
            row = {
                "bits": bits,
                "n": n,
                "p": p,
                "q": q,
                "actual_bits": n.bit_length()
            }
            
            # Rho-метод Полларда
            print(f"{'-'*80}")
            print("RHO-МЕТОД ПОЛЛАРДА")
            print(f"{'-'*80}")
            
            algo_rho = PollardRho()
            start = time.perf_counter()
            
            # Ограничиваем время выполнения 120 секундами
            import signal
            
            def timeout_handler(signum, frame):
                raise TimeoutError()
            
            try:
                # Для Windows используем другой подход - просто запускаем
                factors_rho = algo_rho.factorize(n)
                time_rho = (time.perf_counter() - start) * 1000
            except TimeoutError:
                time_rho = 120000  # 120 секунд
                factors_rho = [n]  # Не справился
            
            success_rho = verify_factorization(n, factors_rho, p, q)
            
            print(f"Результат: {factors_rho}")
            print(f"Время: {time_rho:.2f} мс ({time_rho/1000:.2f} сек)")
            
            if success_rho:
                row["rho_status"] = "[OK] Успешно"
                row["rho_time"] = f"{time_rho:.2f} мс"
                row["rho_result"] = str(factors_rho)
                print(f"Статус: [OK] УСПЕШНО")
            else:
                row["rho_status"] = "[X] Не справился"
                row["rho_time"] = f"{time_rho:.2f} мс"
                row["rho_result"] = str(factors_rho)
                if time_rho >= 120000:
                    print(f"Статус: [X] ТАЙМ-АУТ (120 сек)")
                elif factors_rho == [n]:
                    print(f"Статус: [X] НЕ СПРАВИЛСЯ (вернул исходное число)")
                else:
                    print(f"Статус: [X] НЕ СПРАВИЛСЯ (вернул {factors_rho})")
            
            # (p-1)-метод Полларда
            print(f"\n{'-'*80}")
            print("(p-1)-МЕТОД ПОЛЛАРДА")
            print(f"{'-'*80}")
            
            algo_p1 = PollardP1()
            start = time.perf_counter()
            factors_p1 = algo_p1.factorize(n)
            time_p1 = (time.perf_counter() - start) * 1000
            
            success_p1 = verify_factorization(n, factors_p1, p, q)
            
            print(f"Результат: {factors_p1}")
            print(f"Время: {time_p1:.2f} мс ({time_p1/1000:.2f} сек)")
            
            if success_p1:
                row["p1_status"] = "[OK] Успешно"
                row["p1_time"] = f"{time_p1:.2f} мс"
                row["p1_result"] = str(factors_p1)
                print(f"Статус: [OK] УСПЕШНО")
            else:
                row["p1_status"] = "[X] Не справился"
                row["p1_time"] = f"{time_p1:.2f} мс"
                row["p1_result"] = str(factors_p1)
                if factors_p1 == [n]:
                    print(f"Статус: [X] НЕ СПРАВИЛСЯ (вернул исходное число - (p-1) не B-гладкое)")
                else:
                    print(f"Статус: [X] НЕ СПРАВИЛСЯ (вернул {factors_p1})")
            
            # Алгоритм Диксона
            print(f"\n{'-'*80}")
            print("АЛГОРИТМ ДИКСОНА (Квадратичное решето)")
            print(f"{'-'*80}")
            
            algo_qs = QuadraticSieveBasic()
            start = time.perf_counter()
            factors_qs = algo_qs.factorize(n)
            time_qs = (time.perf_counter() - start) * 1000
            
            success_qs = verify_factorization(n, factors_qs, p, q)
            
            print(f"Результат: {factors_qs}")
            print(f"Время: {time_qs:.2f} мс ({time_qs/1000:.1f} сек = {time_qs/60000:.1f} мин)")
            
            if success_qs:
                row["qs_status"] = "[OK] Успешно"
                row["qs_time"] = f"{time_qs:.2f} мс"
                row["qs_result"] = str(factors_qs)
                print(f"Статус: [OK] УСПЕШНО")
            else:
                row["qs_status"] = "[X] Не справился"
                row["qs_time"] = f"{time_qs:.2f} мс"
                row["qs_result"] = str(factors_qs)
                if factors_qs == [n]:
                    print(f"Статус: [X] НЕ СПРАВИЛСЯ (вернул исходное число)")
                else:
                    print(f"Статус: [X] НЕ СПРАВИЛСЯ (вернул {factors_qs})")
            
            results.append(row)
            
            # Промежуточная сводка
            print(f"\n{'='*80}")
            print(f"СВОДКА ДЛЯ {bits} БИТ:")
            print(f"  Rho-Поллард:     {row['rho_status']:<20} {row['rho_time']}")
            print(f"  (p-1)-Поллард: {row['p1_status']:<20} {row['p1_time']}")
            print(f"  Диксон:        {row['qs_status']:<20} {row['qs_time']}")
            print(f"{'='*80}")
        
        # Итоговая таблица
        print(f"\n\n{'='*80}")
        print("ИТОГОВАЯ ТАБЛИЦА РЕЗУЛЬТАТОВ")
        print(f"{'='*80}")
        print(f"{'Биты':<6} {'Rho-Поллард':<30} {'(p-1)-Поллард':<30} {'Диксон':<30}")
        print(f"{'-'*80}")
        
        for row in results:
            rho_str = f"{row['rho_status'][:3]} {row['rho_time']}"
            p1_str = f"{row['p1_status'][:3]} {row['p1_time']}"
            qs_str = f"{row['qs_status'][:3]} {row['qs_time']}"
            print(f"{row['bits']:<6} {rho_str:<30} {p1_str:<30} {qs_str:<30}")
        
        print(f"{'='*80}")
        
        # Детальная таблица с результатами
        print(f"\n\n{'='*80}")
        print("ДЕТАЛЬНАЯ ТАБЛИЦА (ЧТО ВЕРНУЛИ АЛГОРИТМЫ)")
        print(f"{'='*80}")
        
        for row in results:
            print(f"\n{row['bits']} БИТ:")
            print(f"  Число n = {row['n']}")
            print(f"  Истинная факторизация: {row['p']} * {row['q']}")
            print(f"  ")
            print(f"  Rho-метод вернул:  {row.get('rho_result', 'N/A')}")
            print(f"  (p-1)-метод вернул: {row.get('p1_result', 'N/A')}")
            print(f"  Диксон вернул:     {row.get('qs_result', 'N/A')}")
        
        print(f"\n{'='*80}")
        
        # Анализ результатов
        print("\n" + "="*80)
        print("АНАЛИЗ РЕЗУЛЬТАТОВ")
        print("="*80)
        
        rho_success_count = sum(1 for r in results if r['rho_status'].startswith('✓'))
        p1_success_count = sum(1 for r in results if r['p1_status'].startswith('✓'))
        qs_success_count = sum(1 for r in results if r['qs_status'].startswith('✓'))
        
        print(f"\nУспешных факторизаций из {len(results)} тестов:")
        print(f"  Rho-метод Полларда:     {rho_success_count}/{len(results)} ({rho_success_count/len(results)*100:.0f}%)")
        print(f"  (p-1)-метод Полларда: {p1_success_count}/{len(results)} ({p1_success_count/len(results)*100:.0f}%)")
        print(f"  Алгоритм Диксона:     {qs_success_count}/{len(results)} ({qs_success_count/len(results)*100:.0f}%)")
        
        print("\n" + "="*80)
        print("ВЫВОДЫ:")
        print("="*80)
        
        # Проверяем утверждения научного руководителя
        rho_failed_large = all(not r['rho_status'].startswith('✓') for r in results if r['bits'] >= 100)
        p1_failed_large = all(not r['p1_status'].startswith('✓') for r in results if r['bits'] >= 100)
        qs_success_large = any(r['qs_status'].startswith('✓') for r in results if r['bits'] >= 80)
        
        print("\n1. Методы Полларда на больших числах (>=100 бит):")
        if rho_failed_large:
            print("   OK ПОДТВЕРЖДЕНО: Rho-метод НЕ СПОСОБЕН разложить 100+ битные числа")
        else:
            print("   X Неожиданно: Rho-метод справился с некоторыми 100+ битными числами")
        
        if p1_failed_large:
            print("   OK ПОДТВЕРЖДЕНО: (p-1)-метод НЕ СПОСОБЕН разложить 100+ битные числа")
        else:
            print("   X Неожиданно: (p-1)-метод справился с некоторыми 100+ битными числами")
        
        print("\n2. Алгоритм Диксона на больших числах (>=80 бит):")
        if qs_success_large:
            print("   OK ПОДТВЕРЖДЕНО: Диксон СПОСОБЕН разложить 80+ битные числа")
        else:
            print("   X Диксон не справился с 80+ битными числами")
        
        print("\n3. Сравнение времени выполнения:")
        for row in results:
            if row['qs_status'].startswith('✓'):
                qs_time_sec = float(row['qs_time'].split()[0]) / 1000
                print(f"   {row['bits']} бит: Диксон = {qs_time_sec:.1f} сек")
        
        print("\n" + "="*80)
        print("ЗАКЛЮЧЕНИЕ:")
        print("="*80)
        print("Утверждения научного руководителя:")
        print("  1. Методы Полларда НЕ СПОСОБНЫ разложить 112-битные числа")
        print("  2. Только Диксон способен разложить 80 и 112-битные числа")
        print()
        
        if rho_failed_large and p1_failed_large and qs_success_large:
            print("*** ВСЕ УТВЕРЖДЕНИЯ ПОДТВЕРЖДЕНЫ ***")
        else:
            print("! Некоторые утверждения требуют дополнительной проверки")
        
        print("="*80)


if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])
