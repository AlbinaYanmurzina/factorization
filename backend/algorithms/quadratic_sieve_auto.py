"""
backend/algorithms/quadratic_sieve_auto.py

Квадратичное решето с автоматическим подбором параметров (раздел 6.5 учебного пособия).

Реализует оценку сложности метода квадратичного решета и автоматический
подбор оптимальных параметров B и M на основе L-нотации.

Раздел 6.5: "Оценка сложности метода квадратичного решета"

L-нотация для субэкспоненциальных алгоритмов:
L_n[α, c] = exp(c · (ln n)^α · (ln ln n)^(1-α))

Для квадратичного решета оптимальная сложность:
L_n[1/2, 1] = exp(√(ln n · ln ln n))

Оптимальный размер факторной базы:
B = L_n[1/2, 1/√2] = exp(√(ln n · ln ln n) / √2)

Проблема базовой и оптимизированной версий:
Используют фиксированные формулы для B, что не всегда оптимально:
- Для малых n (< 40 бит) можно использовать меньший B
- Для больших n (> 70 бит) нужен больший B и M

Решение AUTO версии:
Автоматически подбирает B и M в зависимости от размера n,
используя эмпирические коэффициенты для Python реализации.

Эмпирические коэффициенты:
Теоретический оптимум предполагает реализацию на C/C++.
Python медленнее, поэтому используются меньшие коэффициенты α:

Разрядность n | Коэффициент α | Интервал M
--------------|---------------|------------
< 35 бит      | 0.15          | 20000
35-50 бит     | 0.25          | 60000
50-70 бит     | 0.35          | 150000
> 70 бит      | 0.45          | 400000

Формула: B = L · α, где L = exp(√(ln n · ln ln n))

Преимущества:
- Адаптивность к размеру числа
- Избегает избыточных вычислений для малых n
- Масштабируется для больших n
- Прозрачность выбора параметров

Применимость:
Универсальная версия, подходящая для чисел любого размера
в диапазоне применимости квадратичного решета (до 100 бит).
"""

import math
from typing import List, Dict
from .quadratic_sieve_optimized import QuadraticSieveOptimized
from .math_utils import is_prime

class QuadraticSieveAuto(QuadraticSieveOptimized):
    """
    Квадратичное решето с автоматическим подбором параметров.
    
    Наследуется от QuadraticSieveOptimized и переопределяет логику
    выбора параметров B и M на основе оценки сложности (раздел 6.5).
    
    Attributes:
        steps_log: Лог шагов работы алгоритма (наследуется от базового класса)
    """

    def __init__(self):
        """Инициализация AUTO версии квадратичного решета."""
        super().__init__()

    def _calculate_params(self, n: int) -> tuple[int, int]:
        """
        Автоматический расчёт параметров B и M на основе L-нотации.
        
        Реализует оценку сложности из раздела 6.5 учебного пособия.
        
        Алгоритм:
        1. Вычисляем ln n и ln ln n
        2. Вычисляем L = exp(√(ln n · ln ln n))
        3. Определяем α и M по таблице в зависимости от разрядности
        4. Вычисляем B = L · α
        5. Ограничиваем: 100 ≤ B ≤ 15000, 10000 ≤ M ≤ 1000000
        
        Args:
            n (int): Число для факторизации
        
        Returns:
            tuple[int, int]: (B, M) - граница гладкости и размер интервала
        
        Примечание:
            Эмпирические коэффициенты α подобраны для Python реализации.
            Для C/C++ можно использовать большие значения α.
        """
        # Вычисляем логарифмы
        ln_n = math.log(n)
        lnln_n = math.log(ln_n)
        # Вычисляем L-функцию
        L = math.exp(math.sqrt(ln_n * lnln_n))
        # Определяем разрядность числа
        bit_len = n.bit_length()

        # Эмпирические коэффициенты: Python медленнее C++,
        # поэтому B растёт медленнее теоретического оптимума
        if bit_len < 35:
            alpha, M = 0.15, 20000
        elif bit_len < 50:
            alpha, M = 0.25, 60000
        elif bit_len < 70:
            alpha, M = 0.35, 150000
        else:
            alpha, M = 0.45, 400000

        # Вычисляем B по формуле B = L · α
        B = int(L * alpha)
        # Ограничиваем B разумными пределами
        B = max(100, min(B, 15000))
        # Ограничиваем M разумными пределами
        M = max(10000, min(M, 1000000))

        self.log_step("Расчёт параметров по L-нотации", {
            "message": (
                f"n = {n} ({bit_len} бит)\n"
                f"ln n = {ln_n:.3f},  ln ln n = {lnln_n:.3f}\n"
                f"L(n, 1/2) = exp(√(ln n · ln ln n)) = {L:.1f}\n"
                f"Коэффициент α = {alpha}  (зависит от разрядности)\n"
                f"B = L · α = {B}  — граница факторной базы\n"
                f"M = {M}  — размер интервала просеивания\n"
                f"Смысл: при B ≈ L^(1/2) достигается теоретический оптимум QS.\n"
                f"Для Python используем меньший α, чтобы не превысить разумное время."
            )
        })
        return B, M

    def factorize(self, n: int) -> List[int]:
        self.clear_logs()

        if n <= 1: return [n]
        if n % 2 == 0: return [2, n // 2]
        if is_prime(n): return [n]

        # ГЛАВНОЕ ОТЛИЧИЕ: Автоматический расчет
        B, M = self._calculate_params(n)

        self.log_step("Запуск: Квадратичное решето (AUTO)", {
            "message": (
                f"n = {n} ({n.bit_length()} бит)\n"
                f"Отличие от Optimized: параметры B и M подбираются автоматически\n"
                f"по формуле L-нотации, а не фиксированы.\n"
                f"Итоговые параметры: B = {B}, M = {M}"
            )
        })

        try:
            # Используем методы из базового класса Optimized, но с новыми B и M
            factor_base = self._get_factor_base(n, B)
            required_smooth = len(factor_base) + 5
            smooth_numbers = self._sieve_and_find_smooth_custom_m(n, factor_base, required_smooth, M)

            matrix_mod2 = [[exp % 2 for exp in sn['exponents']] for sn in smooth_numbers]
            rows, cols = len(matrix_mod2), len(matrix_mod2[0])

            display_matrix = [row[:30] for row in matrix_mod2[:30]]
            self.log_step("Этап 3: Метод Гаусса над GF(2)", {
                "message": (
                    f"Матрица {rows} × {cols}: строки — гладкие числа, столбцы — простые из FB.\n"
                    f"Элемент [i][j] = степень j-го простого в Q(xᵢ) mod 2 (0 = чётная, 1 = нечётная).\n"
                    f"Приводим к ступенчатому виду: XOR строк вместо обычного сложения.\n"
                    f"Нулевые строки = линейные зависимости = подмножества, чьё произведение — полный квадрат."
                ),
                "matrix_data": display_matrix
            })

            dependencies = self._gauss_elimination_gf2(matrix_mod2)

            self.log_step("Этап 4: Проверка зависимостей", {
                "message": (
                    f"Найдено {len(dependencies)} зависимостей.\n"
                    f"Для каждой: X = ∏xᵢ mod n,  Y = ∏pⱼ^(eⱼ/2) mod n\n"
                    f"Ищем НОД(X−Y, n) ∈ (1, n)."
                )
            })

            for idx, dep in enumerate(dependencies):
                X = 1
                exponents_sum = [0] * len(factor_base)
                for i, is_used in enumerate(dep):
                    if is_used:
                        X = (X * smooth_numbers[i]['x']) % n
                        for j, exp in enumerate(smooth_numbers[i]['exponents']):
                            exponents_sum[j] += exp
                Y = 1
                for i, fb in enumerate(factor_base):
                    Y = (Y * pow(fb['p'], exponents_sum[i] // 2, n)) % n

                d = math.gcd(abs(X - Y), n)
                self.log_step(f"Зависимость #{idx + 1}", {
                    "message": (
                        f"Использовано {sum(dep)} гладких чисел.\n"
                        f"X = {X}, Y = {Y}\n"
                        f"НОД(|X−Y|, n) = {d}\n"
                        f"{'✓ Нетривиальный делитель!' if 1 < d < n else '✗ Тривиальный, продолжаем.'}"
                    )
                })
                if 1 < d < n:
                    self.log_step("Факторизация завершена", {
                        "message": f"{n} = {d} × {n // d}"
                    })
                    return sorted([d, n // d])

            self.log_step("Провал", {
                "message": (
                    "Все зависимости дали тривиальные делители.\n"
                    "Попробуйте увеличить коэффициенты α в расчёте параметров."
                )
            })

        except Exception as e:
            self.log_step("Ошибка", {"message": str(e)})


        return [n]

    def _sieve_and_find_smooth_custom_m(self, n: int, factor_base: List[Dict], required_count: int, M: int) -> List[Dict]:
        """
        Просеивание с настраиваемым размером интервала M.
        
        Копия метода из QuadraticSieveOptimized, но принимает M как параметр
        вместо использования фиксированного значения.
        
        Args:
            n (int): Число для факторизации
            factor_base (List[Dict]): Факторная база с корнями
            required_count (int): Требуемое количество гладких чисел
            M (int): Размер интервала просеивания (настраиваемый)
        
        Returns:
            List[Dict]: Список гладких чисел
        
        Raises:
            ValueError: Если найдено недостаточно гладких чисел
        
        Примечание:
            Этот метод позволяет AUTO версии использовать оптимальный M,
            вычисленный на основе размера n, вместо фиксированного значения.
        """
        x_start = math.isqrt(n) + 1
        sieve_array = [0.0] * M

        for fb in factor_base:
            p, log_p = fb['p'], fb['log_p']
            idx1 = (fb['r1'] - x_start) % p
            idx2 = (fb['r2'] - x_start) % p
            for i in range(idx1, M, p):
                sieve_array[i] += log_p
            if p != 2 and idx1 != idx2:
                for i in range(idx2, M, p):
                    sieve_array[i] += log_p

        smooth_numbers = []
        tolerance = math.log2(factor_base[-1]['p']) if factor_base else 5.0

        for i in range(M):
            x = x_start + i
            q_x = x * x - n
            if q_x <= 0: continue
            if sieve_array[i] >= math.log2(q_x) - tolerance:
                temp_q = q_x
                exponents = [0] * len(factor_base)
                for j, fb in enumerate(factor_base):
                    p = fb['p']
                    while temp_q % p == 0:
                        exponents[j] += 1
                        temp_q //= p
                if temp_q == 1:
                    smooth_numbers.append({'x': x, 'q_x': q_x, 'exponents': exponents})
                if len(smooth_numbers) >= required_count:
                    break
        
        if len(smooth_numbers) < required_count:
            raise ValueError(f"Не хватило интервала M={M}. Найдено {len(smooth_numbers)} из {required_count}")
        return smooth_numbers