import math
from typing import List, Dict
from .base import FactorizationAlgorithm
from .math_utils import is_prime, generate_primes, legendre_symbol, tonelli_shanks


class QuadraticSieveOptimized(FactorizationAlgorithm):
    """
    Оптимизированная реализация метода квадратичного решета.
    Использует логарифмическое просеивание (sieving) и алгоритм Тонелли-Шенкса.
    """

    def __init__(self):
        super().__init__()

    def _get_factor_base(self, n: int, B: int) -> List[Dict]:
        primes = generate_primes(B)
        factor_base = []

        for p in primes:
            if p == 2:
                factor_base.append({'p': 2, 'r1': n % 2, 'r2': n % 2, 'log_p': math.log2(2)})
            elif legendre_symbol(n, p) == 1:
                r = tonelli_shanks(n, p)
                if r is not None:
                    factor_base.append({'p': p, 'r1': r, 'r2': p - r, 'log_p': math.log2(p)})

        self.log_step("Этап 1: Факторная база с корнями (Тонелли–Шенкс)", {
            "message": (
                f"Граница B = {B}. Найдено {len(factor_base)} простых.\n"
                f"Для каждого p вычисляем r: r² ≡ n (mod p) — алгоритм Тонелли–Шенкса.\n"
                f"Корни r₁ = r и r₂ = p − r задают стартовые позиции просеивания:\n"
                f"  x ≡ r₁ (mod p) или x ≡ r₂ (mod p)  ⟹  p | Q(x) = x² − n\n"
                f"Это позволяет просеивать без пробного деления — O(M/p) вместо O(M·|FB|)."
            ),
            "FB": [fb['p'] for fb in factor_base[:20]]
        })
        return factor_base

    def _sieve_and_find_smooth(self, n: int, factor_base: List[Dict], required_count: int) -> List[Dict]:
        M = 50000
        x_start = math.isqrt(n) + 1
        sieve_array = [0.0] * M

        self.log_step("Этап 2: Логарифмическое просеивание", {
            "message": (
                f"Интервал просеивания: x ∈ [{x_start}, {x_start + M})\n"
                f"Для каждого p из FB добавляем log₂(p) в позиции i, где p | Q(x_start + i).\n"
                f"После просеивания: sieve[i] ≈ log₂(Q(xᵢ)) означает, что Q(xᵢ) — B-гладкое.\n"
                f"Допуск: ±log₂(p_max) для компенсации ошибок округления."
            )
        })

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
        table_data = []
        candidates_checked = 0
        tolerance = math.log2(factor_base[-1]['p']) if factor_base else 5.0

        for i in range(M):
            x = x_start + i
            q_x = x * x - n
            if q_x <= 0:
                continue
            target_log = math.log2(q_x)

            if sieve_array[i] >= target_log - tolerance:
                candidates_checked += 1
                temp_q = q_x
                exponents = [0] * len(factor_base)

                for j, fb in enumerate(factor_base):
                    p = fb['p']
                    while temp_q % p == 0:
                        exponents[j] += 1
                        temp_q //= p

                if temp_q == 1:
                    smooth_numbers.append({'x': x, 'q_x': q_x, 'exponents': exponents})
                    if len(smooth_numbers) <= 12:
                        table_data.append({
                            "x": x,
                            "Q(x) = x²−n": q_x,
                            "sieve[i]": round(sieve_array[i], 2),
                            "log₂(Q(x))": round(target_log, 2),
                            "Гладкое?": "✓"
                        })

                if len(smooth_numbers) >= required_count:
                    break

        self.log_step("Результат просеивания", {
            "message": (
                f"Просеян интервал M = {M}\n"
                f"Кандидатов прошло порог: {candidates_checked}\n"
                f"Из них реально B-гладких: {len(smooth_numbers)}\n"
                f"Эффективность: {len(smooth_numbers)}/{candidates_checked} = "
                f"{100*len(smooth_numbers)//max(candidates_checked,1)}% кандидатов оказались гладкими."
            ),
            "table": table_data
        })

        if len(smooth_numbers) < required_count:
            raise ValueError(
                f"Найдено только {len(smooth_numbers)} из {required_count} нужных гладких чисел. "
                f"Увеличьте интервал M или границу B."
            )

        return smooth_numbers

    def _gauss_elimination_gf2(self, matrix: List[List[int]]) -> List[List[int]]:
        if not matrix:
            return []

        rows = len(matrix)
        cols = len(matrix[0])
        M = [matrix[i] + [1 if i == j else 0 for j in range(rows)] for i in range(rows)]

        self.log_step("Этап 3: Метод Гаусса над GF(2)", {
            "message": (
                f"Матрица {rows} × {cols}: строки — гладкие числа, столбцы — простые из FB.\n"
                f"Элемент [i][j] = степень j-го простого в Q(xᵢ) mod 2.\n"
                f"Приводим к ступенчатому виду: XOR строк вместо сложения.\n"
                f"Нулевые строки в результате = линейные зависимости = кандидаты на делитель."
            )
        })

        pivot_row = 0
        for c in range(cols):
            pivot = next((r for r in range(pivot_row, rows) if M[r][c] == 1), -1)
            if pivot == -1:
                continue
            M[pivot_row], M[pivot] = M[pivot], M[pivot_row]
            for r in range(rows):
                if r != pivot_row and M[r][c] == 1:
                    M[r] = [M[r][i] ^ M[pivot_row][i] for i in range(len(M[0]))]
            pivot_row += 1

        dependencies = [M[r][cols:] for r in range(pivot_row, rows)]

        self.log_step("Результат Гаусса", {
            "message": (
                f"Ранг матрицы: {pivot_row} из {cols} столбцов.\n"
                f"Линейных зависимостей: {len(dependencies)}\n"
                f"Каждая зависимость задаёт подмножество гладких чисел,\n"
                f"произведение которых является полным квадратом."
            )
        })
        return dependencies

    def factorize(self, n: int) -> List[int]:
        self.clear_logs()

        if n <= 1: return [n]
        if n % 2 == 0: return [2, n // 2]
        if is_prime(n): return [n]

        self.log_step("Запуск: Квадратичное решето (оптимизированный)", {
            "message": (
                f"n = {n} ({n.bit_length()} бит)\n"
                f"Факторизация — нахождение простых p, q таких что p × q = n.\n"
                f"Метод Диксона: ищем x, y такие что x² ≡ y² (mod n).\n"
                f"Тогда НОД(x−y, n) с вероятностью 1/2 даёт нетривиальный делитель.\n"
                f"Улучшения по сравнению с базовым:\n"
                f"  • Тонелли–Шенкс: вычисляем корни x² ≡ n (mod p) заранее\n"
                f"  • Логарифмическое просеивание: O(M log log M) вместо O(M·|FB|)\n"
                f"  • Допуск по логарифму: отсеиваем заведомо негладкие числа"
            )
        })

        L = math.exp(0.5 * math.sqrt(math.log(n) * math.log(math.log(n))))
        B = int(L ** 0.8)
        B = max(B, 100)
        B = min(B, 10000)

        self.log_step("Выбор параметра B", {
            "message": (
                f"L(n, 1/2) = exp(√(ln n · ln ln n)) ≈ {int(L)}\n"
                f"B = L^0.8 ≈ {B}\n"
                f"Коэффициент 0.8 — эмпирически подобранный баланс\n"
                f"между размером базы и количеством нужных гладких чисел."
            )
        })

        try:
            factor_base = self._get_factor_base(n, B)
            required_smooth = len(factor_base) + 5
            smooth_numbers = self._sieve_and_find_smooth(n, factor_base, required_smooth)

            matrix_mod2 = [[exp % 2 for exp in sn['exponents']] for sn in smooth_numbers]
            dependencies = self._gauss_elimination_gf2(matrix_mod2)

            self.log_step("Этап 4: Проверка зависимостей", {
                "message": (
                    f"Проверяем {len(dependencies)} зависимостей.\n"
                    f"Для каждой: X = ∏xᵢ mod n, Y = ∏pⱼ^(eⱼ/2) mod n\n"
                    f"Ищем НОД(X−Y, n) ∈ (1, n)."
                )
            })

            for idx, dep in enumerate(dependencies):
                X = 1
                exponents_sum = [0] * len(factor_base)
                used_xs = []

                for i, is_used in enumerate(dep):
                    if is_used:
                        X = (X * smooth_numbers[i]['x']) % n
                        used_xs.append(smooth_numbers[i]['x'])
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
                "message": "Все зависимости тривиальны. Попробуйте увеличить B."
            })

        except Exception as e:
            self.log_step("Ошибка выполнения", {"message": str(e)})

        return [n]
