import math
from typing import List, Dict
from .base import FactorizationAlgorithm
from .math_utils import is_prime, generate_primes, legendre_symbol


class QuadraticSieveBasic(FactorizationAlgorithm):
    """
    Базовая реализация метода квадратичного решета.
    Поиск гладких чисел — пробным делением (без просеивания).
    """

    def __init__(self):
        super().__init__()

    def _get_factor_base(self, n: int, B: int) -> List[int]:
        primes = generate_primes(B)
        # Включаем только простые p, для которых n является квадратичным вычетом mod p
        # (символ Лежандра = 1), иначе x² ≡ n (mod p) не имеет решений
        factor_base = [2] + [p for p in primes if p > 2 and legendre_symbol(n, p) == 1]

        self.log_step("Этап 1: Построение факторной базы", {
            "message": (
                f"Граница гладкости B = {B}\n"
                f"Включаем простые p ≤ B, для которых символ Лежандра (n/p) = 1,\n"
                f"т.е. уравнение x² ≡ n (mod p) имеет решение.\n"
                f"Это гарантирует, что x² − n делится на такие p.\n"
                f"Итого в базе: {len(factor_base)} простых чисел.\n"
                f"Первые 10: {factor_base[:10]}"
            ),
            "FB": factor_base
        })
        return factor_base

    def _find_smooth_numbers(self, n: int, factor_base: List[int], required_count: int) -> List[Dict]:
        smooth_numbers = []
        x_start = math.isqrt(n) + 1
        x = x_start
        table_data = []
        max_search = 150000
        checked = 0

        self.log_step("Этап 2: Поиск B-гладких чисел (пробное деление)", {
            "message": (
                f"Перебираем x = ⌊√n⌋ + 1, ⌊√n⌋ + 2, ...\n"
                f"Для каждого x вычисляем Q(x) = x² − n.\n"
                f"Q(x) мало по модулю (≈ 2x·Δ), поэтому шансы на гладкость высоки.\n"
                f"Число B-гладкое, если после деления на все p из базы остаток = 1.\n"
                f"Нужно найти: {required_count} гладких чисел (|FB| + 5).\n"
                f"Начинаем с x = {x_start}"
            )
        })

        while len(smooth_numbers) < required_count and (x - x_start) < max_search:
            q_x = x * x - n
            temp_q = q_x
            checked += 1

            # Быстрый отсев: если не делится на малые простые — пропускаем
            if temp_q % 2 != 0 and temp_q % 3 != 0 and temp_q % 5 != 0 and temp_q % 7 != 0:
                x += 1
                continue

            exponents = [0] * len(factor_base)
            for i, p in enumerate(factor_base):
                while temp_q % p == 0:
                    exponents[i] += 1
                    temp_q //= p

            if temp_q == 1:
                smooth_numbers.append({'x': x, 'q_x': q_x, 'exponents': exponents})
                if len(smooth_numbers) <= 12:
                    table_data.append({
                        "x": x,
                        "Q(x) = x²−n": q_x,
                        "Вектор степеней": str(exponents),
                        "Гладкое?": "✓"
                    })

            x += 1

        self.log_step("Результат поиска гладких чисел", {
            "message": (
                f"Просмотрено значений x: {checked}\n"
                f"Найдено B-гладких чисел: {len(smooth_numbers)} / {required_count}\n"
                f"Каждая строка таблицы — это соотношение x² ≡ Q(x) (mod n),\n"
                f"где Q(x) полностью раскладывается по факторной базе."
            ),
            "table": table_data
        })
        return smooth_numbers

    def _gauss_elimination_gf2(self, matrix: List[List[int]], smooth_numbers: List[Dict]) -> List[List[int]]:
        rows = len(matrix)
        cols = len(matrix[0])
        # Расширяем матрицу единичной — для отслеживания комбинаций строк
        M = [matrix[i] + [1 if i == j else 0 for j in range(rows)] for i in range(rows)]

        self.log_step("Этап 3: Линейная алгебра над GF(2)", {
            "message": (
                f"Матрица размером {rows} × {cols} (строки = гладкие числа, столбцы = простые из FB).\n"
                f"Каждый элемент = степень простого mod 2 (чётная/нечётная).\n"
                f"Цель: найти подмножество строк с нулевой суммой mod 2 —\n"
                f"это даст x² ≡ y² (mod n), откуда НОД(x−y, n) — делитель.\n"
                f"Метод: Гаусс над GF(2) с расширенной матрицей."
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
                f"Ранг матрицы: {pivot_row}\n"
                f"Найдено линейных зависимостей: {len(dependencies)}\n"
                f"Каждая зависимость — кандидат на нетривиальный делитель.\n"
                f"Чем больше зависимостей, тем выше шанс успеха."
            )
        })
        return dependencies

    def factorize(self, n: int) -> List[int]:
        self.clear_logs()

        if n <= 1: return [n]
        if n % 2 == 0: return [2, n // 2]
        if is_prime(n): return [n]

        self.log_step("Запуск: Квадратичное решето (базовый)", {
            "message": (
                f"n = {n} ({n.bit_length()} бит)\n"
                f"Факторизация — нахождение простых p, q таких что p × q = n.\n"
                f"Метод Диксона (основа квадратичного решета): ищем x, y такие что x² ≡ y² (mod n).\n"
                f"Тогда n | (x−y)(x+y), и НОД(x−y, n) с вероятностью 1/2 даёт нетривиальный делитель.\n"
                f"Поиск гладких чисел: пробное деление (без просеивания)."
            )
        })

        B = int(math.exp(0.5 * math.sqrt(math.log(n) * math.log(math.log(n)))))
        B = max(B, 50)
        B = min(B, 5000)

        self.log_step("Выбор параметра B", {
            "message": (
                f"Оптимальная граница по L-нотации:\n"
                f"B = exp(0.5 · √(ln n · ln ln n)) ≈ {B}\n"
                f"Слишком малый B → мало гладких чисел.\n"
                f"Слишком большой B → большая факторная база, медленный Гаусс."
            )
        })

        factor_base = self._get_factor_base(n, B)
        required_smooth = len(factor_base) + 5
        smooth_numbers = self._find_smooth_numbers(n, factor_base, required_smooth)

        matrix_mod2 = [[exp % 2 for exp in sn['exponents']] for sn in smooth_numbers]
        dependencies = self._gauss_elimination_gf2(matrix_mod2, smooth_numbers)

        self.log_step("Этап 4: Проверка зависимостей", {
            "message": (
                f"Для каждой зависимости вычисляем:\n"
                f"  X = ∏ xᵢ (mod n)  — левая часть\n"
                f"  Y = ∏ pⱼ^(eⱼ/2) (mod n)  — квадратный корень из ∏ Q(xᵢ)\n"
                f"Проверяем: НОД(X−Y, n) и НОД(X+Y, n)"
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
            for i, p in enumerate(factor_base):
                Y = (Y * pow(p, exponents_sum[i] // 2, n)) % n

            d = math.gcd(abs(X - Y), n)

            self.log_step(f"Зависимость #{idx + 1}", {
                "message": (
                    f"Использованы x: {used_xs[:8]}{'...' if len(used_xs) > 8 else ''}\n"
                    f"X = {X}, Y = {Y}\n"
                    f"НОД(|X−Y|, n) = НОД({abs(X-Y)}, {n}) = {d}\n"
                    f"{'✓ Нетривиальный делитель!' if 1 < d < n else '✗ Тривиальный (1 или n), пробуем следующую зависимость.'}"
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
                "Причина: недостаточно гладких чисел или неудачные комбинации.\n"
                "Решение: увеличить B для расширения факторной базы."
            )
        })
        return [n]
