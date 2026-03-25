import math
from typing import List, Dict
from .base import FactorizationAlgorithm
from .math_utils import is_prime, generate_primes, legendre_symbol, tonelli_shanks

class QuadraticSieveLPV(FactorizationAlgorithm):
    """
    Квадратичное решето с Вариацией Большого Множителя (Large Prime Variation).
    Раздел 6.8 учебника.
    """
    def __init__(self):
        super().__init__()

    def factorize(self, n: int) -> List[int]:
        self.clear_logs()
        if n <= 1: return [n]
        if n % 2 == 0: return [2, n // 2]
        if is_prime(n): return [n]

        ln_n = math.log(n)
        L = math.exp(math.sqrt(ln_n * math.log(ln_n)))
        bit_len = n.bit_length()

        if bit_len < 40:
            B, M = int(L * 0.2), 30000
        else:
            B, M = int(L * 0.3), 100000
        B = max(100, min(B, 10000))

        self.log_step("Запуск QS-LPV (Large Prime Variation)", {
            "message": (
                f"n = {n} ({bit_len} бит)\n"
                f"Факторизация — нахождение простых p, q таких что p × q = n.\n"
                f"Идея LPV: принимаем не только B-гладкие числа, но и полугладкие —\n"
                f"числа вида Q(x) = smooth_part × P, где P — один большой простой B < P < B².\n"
                f"Два полугладких с одинаковым P перемножаются → полностью гладкое произведение.\n"
                f"Это резко увеличивает количество используемых соотношений.\n"
                f"Параметры: B = {B}, M = {M}, допустимый P < B² = {B**2}"
            )
        })

        # 1. Факторная база
        primes = generate_primes(B)
        fb = []
        for p in primes:
            if p == 2:
                fb.append({'p': 2, 'r1': n % 2, 'r2': n % 2, 'log_p': math.log2(2)})
            elif legendre_symbol(n, p) == 1:
                r = tonelli_shanks(n, p)
                if r is not None:
                    fb.append({'p': p, 'r1': r, 'r2': p - r, 'log_p': math.log2(p)})

        self.log_step("Этап 1: Факторная база", {
            "message": (
                f"Построена факторная база: {len(fb)} простых чисел до B={B}.\n"
                f"Для каждого p вычислены корни r₁, r₂: r² ≡ n (mod p) (Тонелли–Шенкс).\n"
                f"Первые 10 простых: {[f['p'] for f in fb[:10]]}"
            )
        })

        # 2. Просеивание
        x_start = math.isqrt(n) + 1
        sieve_array = [0.0] * M

        for p_data in fb:
            p, log_p = p_data['p'], p_data['log_p']
            idx1 = (p_data['r1'] - x_start) % p
            idx2 = (p_data['r2'] - x_start) % p
            for i in range(idx1, M, p): sieve_array[i] += log_p
            if p != 2 and idx1 != idx2:
                for i in range(idx2, M, p): sieve_array[i] += log_p

        smooth_numbers = []
        partials = {}
        lpv_matches = 0
        req_count = len(fb) + 5
        # Расширенный допуск: принимаем числа с остатком до B²
        tolerance = math.log2(fb[-1]['p']) + 2.0

        self.log_step("Этап 2: Просеивание с LPV", {
            "message": (
                f"Просеиваем интервал M = {M} значений x.\n"
                f"Расширенный допуск: принимаем кандидатов с sieve[i] ≥ log₂(Q(x)) − {tolerance:.1f}.\n"
                f"После пробного деления на FB:\n"
                f"  • остаток = 1 → полностью гладкое, сразу в базу\n"
                f"  • 1 < остаток < B² → полугладкое, сохраняем в словарь partials[P]\n"
                f"  • два полугладких с одним P → склеиваем в одно соотношение"
            )
        })

        smooth_table = []
        partial_table = []

        for i in range(M):
            x = x_start + i
            q_x = x * x - n
            if q_x <= 0: continue

            if sieve_array[i] >= math.log2(q_x) - (2 * math.log2(B)):
                temp_q = q_x
                exps = [0] * len(fb)
                for j, p_data in enumerate(fb):
                    p = p_data['p']
                    while temp_q % p == 0:
                        exps[j] += 1
                        temp_q //= p

                if temp_q == 1:
                    smooth_numbers.append({'x': x, 'exponents': exps, 'extra_y': 1})
                    if len(smooth_table) < 6:
                        smooth_table.append({"x": x, "Q(x)": q_x, "тип": "полностью гладкое"})
                elif 1 < temp_q < B * B:
                    large_prime = temp_q
                    if large_prime in partials:
                        lpv_matches += 1
                        match = partials.pop(large_prime)
                        comb_x = (x * match['x']) % n
                        comb_exps = [e1 + e2 for e1, e2 in zip(exps, match['exponents'])]
                        smooth_numbers.append({'x': comb_x, 'exponents': comb_exps, 'extra_y': large_prime})
                        if len(partial_table) < 6:
                            partial_table.append({
                                "P": large_prime,
                                "x₁": match['x'], "x₂": x,
                                "склеено": "✓"
                            })
                    else:
                        partials[large_prime] = {'x': x, 'exponents': exps}

            if len(smooth_numbers) >= req_count: break

        self.log_step("Результат просеивания LPV", {
            "message": (
                f"Всего соотношений: {len(smooth_numbers)} / {req_count}\n"
                f"  • Полностью гладких: {len(smooth_numbers) - lpv_matches}\n"
                f"  • Склеено из полугладких пар: {lpv_matches}\n"
                f"  • Неиспользованных полугладких в ожидании: {len(partials)}\n"
                f"LPV дал прирост: +{lpv_matches} соотношений"
            ),
            "Примеры гладких": smooth_table,
            "Примеры склеенных пар": partial_table
        })

        if len(smooth_numbers) < req_count:
            self.log_step("Недостаточно соотношений", {
                "message": f"Найдено {len(smooth_numbers)}, нужно {req_count}. Увеличьте B или M."
            })
            return [n]

        # 3. Гаусс над GF(2)
        matrix_mod2 = [[exp % 2 for exp in sn['exponents']] for sn in smooth_numbers]
        rows, cols = len(matrix_mod2), len(matrix_mod2[0])
        M_mat = [matrix_mod2[i] + [1 if i == j else 0 for j in range(rows)] for i in range(rows)]

        pivot_row = 0
        for c in range(cols):
            pivot = next((r for r in range(pivot_row, rows) if M_mat[r][c] == 1), -1)
            if pivot == -1: continue
            M_mat[pivot_row], M_mat[pivot] = M_mat[pivot], M_mat[pivot_row]
            for r in range(rows):
                if r != pivot_row and M_mat[r][c] == 1:
                    M_mat[r] = [M_mat[r][i] ^ M_mat[pivot_row][i] for i in range(len(M_mat[0]))]
            pivot_row += 1
        dependencies = [M_mat[r][cols:] for r in range(pivot_row, rows)]

        self.log_step("Этап 3: Метод Гаусса над GF(2)", {
            "message": (
                f"Матрица {rows} × {cols}: строки — гладкие числа, столбцы — простые из FB.\n"
                f"Элемент [i][j] = степень j-го простого в Q(xᵢ) mod 2 (0 = чётная, 1 = нечётная).\n"
                f"Приводим к ступенчатому виду методом Гаусса, но сложение заменено на XOR.\n"
                f"Нулевые строки в результате = линейные зависимости:\n"
                f"  подмножество гладких чисел, чьё произведение — полный квадрат.\n"
                f"Ранг: {pivot_row}. Зависимостей: {len(dependencies)}.\n"
                f"Важно для LPV: соотношения со склеенными парами содержат P² под произведением,\n"
                f"поэтому при вычислении Y дополнительно умножаем на P — он уходит под корень."
            )
        })

        # 4. Проверка зависимостей
        self.log_step("Этап 4: Проверка зависимостей", {
            "message": f"Проверяем {len(dependencies)} зависимостей на нетривиальность НОД."
        })

        for idx, dep in enumerate(dependencies):
            X = 1
            exponents_sum = [0] * len(fb)
            Y = 1

            for i, is_used in enumerate(dep):
                if is_used:
                    X = (X * smooth_numbers[i]['x']) % n
                    Y = (Y * smooth_numbers[i]['extra_y']) % n
                    for j, exp in enumerate(smooth_numbers[i]['exponents']):
                        exponents_sum[j] += exp

            for i, p_data in enumerate(fb):
                Y = (Y * pow(p_data['p'], exponents_sum[i] // 2, n)) % n

            d = math.gcd(abs(X - Y), n)

            self.log_step(f"Зависимость #{idx + 1}", {
                "message": (
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

        self.log_step("Провал", {"message": "Все зависимости тривиальны. Увеличьте B или M."})
        return [n]
