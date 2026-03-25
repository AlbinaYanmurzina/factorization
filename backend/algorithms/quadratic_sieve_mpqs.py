import math
from typing import List, Dict
from .base import FactorizationAlgorithm
from .math_utils import is_prime, generate_primes, legendre_symbol, tonelli_shanks

def next_prime(n):
    n = n + 1 if n % 2 == 0 else n + 2
    while not is_prime(n): n += 2
    return n

class QuadraticSieveMPQS(FactorizationAlgorithm):
    """
    Метод квадратичного решета с множеством полиномов (MPQS).
    Раздел 6.9 учебника. Полиномы Q(x) = ax² + 2bx + c, где a = t².
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
        B = max(100, int(L * 0.3))
        M_half = 5000

        self.log_step("Запуск MPQS (Multiple Polynomial QS)", {
            "message": (
                f"n = {n} ({n.bit_length()} бит)\n"
                f"Факторизация — нахождение простых p, q таких что p × q = n.\n"
                f"Идея MPQS: вместо одного полинома Q(x) = x² − n используем семейство:\n"
                f"  Q(x) = ax² + 2bx + c,  где a = t², b² ≡ n (mod a)\n"
                f"Тогда a·Q(x) = (ax + b)² − n, и значения Q(x) малы на интервале [-M, M].\n"
                f"Малые значения → больше гладких чисел → быстрее набираем базу.\n"
                f"Параметры: B = {B}, интервал [-{M_half}, {M_half}]"
            )
        })

        # 1. Факторная база
        primes = generate_primes(B)
        fb = []
        for p in primes:
            if p == 2:
                fb.append({'p': 2, 'sqrt_n': n % 2, 'log_p': math.log2(2)})
            elif legendre_symbol(n, p) == 1:
                r = tonelli_shanks(n, p)
                if r is not None:
                    fb.append({'p': p, 'sqrt_n': r, 'log_p': math.log2(p)})

        self.log_step("Этап 1: Факторная база", {
            "message": (
                f"Построена факторная база: {len(fb)} простых до B={B}.\n"
                f"Первые 10: {[f['p'] for f in fb[:10]]}"
            )
        })

        smooth_numbers = []
        req_count = len(fb) + 5

        target_t = int((2 * n / M_half)**0.25)
        if target_t % 2 == 0: target_t += 1
        t = target_t

        self.log_step("Выбор начального t", {
            "message": (
                f"Параметр a = t² должен быть ≈ √(2n/M).\n"
                f"Целевое t ≈ (2n/M)^(1/4) = {target_t}\n"
                f"Это обеспечивает |Q(x)| ≈ M·√n/2 — оптимальный размер для гладкости."
            )
        })

        poly_count = 0
        poly_table = []

        # 2. Цикл смены полиномов
        while len(smooth_numbers) < req_count:
            while legendre_symbol(n, t) != 1:
                t = next_prime(t)

            b0 = tonelli_shanks(n, t)

            # Поднятие корня по лемме Гензеля: b² ≡ n (mod t²)
            inv_2b0 = pow(2 * b0, -1, t)
            k = (inv_2b0 * ((n - b0**2) // t)) % t
            b = b0 + t * k
            a = t**2
            c = (b**2 - n) // a

            poly_count += 1
            if len(poly_table) < 5:
                poly_table.append({
                    "t": t, "a = t²": a,
                    "b (Гензель)": b,
                    "c = (b²−n)/a": c,
                    "Q(0) = c": c
                })

            # Просеивание полинома на [-M_half, M_half]
            interval_size = 2 * M_half
            sieve_array = [0.0] * interval_size

            for p_data in fb:
                p, log_p, root = p_data['p'], p_data['log_p'], p_data['sqrt_n']
                if a % p == 0: continue

                inv_a = pow(a, -1, p)
                idx1 = ((-b + root) * inv_a) % p
                idx2 = ((-b - root) * inv_a) % p
                idx1 = (idx1 - (-M_half)) % p
                idx2 = (idx2 - (-M_half)) % p

                for i in range(idx1, interval_size, p): sieve_array[i] += log_p
                if p != 2 and idx1 != idx2:
                    for i in range(idx2, interval_size, p): sieve_array[i] += log_p

            tolerance = math.log2(fb[-1]['p'])

            for i in range(interval_size):
                x = i - M_half
                q_x = a * x**2 + 2 * b * x + c
                if q_x <= 0: continue

                if sieve_array[i] >= math.log2(q_x) - tolerance:
                    temp_q = q_x
                    exps = [0] * len(fb)
                    for j, p_data in enumerate(fb):
                        p = p_data['p']
                        while temp_q % p == 0:
                            exps[j] += 1
                            temp_q //= p
                    if temp_q == 1:
                        # Ключевой трюк: (ax + b)² ≡ a·Q(x) ≡ t²·Q(x) (mod n)
                        # Поэтому X = ax + b, а в Y добавляем t (из a = t²)
                        real_x = (a * x + b) % n
                        smooth_numbers.append({'x': real_x, 'exponents': exps, 'extra_y': t})

            t = next_prime(t)

        self.log_step("Этап 2: Просеивание MPQS завершено", {
            "message": (
                f"Сгенерировано полиномов: {poly_count}\n"
                f"Найдено гладких чисел: {len(smooth_numbers)}\n"
                f"Каждый полином даёт независимые соотношения — нет повторений.\n"
                f"Трюк с X: вместо x берём (ax + b) mod n, т.к. (ax+b)² = a·Q(x) + n·(...)"
            ),
            "Примеры полиномов": poly_table
        })

        # 3. Гаусс
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
                f"Приводим к ступенчатому виду: XOR строк вместо обычного сложения.\n"
                f"Нулевые строки = линейные зависимости = подмножества, чьё произведение — полный квадрат.\n"
                f"Ранг: {pivot_row}. Зависимостей: {len(dependencies)}."
            )
        })

        # 4. Проверка
        self.log_step("Этап 4: Проверка зависимостей", {
            "message": (
                f"Для MPQS: Y = (∏ tᵢ) · ∏ pⱼ^(eⱼ/2) mod n\n"
                f"Откуда берётся tᵢ: мы использовали полином Q(x) = ax² + 2bx + c, где a = t².\n"
                f"Тогда (ax + b)² = a·Q(x) + (b² − n) = t²·Q(x) (mod n).\n"
                f"Значит X = ax + b, а в Y нужно учесть t² — он уходит под корень как t.\n"
                f"Проверяем {len(dependencies)} зависимостей: ищем НОД(X−Y, n) ∈ (1, n)."
            )
        })

        for idx, dep in enumerate(dependencies):
            X, Y = 1, 1
            exponents_sum = [0] * len(fb)
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

        self.log_step("Провал", {"message": "Все зависимости тривиальны."})
        return [n]
