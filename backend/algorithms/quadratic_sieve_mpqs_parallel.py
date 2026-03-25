import math
import os
import array
import concurrent.futures
from typing import List, Dict
from .base import FactorizationAlgorithm
from .math_utils import is_prime, generate_primes, legendre_symbol, tonelli_shanks

def next_prime(n):
    n = n + 1 if n % 2 == 0 else n + 2
    while not is_prime(n): n += 2
    return n

def _worker_sieve_polynomials(n: int, fb_compact: List[tuple], M_half: int, t_values: List[int], tolerance: float) -> List[Dict]:
    """
    Рабочий процесс: просеивает пачку полиномов и возвращает гладкие числа.
    Запускается в отдельном процессе через ProcessPoolExecutor.
    fb_compact — список кортежей (p, log_p, sqrt_n) вместо словарей, чтобы
    уменьшить объём pickle-сериализации при передаче между процессами.
    """
    smooth_numbers = []
    interval_size = 2 * M_half

    for t in t_values:
        b0 = tonelli_shanks(n, t)
        if b0 is None: continue

        inv_2b0 = pow(2 * b0, -1, t)
        k = (inv_2b0 * ((n - b0**2) // t)) % t
        b = b0 + t * k
        a = t**2
        c = (b**2 - n) // a

        # array.array('f') занимает в ~8 раз меньше памяти, чем list[float]
        sieve_array = array.array('f', bytes(4 * interval_size))

        for p, log_p, root in fb_compact:
            if a % p == 0: continue

            inv_a = pow(a, -1, p)
            idx1 = ((-b + root) * inv_a) % p
            idx2 = ((-b - root) * inv_a) % p
            idx1 = (idx1 - (-M_half)) % p
            idx2 = (idx2 - (-M_half)) % p

            for i in range(idx1, interval_size, p): sieve_array[i] += log_p
            if p != 2 and idx1 != idx2:
                for i in range(idx2, interval_size, p): sieve_array[i] += log_p

        for i in range(interval_size):
            x = i - M_half
            q_x = a * x**2 + 2 * b * x + c
            if q_x <= 0: continue

            if sieve_array[i] >= math.log2(q_x) - tolerance:
                temp_q = q_x
                exps = [0] * len(fb_compact)
                for j, (p, _, _) in enumerate(fb_compact):
                    while temp_q % p == 0:
                        exps[j] += 1
                        temp_q //= p
                if temp_q == 1:
                    real_x = (a * x + b) % n
                    smooth_numbers.append({'x': real_x, 'exponents': exps, 'extra_y': t})

    return smooth_numbers


class QuadraticSieveMPQSParallel(FactorizationAlgorithm):
    """
    Параллельная версия MPQS. Распределяет просеивание полиномов
    на все доступные ядра процессора (Master-Worker паттерн).
    """
    def __init__(self):
        super().__init__()

    def factorize(self, n: int) -> List[int]:
        self.clear_logs()
        if n <= 1: return [n]
        if n % 2 == 0: return [2, n // 2]
        if is_prime(n): return [n]

        num_cores = os.cpu_count() or 1

        ln_n = math.log(n)
        L = math.exp(math.sqrt(ln_n * math.log(ln_n)))
        B = max(100, int(L * 0.35))
        M_half = 10000

        self.log_step("Запуск MPQS (Параллельный)", {
            "message": (
                f"n = {n} ({n.bit_length()} бит)\n"
                f"Факторизация — нахождение простых p, q таких что p × q = n.\n"
                f"Архитектура: Master-Worker с ProcessPoolExecutor.\n"
                f"  • Master (главный процесс): генерирует t-значения, собирает результаты, запускает Гаусс.\n"
                f"  • Workers ({num_cores} процессов): каждый просеивает свою пачку полиномов независимо.\n"
                f"Параллелизм возможен, т.к. полиномы независимы друг от друга.\n"
                f"Параметры: B = {B}, M_half = {M_half}, ядер CPU: {num_cores}"
            )
        })

        # 1. Факторная база (в главном процессе)
        primes = generate_primes(B)
        fb = []
        for p in primes:
            if p == 2:
                fb.append({'p': 2, 'sqrt_n': n % 2, 'log_p': math.log2(2)})
            elif legendre_symbol(n, p) == 1:
                r = tonelli_shanks(n, p)
                if r is not None:
                    fb.append({'p': p, 'sqrt_n': r, 'log_p': math.log2(p)})

        self.log_step("Этап 1: Факторная база (главный процесс)", {
            "message": (
                f"Построена в главном процессе: {len(fb)} простых до B={B}.\n"
                f"Передаётся воркерам как список кортежей (p, log_p, sqrt_n) — меньше pickle-overhead.\n"
                f"Первые 10: {[f['p'] for f in fb[:10]]}"
            )
        })

        # Компактное представление для передачи между процессами (меньше pickle-overhead)
        fb_compact = [(f['p'], f['log_p'], f['sqrt_n']) for f in fb]

        req_count = len(fb) + 5
        tolerance = math.log2(fb[-1]['p'])

        target_t = int((2 * n / M_half)**0.25)
        if target_t % 2 == 0: target_t += 1
        current_t = target_t

        smooth_numbers = []
        poly_generated = 0
        batch_num = 0

        # 2. Параллельное просеивание
        # batch_size = num_cores: один полином на ядро за раз,
        # чтобы не накапливать лишние futures и результаты в памяти
        with concurrent.futures.ProcessPoolExecutor(max_workers=num_cores) as executor:
            while len(smooth_numbers) < req_count:
                batch_size = num_cores
                t_batch = []
                while len(t_batch) < batch_size:
                    if legendre_symbol(n, current_t) == 1:
                        t_batch.append(current_t)
                    current_t = next_prime(current_t)

                poly_generated += len(t_batch)
                batch_num += 1

                # Каждый воркер получает ровно 1 полином
                chunks = [[t] for t in t_batch]

                self.log_step(f"Батч #{batch_num}: раздача задач воркерам", {
                    "message": (
                        f"Сгенерировано {len(t_batch)} значений t (полиномов).\n"
                        f"Разбито на {len(chunks)} чанков по 1 полиному.\n"
                        f"Каждый чанк отправлен в отдельный процесс.\n"
                        f"t от {t_batch[0]} до {t_batch[-1]}"
                    )
                })

                futures = [
                    executor.submit(_worker_sieve_polynomials, n, fb_compact, M_half, chunk, tolerance)
                    for chunk in chunks
                ]

                batch_results = 0
                for future in concurrent.futures.as_completed(futures):
                    result = future.result()
                    smooth_numbers.extend(result)
                    batch_results += len(result)

                self.log_step(f"Батч #{batch_num}: результаты собраны", {
                    "message": (
                        f"Воркеры вернули {batch_results} новых гладких чисел.\n"
                        f"Итого собрано: {len(smooth_numbers)} / {req_count}\n"
                        f"Всего обработано полиномов: {poly_generated}"
                    )
                })

        self.log_step("Этап 2: Параллельное просеивание завершено", {
            "message": (
                f"Итого батчей: {batch_num}\n"
                f"Итого полиномов: {poly_generated}\n"
                f"Итого гладких чисел: {len(smooth_numbers)}\n"
                f"Среднее гладких на полином: {len(smooth_numbers)/max(poly_generated,1):.2f}"
            )
        })

        # 3. Гаусс (в главном процессе)
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

        self.log_step("Этап 3: Метод Гаусса (главный процесс)", {
            "message": (
                f"Матрица {rows} × {cols}: строки — гладкие числа, столбцы — простые из FB.\n"
                f"Элемент [i][j] = степень j-го простого в Q(xᵢ) mod 2 (0 = чётная, 1 = нечётная).\n"
                f"Приводим к ступенчатому виду: XOR строк вместо обычного сложения.\n"
                f"Нулевые строки = линейные зависимости = подмножества, чьё произведение — полный квадрат.\n"
                f"Ранг: {pivot_row}. Зависимостей: {len(dependencies)}.\n"
                f"Гаусс выполняется в главном процессе — все данные уже собраны от воркеров."
            )
        })

        # 4. Проверка
        self.log_step("Этап 4: Проверка зависимостей", {
            "message": (
                f"Для MPQS: Y = (∏ tᵢ) · ∏ pⱼ^(eⱼ/2) mod n\n"
                f"Откуда берётся tᵢ: мы использовали полином Q(x) = ax² + 2bx + c, где a = t².\n"
                f"Тогда (ax + b)² = t²·Q(x) (mod n), поэтому X = ax + b, а t входит в Y под корнем.\n"
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
