# backend/algorithms/cfrac.py
import math
from typing import List, Dict, Tuple, Optional
from .base import FactorizationAlgorithm
from .math_utils import is_prime, generate_primes, legendre_symbol


class CFRAC(FactorizationAlgorithm):
    """
    Факторизация непрерывными дробями (CFRAC, раздел 3.6).

    Первый субэкспоненциальный алгоритм факторизации (Morrison & Brillhart, 1975).

    Идея: разложение √n в непрерывную дробь даёт подходящие дроби A_k/B_k,
    для которых A_k² ≡ Q_k (mod n), где |Q_k| < 2√n — малые числа.
    Если Q_k раскладывается по факторной базе, получаем соотношение.
    Набрав достаточно соотношений, находим подмножество с квадратным
    произведением через Гаусс над GF(2) → НОД(X−Y, n).

    Отличие от QS: числа Q_k берутся из цепной дроби, а не из полинома x²−n.
    """

    def __init__(self):
        super().__init__()

    # ------------------------------------------------------------------
    # Разложение √n в непрерывную дробь
    # ------------------------------------------------------------------

    def _cfrac_iter(self, n: int):
        """
        Генератор итераций разложения √n в непрерывную дробь.

        На каждом шаге k возвращает (a_k, A_k mod n, Q_k, sign_k),
        где sign_k = (-1)^k (знак Q_k с учётом чередования).

        Рекуррентные формулы (стандартный алгоритм SQUFOF/CFRAC):
            m_0 = 0,  d_0 = 1,  a_0 = floor(√n)
            m_{k+1} = d_k * a_k − m_k
            d_{k+1} = (n − m_{k+1}²) / d_k
            a_{k+1} = floor((a_0 + m_{k+1}) / d_{k+1})

        Числитель подходящей дроби:
            A_{-1} = 1, A_0 = a_0
            A_k = a_k * A_{k-1} + A_{k-2}  (mod n)

        Q_k = d_k  (с точностью до знака (-1)^k).
        """
        a0 = math.isqrt(n)
        if a0 * a0 == n:
            return  # n — точный квадрат, обрабатывается отдельно

        m, d, a = 0, 1, a0

        # Числители подходящих дробей mod n
        A_prev2 = 1      # A_{-1}
        A_prev1 = a0 % n  # A_0

        k = 0
        yield (a, A_prev1, d, 1)  # k=0: Q_0 = d_0 = 1 (тривиально, но для полноты)

        while True:
            m = d * a - m
            d = (n - m * m) // d
            if d == 0:
                return
            a = (a0 + m) // d

            A_new = (a * A_prev1 + A_prev2) % n
            A_prev2, A_prev1 = A_prev1, A_new

            k += 1
            sign = 1 if k % 2 == 0 else -1
            yield (a, A_new, d, sign)

    # ------------------------------------------------------------------
    # Проверка на B-гладкость
    # ------------------------------------------------------------------

    def _try_smooth(self, val: int, factor_base: List[int]) -> Optional[List[int]]:
        """
        Пробует разложить |val| по factor_base.
        Возвращает вектор показателей или None если не гладкое.
        """
        v = abs(val)
        if v == 0:
            return None
        exponents = []
        for p in factor_base:
            e = 0
            while v % p == 0:
                e += 1
                v //= p
            exponents.append(e)
        return exponents if v == 1 else None

    # ------------------------------------------------------------------
    # Гаусс над GF(2)
    # ------------------------------------------------------------------

    def _gauss_gf2(self, matrix: List[List[int]]) -> List[List[int]]:
        rows = len(matrix)
        if rows == 0:
            return []
        cols = len(matrix[0])
        # Расширяем единичной матрицей для отслеживания комбинаций
        M = [row[:] + [1 if i == j else 0 for j in range(rows)]
             for i, row in enumerate(matrix)]

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

        return [M[r][cols:] for r in range(pivot_row, rows)]

    # ------------------------------------------------------------------
    # Основной алгоритм
    # ------------------------------------------------------------------

    def factorize(self, n: int) -> List[int]:
        self.clear_logs()

        if n <= 1:
            return [n]
        if n % 2 == 0:
            return sorted([2, n // 2])

        # Проверка на точный квадрат
        sq = math.isqrt(n)
        if sq * sq == n:
            self.log_step("Точный квадрат", {
                "message": f"n = {sq}² → делители {sq} × {sq}"
            })
            return sorted([sq, sq])

        if is_prime(n):
            self.log_step("Число простое", {
                "message": f"n = {n} — простое (тест Миллера–Рабина)."
            })
            return [n]

        # --- Параметры ---
        B = int(math.exp(0.5 * math.sqrt(math.log(n) * math.log(math.log(n)))))
        B = max(B, 50)
        B = min(B, 8000)

        self.log_step("Запуск: CFRAC (факторизация непрерывными дробями)", {
            "message": (
                f"n = {n} ({n.bit_length()} бит)\n"
                f"CFRAC — первый субэкспоненциальный алгоритм факторизации\n"
                f"(Morrison & Brillhart, 1975).\n\n"
                f"Идея: разложение √n в непрерывную дробь даёт подходящие дроби A_k/B_k.\n"
                f"Для них: A_k² ≡ (−1)^k · d_k  (mod n), где d_k < 2√n — малое число.\n"
                f"Если d_k раскладывается по факторной базе → получаем соотношение.\n"
                f"Набрав |FB|+extra соотношений, Гаусс над GF(2) даёт квадратный корень.\n\n"
                f"Отличие от QS: числа Q_k берутся из цепной дроби, а не из полинома x²−n.\n"
                f"Граница гладкости B = {B}"
            )
        })

        # --- Факторная база ---
        # Включаем -1 (для знака) и простые p с символом Лежандра (n/p) = 1
        raw_primes = generate_primes(B)
        factor_base = [-1] + [2] + [p for p in raw_primes if p > 2 and legendre_symbol(n, p) == 1]

        self.log_step("Этап 1: Построение факторной базы", {
            "message": (
                f"Включаем −1 (для учёта знака Q_k) и простые p ≤ B = {B},\n"
                f"для которых символ Лежандра (n/p) = 1.\n"
                f"Итого в базе: {len(factor_base)} элементов.\n"
                f"Первые 12: {factor_base[:12]}"
            ),
            "FB": factor_base
        })

        required = len(factor_base) + 10
        relations: List[Dict] = []   # {'Ak': int, 'Qk': int, 'exponents': List[int]}
        cfrac_table = []
        max_table_rows = 20
        max_iterations = 200_000

        self.log_step("Этап 2: Сбор соотношений из цепной дроби", {
            "message": (
                f"Разворачиваем √{n} в непрерывную дробь:\n"
                f"  a_0 = ⌊√n⌋ = {math.isqrt(n)}\n"
                f"  m_{{k+1}} = d_k·a_k − m_k\n"
                f"  d_{{k+1}} = (n − m_{{k+1}}²) / d_k\n"
                f"  a_{{k+1}} = ⌊(a_0 + m_{{k+1}}) / d_{{k+1}}⌋\n"
                f"  A_k = a_k·A_{{k-1}} + A_{{k-2}}  (mod n)\n"
                f"  Q_k = (−1)^k · d_k\n\n"
                f"Проверяем Q_k на B-гладкость пробным делением.\n"
                f"Нужно собрать: {required} соотношений."
            )
        })

        for k, (ak, Ak, dk, sign) in enumerate(self._cfrac_iter(n)):
            if k > max_iterations:
                break

            Qk = sign * dk  # знаковое значение

            # Вектор показателей: первый элемент — знак (-1 если Qk < 0)
            sign_exp = 1 if Qk < 0 else 0
            exps = self._try_smooth(dk, factor_base[1:])  # проверяем |Qk| = dk

            if exps is not None:
                full_exp = [sign_exp] + exps
                relations.append({'Ak': Ak, 'Qk': Qk, 'dk': dk, 'exponents': full_exp})

                if len(relations) <= max_table_rows:
                    cfrac_table.append({
                        "k": k,
                        "a_k": ak,
                        "A_k mod n": Ak,
                        "d_k (|Q_k|)": dk,
                        "знак": "−" if Qk < 0 else "+",
                        "Вектор степеней": str(full_exp),
                        "Гладкое?": "✓"
                    })

            if len(relations) >= required:
                break

        self.log_step("Результат сбора соотношений", {
            "message": (
                f"Итераций цепной дроби: {k + 1}\n"
                f"Найдено B-гладких соотношений: {len(relations)} / {required}\n"
                f"Каждая строка: A_k² ≡ Q_k (mod n), Q_k — B-гладкое."
            ),
            "table": cfrac_table
        })

        if len(relations) < 2:
            self.log_step("Провал: недостаточно соотношений", {
                "message": (
                    f"Найдено только {len(relations)} соотношений.\n"
                    f"Попробуйте увеличить B или число итераций."
                )
            })
            return [n]

        # --- Гаусс над GF(2) ---
        matrix_mod2 = [[e % 2 for e in rel['exponents']] for rel in relations]

        display_matrix = [row[:30] for row in matrix_mod2[:30]]
        self.log_step("Этап 3: Линейная алгебра над GF(2)", {
            "message": (
                f"Матрица {len(matrix_mod2)} × {len(factor_base)} "
                f"(строки = соотношения, столбцы = элементы FB).\n"
                f"Элемент = показатель степени mod 2 (чётный/нечётный).\n"
                f"Цель: найти подмножество строк с нулевой суммой mod 2 —\n"
                f"тогда произведение Q_k является точным квадратом.\n"
                f"Метод: Гаусс над GF(2) с расширенной матрицей."
            ),
            "matrix_data": display_matrix
        })

        dependencies = self._gauss_gf2(matrix_mod2)

        self.log_step("Результат Гаусса", {
            "message": (
                f"Ранг матрицы: {len(matrix_mod2) - len(dependencies)}\n"
                f"Найдено линейных зависимостей: {len(dependencies)}\n"
                f"Каждая зависимость — кандидат на нетривиальный делитель."
            )
        })

        # --- Проверка зависимостей ---
        self.log_step("Этап 4: Проверка зависимостей", {
            "message": (
                f"Для каждой зависимости:\n"
                f"  X = ∏ A_k  (mod n)  — левая часть\n"
                f"  Y = √(∏ |Q_k|)  (mod n)  — правая часть\n"
                f"  X² ≡ Y²  (mod n)  →  НОД(X−Y, n)"
            )
        })

        for idx, dep in enumerate(dependencies):
            X = 1
            exp_sum = [0] * len(factor_base)
            used_k = []

            for i, used in enumerate(dep):
                if used:
                    X = (X * relations[i]['Ak']) % n
                    used_k.append(i)
                    for j, e in enumerate(relations[i]['exponents']):
                        exp_sum[j] += e

            # Вычисляем Y = ∏ p_j^(e_j/2) mod n
            # exp_sum[0] — суммарный показатель -1 (должен быть чётным)
            if exp_sum[0] % 2 != 0:
                continue  # знак не сошёлся — пропускаем

            Y = 1
            for j, p in enumerate(factor_base[1:], start=1):
                Y = (Y * pow(p, exp_sum[j] // 2, n)) % n

            d = math.gcd(abs(X - Y), n)

            self.log_step(f"Зависимость #{idx + 1}", {
                "message": (
                    f"Использованы соотношения с индексами: {used_k[:10]}"
                    f"{'...' if len(used_k) > 10 else ''}\n"
                    f"X = {X},  Y = {Y}\n"
                    f"НОД(|X−Y|, n) = НОД({abs(X - Y)}, {n}) = {d}\n"
                    f"{'✓ Нетривиальный делитель!' if 1 < d < n else '✗ Тривиальный, пробуем следующую зависимость.'}"
                )
            })

            if 1 < d < n:
                q = n // d
                self.log_step("Факторизация завершена", {
                    "message": f"{n} = {d} × {q}"
                })
                return sorted([d, q])

        self.log_step("Провал: все зависимости тривиальны", {
            "message": (
                "Все зависимости дали НОД = 1 или n.\n"
                "Причина: недостаточно соотношений или неудачные комбинации.\n"
                "Решение: увеличить B для расширения факторной базы."
            )
        })
        return [n]
