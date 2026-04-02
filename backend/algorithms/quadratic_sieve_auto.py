import math
from typing import List, Dict
from .quadratic_sieve_optimized import QuadraticSieveOptimized
from .math_utils import is_prime

class QuadraticSieveAuto(QuadraticSieveOptimized):
    """
    Версия 'Auto': интеллектуальный подбор параметров B и M 
    на основе теоретической сложности L-нотации.
    """

    def __init__(self):
        super().__init__()

    def _calculate_params(self, n: int) -> tuple[int, int]:
        """
        Теоретический расчет параметров на основе L-нотации:
        L_n[1/2, alpha] = exp(alpha * sqrt(ln n * ln ln n))
        """
        ln_n = math.log(n)
        lnln_n = math.log(ln_n)
        L = math.exp(math.sqrt(ln_n * lnln_n))
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

        B = int(L * alpha)
        B = max(100, min(B, 15000))
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
        Копия метода из Optimized, но с внешним параметром M.
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