import math
from typing import List
from .base import FactorizationAlgorithm
from .math_utils import is_prime, generate_primes

class PollardP1(FactorizationAlgorithm):
    def __init__(self):
        super().__init__()

    def _p1_step(self, n: int, B: int = 2000) -> int:
        primes = generate_primes(B)

        self.log_step("Инициализация (p−1)-метода", {
            "message": (
                f"Число n = {n}, граница гладкости B = {B}\n"
                f"Идея: если p | n и (p−1) является B-гладким числом,\n"
                f"то по МТФ: a^(p−1) ≡ 1 (mod p), значит p | НОД(a^M − 1, n),\n"
                f"где M = НОК всех простых степеней ≤ B.\n"
                f"База a = 2. Простых в базе: {len(primes)} (до {primes[-1] if primes else '—'})"
            )
        })

        a = 2
        table_data = []
        log_interval = max(1, len(primes) // 8)

        for idx, p in enumerate(primes):
            # Максимальная степень p^k ≤ B
            p_pow = p
            while p_pow * p <= B:
                p_pow *= p

            a_prev = a
            a = pow(a, p_pow, n)
            d = math.gcd(a - 1, n)

            if idx % log_interval == 0 or (1 < d < n):
                table_data.append({
                    "p": p,
                    "p^k ≤ B": p_pow,
                    "a = a^(p^k) mod n": a,
                    "НОД(a−1, n)": d
                })

            if 1 < d < n:
                self.log_step("Промежуточные шаги", {
                    "message": (
                        f"На простом p={p} (степень p^k={p_pow}):\n"
                        f"a = {a_prev}^{p_pow} mod {n} = {a}\n"
                        f"НОД(a−1, n) = НОД({a}−1, {n}) = {d}\n"
                        f"Нетривиальный делитель найден!"
                    ),
                    "table": table_data
                })
                self.log_step("Найден делитель", {
                    "message": f"p−1 оказался B-гладким при B={B}. Делитель: {d}"
                })
                return d

        d = math.gcd(a - 1, n)
        if 1 < d < n:
            self.log_step("Успех на финальном шаге", {
                "message": f"НОД(a−1, n) = {d} после обработки всех простых до B={B}.",
                "table": table_data
            })
            return d

        self.log_step("Неудача при B=" + str(B), {
            "message": (
                f"НОД(a−1, n) = {d} — тривиальный результат.\n"
                f"Вероятная причина: (p−1) содержит простой множитель > B.\n"
                f"Решение: увеличить границу B и повторить."
            ),
            "table": table_data
        })
        return n

    def factorize(self, n: int) -> List[int]:
        self.clear_logs()

        if n <= 1: return [n]
        if n % 2 == 0:
            self.log_step("Тривиальный делитель", {"message": f"{n} чётное → делитель 2"})
            return [2, n // 2]
        if is_prime(n):
            self.log_step("Число простое", {"message": f"{n} — простое, факторизация тривиальна."})
            return [n]

        self.log_step("Стратегия факторизации", {
            "message": (
                f"n = {n}\n"
                f"Факторизация — это нахождение простых чисел, произведение которых равно n.\n"
                f"Например: 15 = 3 × 5, 84 = 2 × 2 × 3 × 7.\n"
                f"(p−1)-метод Полларда работает, когда у делителя p число (p−1) раскладывается\n"
                f"только на малые простые — тогда его можно найти через теорему Ферма.\n"
                f"Пробуем с возрастающими границами B: 100 → 1000 → 10000"
            )
        })

        B_levels = [100, 1000, 10000]
        
        for B in B_levels:
            self.log_step(f"Попытка с B = {B}", {
                "message": f"Запускаем (p−1)-шаг с границей B={B}."
            })
            divisor = self._p1_step(n, B)
            if divisor != n:
                factors = sorted([divisor, n // divisor])
                self.log_step("Факторизация завершена", {
                    "message": f"{n} = {factors[0]} × {factors[1]}"
                })
                return factors

        self.log_step("Алгоритм не справился", {
            "message": (
                f"Ни при одном B делитель не найден.\n"
                f"Вероятно, оба простых множителя имеют негладкое (p−1).\n"
                f"Рекомендуется использовать ρ-метод или квадратичное решето."
            )
        })
        return [n]
