# backend/algorithms/fermat.py
import math
from typing import List
from .base import FactorizationAlgorithm
from .math_utils import is_prime


class FermatFactorization(FactorizationAlgorithm):
    """
    Метод факторизации Ферма (раздел 3.1).

    Идея: любое нечётное n можно представить как n = x² − y² = (x−y)(x+y).
    Ищем x ≥ ⌈√n⌉ такое, что w = x² − n является полным квадратом.
    Тогда p = x − √w, q = x + √w.

    Работает мгновенно, если делители близки к √n (|p − q| мало).
    Катастрофически медленен, если один делитель мал (тогда x ≈ n/2).
    """

    def __init__(self):
        super().__init__()

    def _isqrt_exact(self, w: int):
        """Возвращает √w если w — полный квадрат, иначе None."""
        if w < 0:
            return None
        s = math.isqrt(w)
        if s * s == w:
            return s
        return None

    def _fermat_step(self, n: int):
        """
        Один проход метода Ферма для числа n.
        Возвращает нетривиальный делитель или n при неудаче.
        """
        if n % 2 == 0:
            self.log_step("Тривиальный делитель", {
                "message": f"{n} чётное → делитель = 2"
            })
            return 2

        x = math.isqrt(n)
        if x * x < n:
            x += 1  # x = ⌈√n⌉

        self.log_step("Инициализация метода Ферма", {
            "message": (
                f"Число n = {n}\n"
                f"Идея: ищем x, y такие что n = x² − y² = (x−y)(x+y).\n"
                f"Тогда делители: p = x − y, q = x + y.\n"
                f"Начинаем с x = ⌈√{n}⌉ = {x}.\n"
                f"На каждом шаге проверяем: является ли w = x² − n полным квадратом.\n"
                f"Если да — нашли y = √w и делители."
            )
        })

        max_log_rows = 20
        table_data = []
        iteration = 0

        while x <= (n + 9) // 2:  # x не может превышать (n+1)/2
            w = x * x - n
            y_exact = self._isqrt_exact(w)
            iteration += 1

            is_perfect = y_exact is not None
            row = {
                "Итерация": iteration,
                "x": x,
                "w = x²−n": w,
                "√w": y_exact if is_perfect else f"≈{math.isqrt(w):.0f} (не целое)",
                "Полный квадрат?": "✓" if is_perfect else "✗",
            }

            if iteration <= max_log_rows or is_perfect:
                table_data.append(row)

            if is_perfect:
                y = y_exact
                p = x - y
                q = x + y

                self.log_step(
                    f"Итерации алгоритма (всего: {iteration})",
                    {
                        "message": (
                            f"Перебираем x = {math.isqrt(n) if math.isqrt(n)*math.isqrt(n) >= n else math.isqrt(n)+1}, "
                            f"{math.isqrt(n)+1 if math.isqrt(n)*math.isqrt(n) >= n else math.isqrt(n)+2}, ...\n"
                            f"Для каждого x вычисляем w = x² − n и проверяем, является ли w полным квадратом.\n"
                            f"Показаны первые {min(iteration, max_log_rows)} итераций + финальная."
                        ),
                        "table": table_data,
                    },
                )

                self.log_step("Найдено представление n = x² − y²", {
                    "message": (
                        f"На итерации {iteration}:\n"
                        f"  x = {x},  w = x² − n = {x}² − {n} = {w}\n"
                        f"  √w = {y} (целое число!) ✓\n"
                        f"  n = {x}² − {y}² = ({x}−{y})·({x}+{y})\n"
                        f"  p = x − y = {x} − {y} = {p}\n"
                        f"  q = x + y = {x} + {y} = {q}\n"
                        f"  Проверка: {p} × {q} = {p * q} {'= n ✓' if p * q == n else '≠ n ✗'}"
                    )
                })

                if p == 1:
                    # Нашли тривиальное разложение n = 1 × n
                    self.log_step("Тривиальное разложение", {
                        "message": (
                            f"p = 1, q = {q} — тривиальное разложение.\n"
                            f"Это означает, что n = {n} является простым числом\n"
                            f"(или алгоритм дошёл до x = (n+1)/2 без успеха)."
                        )
                    })
                    return n

                return p  # возвращаем меньший делитель

            x += 1

        # Дошли до конца — n простое
        self.log_step(f"Итерации алгоритма (всего: {iteration})", {
            "message": (
                f"Перебрали все x от ⌈√n⌉ до (n+1)/2.\n"
                f"Полный квадрат не найден → n простое."
            ),
            "table": table_data,
        })
        return n

    def factorize(self, n: int) -> List[int]:
        self.clear_logs()

        if n <= 1:
            return [n]
        if n % 2 == 0:
            return sorted([2, n // 2])
        if is_prime(n):
            self.log_step("Число простое", {
                "message": f"n = {n} — простое (тест Миллера–Рабина). Факторизация не требуется."
            })
            return [n]

        self.log_step("Запуск: Метод Ферма", {
            "message": (
                f"n = {n} ({n.bit_length()} бит)\n"
                f"Метод Ферма — «дедушка» квадратичного решета.\n"
                f"Основан на представлении n = x² − y² = (x−y)(x+y).\n"
                f"Сложность: O(|p−q|/2) итераций, где p, q — делители.\n"
                f"• Если p ≈ q (делители близки к √n) → работает мгновенно.\n"
                f"• Если p ≪ q (один делитель мал) → катастрофически медленно.\n"
                f"Именно этот недостаток привёл к созданию квадратичного решета."
            )
        })

        factors = []
        stack = [n]

        while stack:
            current = stack.pop()

            if current <= 1:
                continue

            if is_prime(current):
                factors.append(current)
                self.log_step("Простое число найдено", {
                    "message": f"{current} — простое (тест Миллера–Рабина). Добавляем в результат."
                })
                continue

            self.log_step("Составное число", {
                "message": f"{current} — составное. Запускаем метод Ферма."
            })

            divisor = self._fermat_step(current)

            if divisor == current:
                # Не удалось разложить (не должно случиться для составных, но на всякий случай)
                factors.append(current)
            else:
                quotient = current // divisor
                self.log_step("Разбиение числа", {
                    "message": (
                        f"{current} = {divisor} × {quotient}\n"
                        f"Оба числа отправляются на дальнейшую проверку."
                    )
                })
                stack.append(divisor)
                stack.append(quotient)

        factors.sort()
        self.log_step("Факторизация завершена", {
            "message": f"Итоговое разложение: {n} = {' × '.join(map(str, factors))}"
        })
        return factors
