# backend/algorithms/pollard.py
import math
import random
from typing import List
from .base import FactorizationAlgorithm
from .math_utils import is_prime

class PollardRho(FactorizationAlgorithm):
    
    def __init__(self):
        super().__init__()

    def _rho_step(self, n: int) -> int:
        if n % 2 == 0:
            self.log_step("Тривиальный делитель", {"message": f"{n} чётное → делитель = 2"})
            return 2

        x = 2
        y = 2
        d = 1
        c = 1
        f = lambda val: (pow(val, 2, n) + c) % n

        self.log_step("Инициализация ρ-метода", {
            "message": (
                f"Число n = {n}\n"
                f"Полином f(x) = (x² + {c}) mod {n}\n"
                f"Черепаха x₀ = {x}, Заяц y₀ = {y}\n"
                f"Идея: ищем цикл в последовательности f(x) mod p для неизвестного делителя p.\n"
                f"Когда x ≡ y (mod p), то НОД(|x−y|, n) даст нетривиальный делитель."
            )
        })

        iteration = 0
        max_logs = 15
        table_data = []

        while d == 1:
            x = f(x)
            y = f(f(y))
            d = math.gcd(abs(x - y), n)
            iteration += 1

            if iteration <= max_logs or d > 1:
                table_data.append({
                    "i": iteration,
                    "x (черепаха)": x,
                    "y (заяц)": y,
                    "|x − y|": abs(x - y),
                    "НОД(|x−y|, n)": d
                })

        self.log_step(f"Итерации алгоритма (всего: {iteration})", {
            "message": (
                f"Черепаха делает 1 шаг: x = f(x)\n"
                f"Заяц делает 2 шага: y = f(f(y))\n"
                f"На каждом шаге проверяем НОД(|x−y|, n).\n"
                f"Показаны первые {min(iteration, max_logs)} итераций + финальная."
            ),
            "table": table_data
        })

        if d == n:
            self.log_step("Вырождение алгоритма", {
                "message": (
                    f"НОД = n = {n} — нашли тривиальный делитель.\n"
                    f"Это значит x ≡ y (mod n), а не только mod p — цикл замкнулся слишком рано.\n"
                    f"Причина: неудачный выбор c={c} или начальной точки x₀=2.\n"
                    f"Примечание: в этой учебной реализации перезапуск с другим c не реализован.\n"
                    f"В полноценном алгоритме здесь меняют c (например, c=2, c=3, ...) и повторяют."
                )
            })
            return n

        self.log_step("Найден нетривиальный делитель", {
            "message": (
                f"На итерации {iteration}: x={x}, y={y}\n"
                f"НОД(|{x} − {y}|, {n}) = НОД({abs(x-y)}, {n}) = {d}\n"
                f"Делитель {d} нетривиален: 1 < {d} < {n} ✓"
            )
        })
        return d

    def factorize(self, n: int) -> List[int]:
        self.clear_logs()
        
        if n <= 1:
            return [n]

        factors = []
        stack = [n]

        self.log_step("Начало факторизации", {
            "message": (
                f"Раскладываем n = {n} на простые множители.\n"
                f"Факторизация — это нахождение простых чисел, произведение которых равно n.\n"
                f"Например: 15 = 3 × 5, 84 = 2 × 2 × 3 × 7.\n"
                f"Стратегия ρ-метода: рекурсивно ищем делители через алгоритм Флойда.\n"
                f"Каждое составное число разбивается на два, пока все не станут простыми."
            )
        })

        while stack:
            current = stack.pop()
            
            if is_prime(current):
                factors.append(current)
                self.log_step("Простое число найдено", {
                    "message": f"{current} — простое (тест Миллера–Рабина). Добавляем в результат."
                })
                continue

            self.log_step("Составное число", {
                "message": f"{current} — составное. Запускаем ρ-шаг для поиска делителя."
            })

            divisor = self._rho_step(current)
            
            if divisor == current:
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
