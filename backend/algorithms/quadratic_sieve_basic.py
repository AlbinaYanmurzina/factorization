"""
backend/algorithms/quadratic_sieve_basic.py

Базовая реализация метода квадратичного решета (раздел 6.1 учебного пособия).

Реализует идею Мориса Крейтчика и алгоритм Диксона - основу для всех
вариантов квадратичного решета.

Метод Диксона (1981):
Ищем пары (x, y) такие, что x² ≡ y² (mod n), но x ≢ ±y (mod n).
Тогда НОД(x - y, n) с вероятностью ≥ 1/2 даёт нетривиальный делитель.

Алгоритм:
1. Построение факторной базы FB (раздел 6.3)
2. Поиск B-гладких чисел Q(x) = x² - n через пробное деление
3. Решение системы линейных уравнений методом Гаусса над GF(2) (раздел 6.4)
4. Проверка зависимостей для нахождения делителя

Сложность: L_n[1/2, 1] = exp(√(ln n · ln ln n))
Память: O(M), где M - размер интервала просеивания

Особенности базовой версии:
- Использует пробное деление (медленно, но просто для понимания)
- Не использует просеивание (в отличие от метода Померанца)
- Идеальна для обучения и понимания основ метода

Исторический контекст:
Квадратичное решето было первым алгоритмом, превзошедшим экспоненциальную
сложность для общей факторизации. До появления решета числового поля (NFS)
оставался самым быстрым методом для чисел до 100 десятичных знаков.
"""

import math
import time
from typing import List, Dict
from .base import FactorizationAlgorithm
from .math_utils import is_prime, generate_primes, legendre_symbol


class QuadraticSieveBasic(FactorizationAlgorithm):
    """
    Реализация метода квадратичного решета на основе алгоритма.
    
    Реализует классический алгоритм Диксона с пробным делением для
    поиска B-гладких чисел.
    
    """

    def __init__(self):
        """Инициализация базового квадратичного решета."""
        super().__init__()

    def _get_factor_base(self, n: int, B: int) -> tuple[List[int], float]:
        """
        Построение факторной базы (раздел 6.3 учебного пособия).
        
        Возвращает: (factor_base, time_ms)
        """
        _fb_start = time.perf_counter()
        
        # Генерируем все простые числа до границы B с помощью решета Эратосфена.
        # Это кандидаты для факторной базы — но не все из них подойдут.
        primes = generate_primes(B)

        # Фильтруем: оставляем только простые p, для которых уравнение x² ≡ n (mod p)
        # имеет решение. Это проверяется символом Лежандра: (n/p) = 1.
        # Смысл: если (n/p) ≠ 1, то Q(x) = x² − n никогда не делится на p,
        # и включать p в базу бессмысленно — гладких чисел с таким p не найти.
        # p = 2 включаем всегда (символ Лежандра для p=2 не определён стандартно).
        factor_base = [2] + [p for p in primes if p > 2 and legendre_symbol(n, p) == 1]

        _fb_ms = (time.perf_counter() - _fb_start) * 1000
        
        pi_B_all = len(primes) + 1  # +1 за двойку
        self.log_step("Этап 1: Построение факторной базы", {
            "message": (
                f"Граница гладкости B = {B}\n"
                f"Всего простых p ≤ B: π({B}) = {pi_B_all}\n"
                f"Из них включаем только p, для которых символ Лежандра (n/p) = 1,\n"
                f"т.е. уравнение x² ≡ n (mod p) имеет решение.\n"
                f"Это гарантирует, что x² − n делится на такие p.\n\n"
                f"Размер факторной базы |FB| = {len(factor_base)} элементов\n"
                f"(из π({B}) = {pi_B_all} простых чисел ≤ B)\n"
                f"Первые 10: {factor_base[:10]}\n\n"
                f"Время построения: {_fb_ms:.3f} мс"
            ),
            "FB": factor_base,
            "FB_size": len(factor_base),
            "pi_B": pi_B_all,
            "B": B,
            "fb_time_ms": round(_fb_ms, 3),
        })
        return factor_base, _fb_ms

    def _find_smooth_numbers(self, n: int, factor_base: List[int], required_count: int) -> List[Dict]:
        """
        Поиск B-гладких чисел методом пробного деления.
        """
        smooth_numbers = []
        # Начинаем с x = ⌊√n⌋ + 1, чтобы Q(x) = x² − n было положительным и малым.
        # Чем ближе x к √n, тем меньше Q(x) и тем выше вероятность гладкости.
        x_start = math.isqrt(n) + 1
        x = x_start
        table_data = []
        # Масштабируем лимит поиска под размер числа
        # Для больших чисел увеличиваем лимит поиска
        bit_length = n.bit_length()
        if bit_length <= 40:
            max_search = max(500_000, required_count * 5000)
        elif bit_length <= 60:
            max_search = max(1_000_000, required_count * 10000)
        elif bit_length <= 80:
            max_search = max(2_000_000, required_count * 20000)
        else:
            max_search = max(5_000_000, required_count * 50000)
        
        checked = 0
        # Множество простых из базы для быстрой проверки делимости
        fb_set = set(factor_base)
        # Берём только малые простые (≤ 13) для быстрого предфильтра
        small_primes = [p for p in factor_base if p <= 13]

        self.log_step("Этап 2: Поиск B-гладких чисел (пробное деление)", {
            "message": (
                f"Перебираем x = ⌊√n⌋ + 1, ⌊√n⌋ + 2, ...\n"
                f"Для каждого x вычисляем Q(x) = x² − n.\n"
                f"Q(x) мало по модулю (≈ 2x·Δ), поэтому шансы на гладкость высоки.\n"
                f"Число B-гладкое, если после деления на все p из базы остаток = 1.\n"
                f"Нужно найти: {required_count} гладких чисел (|FB| + 5).\n"
                f"Начинаем с x = {x_start}, лимит поиска: {max_search}"
            )
        })

        _sieve_start = time.perf_counter()
        # Логируем прогресс каждые 10%
        log_interval = max_search // 10
        last_log = 0
        
        while len(smooth_numbers) < required_count and (x - x_start) < max_search:
            # Вычисляем Q(x) = x² − n
            q_x = x * x - n
            temp_q = q_x
            checked += 1

            # Быстрый предфильтр
            if small_primes and not any(temp_q % p == 0 for p in small_primes):
                x += 1
                continue

            # Пробное деление
            exponents = [0] * len(factor_base)
            for i, p in enumerate(factor_base):
                while temp_q % p == 0:
                    exponents[i] += 1
                    temp_q //= p

            # Если B-гладкое
            if temp_q == 1:
                smooth_numbers.append({'x': x, 'q_x': q_x, 'exponents': exponents})
                # Логируем только первые 12 для экономии времени
                if len(smooth_numbers) <= 12:
                    table_data.append({
                        "x": str(x),
                        "Q(x) = x²−n": str(q_x),
                        "Вектор степеней": str(exponents[:8]) + ("..." if len(exponents) > 8 else ""),
                        "Гладкое?": "✓"
                    })
            
            # Логируем прогресс
            if checked - last_log >= log_interval:
                progress = (checked / max_search) * 100
                found_rate = (len(smooth_numbers) / checked) * 100 if checked > 0 else 0
                print(f"  Прогресс: {progress:.0f}% | Проверено: {checked} | Найдено: {len(smooth_numbers)} ({found_rate:.2f}%)")
                last_log = checked

            x += 1

        _sieve_ms = (time.perf_counter() - _sieve_start) * 1000
        self.log_step("Результат поиска гладких чисел", {
            "message": (
                f"Просмотрено значений x: {checked}\n"
                f"Найдено B-гладких чисел: {len(smooth_numbers)} / {required_count}\n"
                f"Время поиска гладких чисел: {_sieve_ms:.3f} мс ({_sieve_ms/1000:.1f} сек)\n"
                f"Каждая строка таблицы — это соотношение x² ≡ Q(x) (mod n),\n"
                f"где Q(x) полностью раскладывается по факторной базе."
            ),
            "table": table_data if table_data else None,
            "sieve_time_ms": round(_sieve_ms, 3),
            "checked": checked,
            "smooth_found": len(smooth_numbers),
            "required_smooth": required_count,
        })
        return smooth_numbers

    def _gauss_elimination_gf2(self, matrix: List[List[int]], smooth_numbers: List[Dict], factor_base: List[int]) -> tuple[List[List[int]], float]:
        """
        Решение системы линейных уравнений методом Гаусса над GF(2).
        
        Возвращает: (dependencies, time_ms)
        """
        _gauss_start = time.perf_counter()
        
        rows = len(matrix)
        if rows == 0 or len(matrix[0]) == 0:
            return [], 0.0
        cols = len(matrix[0])

        # Строим расширенную матрицу [A | I], где:
        #   A — исходная матрица над GF(2) (показатели степеней mod 2)
        #   I — единичная матрица того же числа строк
        # Правая часть I нужна для отслеживания, какие строки участвовали
        # в каждой линейной комбинации. После приведения Гаусса правая часть
        # нулевой строки покажет маску гладких чисел, образующих зависимость.
        M = [matrix[i] + [1 if i == j else 0 for j in range(rows)] for i in range(rows)]

        # Передаём срез матрицы (макс 30×30) для визуализации heatmap
        display_matrix = [row[:30] for row in matrix[:30]]

        # Строим таблицу: первые 8 строк исходной матрицы для наглядности
        matrix_preview = []
        for i, sn in enumerate(smooth_numbers[:8]):
            row_mod2 = [exp % 2 for exp in sn['exponents']]
            matrix_preview.append({
                "x": str(sn['x']),
                "Q(x)": str(sn['q_x']),
                "Вектор mod 2": str(row_mod2[:12]) + ("..." if len(row_mod2) > 12 else ""),
            })

        # Формируем примеры уравнений СЛАУ для первых 5 строк
        equations_examples = []
        for i in range(min(5, len(smooth_numbers))):
            sn = smooth_numbers[i]
            row_mod2 = [exp % 2 for exp in sn['exponents']]
            
            # Разложение Q(x) для отображения в LaTeX-формате
            factorization_parts = []
            for j, exp in enumerate(sn['exponents'][:min(6, len(sn['exponents']))]):
                if exp > 0:
                    p = factor_base[j]
                    if exp == 1:
                        factorization_parts.append(f"{p}")
                    else:
                        factorization_parts.append(f"{p}^{{{exp}}}")
            
            if len(sn['exponents']) > 6:
                factorization_parts.append(r"\cdots")
            
            factorization_str = r" \cdot ".join(factorization_parts) if factorization_parts else "1"
            
            equations_examples.append({
                "i": i + 1,
                "x_i": str(sn['x']),
                "q_x": str(sn['q_x']),
                "factorization_latex": factorization_str,
                "row_mod2": row_mod2,  # передаём сырые коэффициенты
            })

        self.log_step("Этап 3: Построение матрицы над GF(2)", {
            "message": (
                f"Каждое B-гладкое число Q(xᵢ) = ∏ pⱼ^eᵢⱼ даёт строку матрицы:\n"
                f"  aᵢⱼ = eᵢⱼ mod 2  (чётная степень → 0, нечётная → 1)\n\n"
                f"Матрица A размером {rows} × {cols}:\n"
                f"  строки  = {rows} найденных B-гладких чисел\n"
                f"  столбцы = {cols} простых из факторной базы\n\n"
                f"Задача: найти ненулевой вектор v ∈ GF(2)^{rows} такой, что\n"
                f"  v · A ≡ 0 (mod 2)\n"
                f"Это означает: произведение выбранных Q(xᵢ) — полный квадрат,\n"
                f"т.е. все показатели простых в нём чётные."
            ),
            "matrix_data": display_matrix,
            "table": matrix_preview,
            "equations_examples": equations_examples,
            "factor_base_first": factor_base[:min(10, len(factor_base))],
        })

        self.log_step("Этап 3: Расширенная матрица [A | I]", {
            "message": (
                f"Для отслеживания комбинаций строк расширяем матрицу единичной:\n\n"
                f"  [A | I]  размер: {rows} × {cols + rows}\n\n"
                f"Правая часть I (единичная матрица {rows}×{rows}) кодирует,\n"
                f"какие исходные строки участвуют в текущей комбинации.\n\n"
                f"После приведения к ступенчатому виду:\n"
                f"  • Ненулевые строки слева → линейно независимые соотношения\n"
                f"  • Нулевые строки слева   → линейные зависимости\n"
                f"    Правая часть нулевой строки = маска участвующих гладких чисел"
            )
        })

        pivot_row = 0
        gauss_steps = []  # лог ключевых шагов исключения для отображения в UI

        # Прямой ход Гаусса: приводим матрицу к ступенчатому виду над GF(2).
        # Обрабатываем столбцы слева направо — каждый столбец соответствует
        # одному простому из факторной базы.
        for c in range(cols):
            # Ищем опорную строку (pivot): первую строку начиная с pivot_row,
            # у которой в столбце c стоит 1. Если таких нет — столбец нулевой,
            # пропускаем (это означает, что данный простой не встречается ни в
            # одном ещё не обработанном гладком числе).
            pivot = next((r for r in range(pivot_row, rows) if M[r][c] == 1), -1)
            if pivot == -1:
                continue
            # Переставляем опорную строку на позицию pivot_row,
            # чтобы она стала "ведущей" для текущего столбца.
            if pivot != pivot_row:
                M[pivot_row], M[pivot] = M[pivot], M[pivot_row]
                if len(gauss_steps) < 10:
                    gauss_steps.append({
                        "Действие": f"Перестановка строк",
                        "Столбец (опорный)": c,
                        "Строка pivot": pivot_row,
                        "Комментарий": f"строка {pivot} ↔ строка {pivot_row}",
                    })
            eliminated = 0
            # XOR-исключение: для каждой строки r (кроме опорной), где M[r][c] = 1,
            # выполняем строка[r] = строка[r] XOR строка[pivot_row].
            # В GF(2) это эквивалентно вычитанию — обнуляет элемент в столбце c.
            # Операция применяется ко всей строке, включая правую часть [I],
            # что позволяет отслеживать, какие исходные строки вошли в комбинацию.
            for r in range(rows):
                if r != pivot_row and M[r][c] == 1:
                    M[r] = [M[r][i] ^ M[pivot_row][i] for i in range(len(M[0]))]
                    eliminated += 1
            if len(gauss_steps) < 10:
                gauss_steps.append({
                    "Действие": "XOR-исключение",
                    "Столбец (опорный)": c,
                    "Строка pivot": pivot_row,
                    "Комментарий": f"обнулено строк: {eliminated}",
                })
            pivot_row += 1

        # Извлекаем линейные зависимости: строки с нулевой левой частью [A|...]
        # означают, что соответствующая комбинация гладких чисел даёт полный квадрат.
        # Правая часть (маска из единичной матрицы) показывает, какие именно строки
        # участвовали в этой комбинации — это и есть зависимость.
        dependencies = [M[r][cols:] for r in range(pivot_row, rows)]

        _gauss_ms = (time.perf_counter() - _gauss_start) * 1000

        # Таблица найденных зависимостей (первые 5)
        dep_table = []
        for i, dep in enumerate(dependencies[:5]):
            used = [j for j, v in enumerate(dep) if v == 1]
            dep_table.append({
                "Зависимость #": i + 1,
                "Индексы строк": str(used),
                "Кол-во слагаемых": len(used),
            })

        self.log_step("Этап 3: Ход метода Гаусса над GF(2)", {
            "message": (
                f"Алгоритм прямого хода (приведение к ступенчатому виду):\n\n"
                f"Для каждого столбца c = 0..{cols-1}:\n"
                f"  1. Ищем опорную строку: первую строку ≥ pivot_row, где M[r][c] = 1\n"
                f"  2. Если не найдена → столбец нулевой, пропускаем\n"
                f"  3. Меняем опорную строку с текущей pivot_row местами\n"
                f"  4. Для всех остальных строк r, где M[r][c] = 1:\n"
                f"       строка[r] = строка[r] XOR строка[pivot_row]\n"
                f"     (это обнуляет элемент в столбце c для строки r)\n"
                f"  5. pivot_row += 1\n\n"
                f"Сложение в GF(2): 0⊕0=0, 0⊕1=1, 1⊕0=1, 1⊕1=0\n"
                f"(эквивалентно XOR — нет переноса, нет деления)\n\n"
                f"Итог: ранг матрицы = {pivot_row}, обработано столбцов = {cols}"
            ),
            "table": gauss_steps if gauss_steps else None,
        })

        self.log_step("Результат Гаусса: линейные зависимости", {
            "message": (
                f"После приведения к ступенчатому виду:\n\n"
                f"  Ранг матрицы A: {pivot_row}\n"
                f"  Строк всего:    {rows}\n"
                f"  Зависимостей:   {len(dependencies)}  (= {rows} − {pivot_row})\n\n"
                f"Каждая нулевая строка в левой части [A|I] означает:\n"
                f"  ∑ aᵢ · (строка i) ≡ 0 (mod 2)  для некоторого набора i\n"
                f"  ⟹  ∏ Q(xᵢ) = Y²  (полный квадрат)\n"
                f"  ⟹  (∏ xᵢ)² ≡ Y² (mod n)\n"
                f"  ⟹  НОД(∏xᵢ − Y, n) — кандидат на делитель\n\n"
                f"Чем больше зависимостей, тем выше шанс найти нетривиальный делитель.\n\n"
                f"Время решения СЛАУ: {_gauss_ms:.3f} мс"
            ),
            "table": dep_table if dep_table else None,
            "gauss_time_ms": round(_gauss_ms, 3),
        })
        return dependencies, _gauss_ms

    def factorize(self, n: int, b_override: int = None) -> List[int]:
        # Очищаем лог шагов от предыдущих запусков
        self.clear_logs()

        # Граничные случаи: числа ≤ 1, чётные и простые не требуют алгоритма
        if n <= 1: return [n]
        if n == 2: return [2]
        
        # Рекурсивная обработка чётных чисел
        if n % 2 == 0:
            # Факторизуем n//2 рекурсивно
            remaining_factors = self.factorize(n // 2, b_override)
            return sorted([2] + remaining_factors)
        
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

        # ── Выбор границы B ──────────────────────────────────────────────────
        # B — граница гладкости: в факторную базу войдут только простые p ≤ B.
        # Оптимальное B по L-нотации: exp(0.5 · √(ln n · ln ln n)).
        # Слишком малый B → мало гладких чисел, алгоритм не найдёт зависимость.
        # Слишком большой B → большая матрица Гаусса, медленное решение СЛАУ.
        B_auto = int(math.exp(0.5 * math.sqrt(math.log(n) * math.log(math.log(n)))))
        B_auto = max(B_auto, 100)    # нижний порог: слишком малый B бесполезен
        B_auto = min(B_auto, 50000)  # верхний порог: ограничиваем память и время

        if b_override is not None:
            # Пользователь задал B вручную — используем его значение
            B = max(10, int(b_override))
            mode_label = "ручной (режим исследователя)"
            mode_note = (
                f"⚠ Задан вручную: B = {B}\n"
                f"Автоматическое значение было бы: B_авто = {B_auto}\n\n"
                f"{'B < B_авто → факторная база мала, гладких чисел может не хватить.' if B < B_auto else ''}"
                f"{'B > B_авто → факторная база велика, матрица Гаусса будет большой.' if B > B_auto else ''}"
            )
        else:
            # Автоматический режим: используем оптимальное B по формуле
            B = B_auto
            mode_label = "автоматический (L-нотация)"
            mode_note = ""

        pi_B = len(generate_primes(B))
        pi_B_approx = B / math.log(B) if B > 1 else 0

        self.log_step("Выбор параметра B", {
            "message": (
                f"Режим: {mode_label}\n\n"
                f"Оптимальная граница по L-нотации:\n"
                f"B = exp(0.5 · √(ln n · ln ln n)) ≈ {B_auto}\n\n"
                + (f"{mode_note}\n\n" if mode_note else "") +
                f"Размер факторной базы определяется функцией π(B) —\n"
                f"количеством простых чисел, не превышающих B.\n\n"
                f"По теореме о простых числах: π(B) ≈ B / ln(B)\n"
                f"π({B}) ≈ {B} / ln({B}) ≈ {B} / {math.log(B):.2f} ≈ {pi_B_approx:.0f}\n"
                f"Точное значение: π({B}) = {pi_B} простых чисел\n\n"
                f"Слишком малый B → мало гладких чисел.\n"
                f"Слишком большой B → большая факторная база, медленный Гаусс."
            ),
            "pi_B": pi_B,
            "B": B,
            "B_auto": B_auto,
            "b_mode": "manual" if b_override is not None else "auto",
            "pi_B_data": {
                "B": B,
                "B_auto": B_auto,
                "pi_B_exact": pi_B,
                "pi_B_approx": round(pi_B_approx, 1),
                "ln_B": round(math.log(B), 4),
                "b_mode": "manual" if b_override is not None else "auto",
            },
        })

        factor_base, fb_time_ms = self._get_factor_base(n, B)
        # Нужно найти на 5 гладких чисел больше, чем размер базы:
        # по теореме линейной алгебры, если строк > столбцов, гарантированно
        # найдётся хотя бы одна линейная зависимость над GF(2).
        required_smooth = len(factor_base) + 5
        smooth_numbers = self._find_smooth_numbers(n, factor_base, required_smooth)
        
        # Извлекаем sieve_time_ms из последнего шага лога (записан внутри _find_smooth_numbers)
        sieve_time_ms = self.steps_log[-1]["details"].get("sieve_time_ms", 0)

        if not smooth_numbers:
            self.log_step("Провал", {
                "message": (
                    "Не найдено ни одного B-гладкого числа.\n"
                    "Решение: увеличить B для расширения факторной базы."
                )
            })
            return [n]

        # Строим матрицу над GF(2): каждая строка — вектор показателей степеней
        # одного гладкого числа Q(xᵢ), взятых по модулю 2.
        # Чётная степень → 0 (не влияет на квадратичность), нечётная → 1.
        matrix_mod2 = [[exp % 2 for exp in sn['exponents']] for sn in smooth_numbers]
        dependencies, gauss_time_ms = self._gauss_elimination_gf2(matrix_mod2, smooth_numbers, factor_base)

        # ── Профилирование времени ───────────────────────────────────────────
        # Суммируем время трёх основных этапов для отображения разбивки
        total_time_ms = fb_time_ms + sieve_time_ms + gauss_time_ms
        other_time_ms = max(0, total_time_ms * 0.02)  # ~2% на остальное (проверка зависимостей и т.д.)
        
        self.log_step("⏱ Профилирование времени выполнения", {
            "message": (
                f"Разбивка времени работы алгоритма Диксона:\n\n"
                f"1. Построение факторной базы:  {fb_time_ms:>8.3f} мс  ({fb_time_ms/total_time_ms*100:>5.1f}%)\n"
                f"2. Поиск гладких чисел (Сито): {sieve_time_ms:>8.3f} мс  ({sieve_time_ms/total_time_ms*100:>5.1f}%)\n"
                f"3. Решение СЛАУ (метод Гаусса): {gauss_time_ms:>8.3f} мс  ({gauss_time_ms/total_time_ms*100:>5.1f}%)\n"
                f"4. Прочее (проверка зависимостей): {other_time_ms:>8.3f} мс  ({other_time_ms/total_time_ms*100:>5.1f}%)\n"
                f"{'─' * 60}\n"
                f"Итого:                          {total_time_ms:>8.3f} мс  (100.0%)\n\n"
                f"Узкое место: {'Сито' if sieve_time_ms == max(fb_time_ms, sieve_time_ms, gauss_time_ms) else 'СЛАУ' if gauss_time_ms > sieve_time_ms else 'FB'}"
            ),
            "profiling": {
                "Построение FB": round(fb_time_ms, 3),
                "Поиск гладких (Сито)": round(sieve_time_ms, 3),
                "Решение СЛАУ (Гаусс)": round(gauss_time_ms, 3),
                "Прочее": round(other_time_ms, 3),
            }
        })

        self.log_step("Этап 4: Проверка зависимостей", {
            "message": (
                f"Для каждой зависимости вычисляем:\n"
                f"  X = ∏ xᵢ (mod n)  — левая часть\n"
                f"  Y = ∏ pⱼ^(eⱼ/2) (mod n)  — квадратный корень из ∏ Q(xᵢ)\n"
                f"Проверяем: НОД(X−Y, n) и НОД(X+Y, n)"
            )
        })

        for idx, dep in enumerate(dependencies):
            # dep — бинарная маска: dep[i] = 1 означает, что i-е гладкое число
            # участвует в данной зависимости.

            X = 1               # X = ∏ xᵢ (mod n) — левая часть конгруэнции
            exponents_sum = [0] * len(factor_base)  # суммарные показатели степеней ∏ Q(xᵢ)
            used_xs = []        # список xᵢ для отображения в логе

            for i, is_used in enumerate(dep):
                if is_used:
                    # Накапливаем произведение xᵢ по модулю n
                    X = (X * smooth_numbers[i]['x']) % n
                    used_xs.append(smooth_numbers[i]['x'])
                    # Суммируем показатели степеней: ∏ Q(xᵢ) = ∏ pⱼ^(∑ eᵢⱼ)
                    # Поскольку зависимость гарантирует чётность всех ∑ eᵢⱼ,
                    # можно взять квадратный корень: Y = ∏ pⱼ^(∑ eᵢⱼ / 2)
                    for j, exp in enumerate(smooth_numbers[i]['exponents']):
                        exponents_sum[j] += exp

            # Y = √(∏ Q(xᵢ)) mod n — правая часть конгруэнции X² ≡ Y² (mod n)
            Y = 1
            for i, p in enumerate(factor_base):
                # Делим каждый показатель на 2 (он чётный по построению зависимости)
                Y = (Y * pow(p, exponents_sum[i] // 2, n)) % n

            # Проверяем НОД(X−Y, n) и НОД(X+Y, n).
            # Если X² ≡ Y² (mod n), то n | (X−Y)(X+Y).
            # Нетривиальный делитель: 1 < НОД < n (т.е. X ≢ ±Y (mod n)).
            d1 = math.gcd(abs(X - Y), n)
            d2 = math.gcd(abs(X + Y), n)
            # Выбираем нетривиальный делитель, если он есть
            d = d1 if 1 < d1 < n else (d2 if 1 < d2 < n else d1)

            self.log_step(f"Зависимость #{idx + 1}", {
                "message": (
                    f"Использованы x: {[str(v) for v in used_xs[:8]]}{'...' if len(used_xs) > 8 else ''}\n"
                    f"X = {X}, Y = {Y}\n"
                    f"НОД(|X−Y|, n) = {d1},  НОД(|X+Y|, n) = {d2}\n"
                    f"{'✓ Нетривиальный делитель!' if 1 < d < n else '✗ Тривиальный (1 или n), пробуем следующую зависимость.'}"
                )
            })

            if 1 < d < n:
                self.log_step("Факторизация завершена", {
                    "message": f"{n} = {d} × {n // d}"
                })
                return sorted([d, n // d])

        # ── Повторная попытка с увеличенным B ───────────────────────────────
        # Если все зависимости дали тривиальные делители (X ≡ ±Y mod n),
        # это не ошибка алгоритма — просто не повезло с комбинациями.
        # Увеличиваем B в 3 раза: больше гладких чисел → больше зависимостей →
        # выше шанс найти нетривиальный делитель.
        if b_override is None and B < 50000:
            B2 = min(B * 3, 50000)
            self.log_step("Повторная попытка с B×3", {
                "message": (
                    f"Все зависимости дали тривиальные делители.\n"
                    f"Автоматически увеличиваем B: {B} → {B2} и повторяем поиск."
                )
            })
            factor_base2, _ = self._get_factor_base(n, B2)
            required2 = len(factor_base2) + 10  # берём чуть больше запаса для надёжности
            smooth2 = self._find_smooth_numbers(n, factor_base2, required2)
            if smooth2:
                matrix2 = [[exp % 2 for exp in sn['exponents']] for sn in smooth2]
                deps2, _ = self._gauss_elimination_gf2(matrix2, smooth2, factor_base2)
                for dep in deps2:
                    X = 1
                    exponents_sum = [0] * len(factor_base2)
                    for i, is_used in enumerate(dep):
                        if is_used:
                            X = (X * smooth2[i]['x']) % n
                            for j, exp in enumerate(smooth2[i]['exponents']):
                                exponents_sum[j] += exp
                    Y = 1
                    for i, p in enumerate(factor_base2):
                        Y = (Y * pow(p, exponents_sum[i] // 2, n)) % n
                    # Проверяем оба кандидата на делитель
                    for candidate in [math.gcd(abs(X - Y), n), math.gcd(abs(X + Y), n)]:
                        if 1 < candidate < n:
                            self.log_step("Факторизация завершена (повторная попытка)", {
                                "message": f"{n} = {candidate} × {n // candidate}"
                            })
                            return sorted([candidate, n // candidate])

        self.log_step("Провал", {
            "message": (
                "Все зависимости дали тривиальные делители.\n"
                "Причина: недостаточно гладких чисел или неудачные комбинации.\n"
                "Решение: увеличить B для расширения факторной базы."
            )
        })
        return [n]
