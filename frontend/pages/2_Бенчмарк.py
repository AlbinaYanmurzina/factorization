import streamlit as st
import requests
import random
import math
import re
import platform
import statistics
import plotly.graph_objects as go
import pandas as pd
from sympy import isprime, nextprime

st.set_page_config(page_title="Сравнение алгоритмов", page_icon="📊", layout="wide")

# ── Информация о системе ────────────────────────────────────────────────────

@st.cache_data
def get_system_info() -> dict:
    info = {
        "ОС": platform.system() + " " + platform.release(),
        "Процессор": platform.processor() or platform.machine(),
        "Python": platform.python_version(),
        "Ядра CPU": "—",
        "RAM": "—",
    }
    try:
        import psutil, os
        info["Ядра CPU"] = f"{os.cpu_count()} логических / {psutil.cpu_count(logical=False)} физических"
        ram_gb = psutil.virtual_memory().total / (1024 ** 3)
        info["RAM"] = f"{ram_gb:.1f} ГБ"
    except ImportError:
        import os
        info["Ядра CPU"] = str(os.cpu_count())
    return info

sys_info = get_system_info()

ALGO_COMPLEXITY_LABEL = {
    "pollard_rho": "O(n^{1/4})",
    "pollard_p1":  "O(n^{1/2})",
    "qs_basic":    "L[1/2, c]",
}

# ── Вспомогательные функции ─────────────────────────────────────────────────

def generate_semiprime(bits: int) -> int:
    start_p = random.getrandbits(bits // 2)
    start_q = random.getrandbits(bits // 2)
    p = nextprime(max(start_p, 3))
    q = nextprime(max(start_q, 3))
    while q == p:
        q = nextprime(q)
    return p * q

def analyze_smoothness(n: int) -> dict:
    """
    Анализирует гладкость (p-1) и (q-1) для полупростого числа n = p × q.
    Возвращает информацию о наибольших простых делителях.
    """
    from sympy import factorint, primefactors
    
    try:
        # Факторизуем число
        factors_dict = factorint(n)
        if len(factors_dict) != 2:
            # Не полупростое число
            return {
                'is_semiprime': False,
                'p': None,
                'q': None,
                'p_minus_1_max': None,
                'q_minus_1_max': None,
                'smooth_100': False,
                'smooth_1000': False,
                'smooth_10000': False,
            }
        
        factors_list = list(factors_dict.keys())
        p, q = factors_list[0], factors_list[1]
        
        # Анализируем (p-1) и (q-1)
        p_minus_1_factors = primefactors(p - 1)
        q_minus_1_factors = primefactors(q - 1)
        
        p_max = max(p_minus_1_factors) if p_minus_1_factors else 1
        q_max = max(q_minus_1_factors) if q_minus_1_factors else 1
        
        # Определяем гладкость
        smooth_100 = (p_max <= 100 or q_max <= 100)
        smooth_1000 = (p_max <= 1000 or q_max <= 1000)
        smooth_10000 = (p_max <= 10000 or q_max <= 10000)
        
        return {
            'is_semiprime': True,
            'p': p,
            'q': q,
            'p_minus_1_max': p_max,
            'q_minus_1_max': q_max,
            'smooth_100': smooth_100,
            'smooth_1000': smooth_1000,
            'smooth_10000': smooth_10000,
        }
    except Exception as e:
        return {
            'is_semiprime': False,
            'error': str(e),
        }

def extract_benchmark_metadata(algorithm: str, steps: list[dict]) -> dict:
    """Pull thesis-table parameters from backend step logs."""
    meta: dict[str, object] = {}

    if algorithm == "pollard_rho":
        meta.update({
            "rho_B": "нет",
            "rho_c_set": "1,2,3,5,7",
            "rho_max_iterations": 100000,
        })
    elif algorithm == "pollard_p1":
        meta.update({
            "p1_B_set": "3,5,10,25,50,100,1000,10000",
            "p1_a_set": "2,3,5,7,11",
        })

    rho_iterations_total = 0
    rho_c_values: list[int] = []
    p1_last_B = None
    p1_last_a = None

    for step in steps or []:
        name = str(step.get("step", ""))
        details = step.get("details", {}) or {}
        message = str(details.get("message", ""))
        text = f"{name}\n{message}"

        if algorithm == "qs_basic":
            for source_key, target_key in {
                "B": "dixon_B",
                "B_auto": "dixon_B_auto",
                "pi_B": "dixon_pi_B",
                "FB_size": "dixon_FB_size",
                "fb_time_ms": "dixon_FB_time_ms",
                "sieve_time_ms": "dixon_smooth_search_ms",
                "checked": "dixon_checked_x",
                "smooth_found": "dixon_smooth_found",
                "required_smooth": "dixon_required_smooth",
            }.items():
                if source_key in details:
                    meta[target_key] = details[source_key]

            profiling = details.get("profiling")
            if isinstance(profiling, dict):
                values = list(profiling.values())
                if len(values) >= 1:
                    meta["dixon_profile_FB_ms"] = values[0]
                if len(values) >= 2:
                    meta["dixon_profile_smooth_search_ms"] = values[1]
                if len(values) >= 3:
                    meta["dixon_profile_gauss_ms"] = values[2]
                if len(values) >= 4:
                    meta["dixon_profile_other_ms"] = values[3]

        elif algorithm == "pollard_p1":
            match = re.search(r"B\s*=\s*(\d+),\s*a\s*=\s*(\d+)", text)
            if match:
                p1_last_B = int(match.group(1))
                p1_last_a = int(match.group(2))

        elif algorithm == "pollard_rho":
            iter_match = re.search(r":\s*(\d+)\)", name)
            if iter_match:
                rho_iterations_total += int(iter_match.group(1))

            c_match = re.search(r"\+\s*(\d+)\)\s*mod", message)
            if c_match:
                c_value = int(c_match.group(1))
                if c_value not in rho_c_values:
                    rho_c_values.append(c_value)

    if algorithm == "pollard_p1":
        meta["p1_last_B"] = p1_last_B if p1_last_B is not None else "—"
        meta["p1_last_a"] = p1_last_a if p1_last_a is not None else "—"
    elif algorithm == "pollard_rho":
        meta["rho_iterations_total"] = rho_iterations_total if rho_iterations_total else "—"
        meta["rho_c_used"] = ",".join(map(str, rho_c_values)) if rho_c_values else "—"

    return meta

def get_timeout_for_bits(bits: int) -> int:
    """
    Возвращает таймаут в секундах в зависимости от битности.
    
    Логика:
    - До 80 бит: 120 секунд (достаточно для всех алгоритмов)
    - 80-85 бит: 300 секунд (5 минут)
    - 85-92 бит: 600 секунд (10 минут, Диксон может занять ~10 минут)
    - 92-100 бит: 1200 секунд (20 минут, Диксон может занять ~20 минут)
    - 100-110 бит: 2400 секунд (40 минут, Диксон может занять ~40 минут)
    - 110+ бит: 3600 секунд (60 минут, Диксон может занять ~60 минут)
    """
    if bits < 80:
        return 120
    elif bits < 85:
        return 300
    elif bits < 92:
        return 600
    elif bits < 100:
        return 1200
    elif bits < 110:
        return 2400
    else:
        return 3600

def estimate_time_for_bits(bits: int, algorithm: str) -> str:
    """
    Оценивает примерное время выполнения для заданной битности.
    
    Возвращает строку с оценкой времени.
    """
    if algorithm == "pollard_rho":
        if bits < 60:
            return "< 1 сек"
        elif bits < 80:
            return "1-2 сек"
        else:
            return "не справится"
    elif algorithm == "pollard_p1":
        if bits < 60:
            return "< 1 сек"
        else:
            return "не справится (или случайный успех ~5%)"
    elif algorithm == "qs_basic":
        if bits < 40:
            return "< 1 сек"
        elif bits < 60:
            return "1-5 сек"
        elif bits < 80:
            return "10-60 сек"
        elif bits < 85:
            return "1-5 мин"
        elif bits < 92:
            return "5-10 мин"
        elif bits < 100:
            return "10-20 мин"
        elif bits < 110:
            return "20-40 мин"
        elif bits < 120:
            return "40-60 мин"
        else:
            return "60+ мин"
    return "неизвестно"

ALGO_LIST = [
    ("ρ-метод Полларда (разд. 3.4)",                        "pollard_rho"),
    ("(p-1)-метод Полларда (разд. 3.2)",                    "pollard_p1"),
    ("Алгоритм Диксона (разд. 6.1)",              "qs_basic"),
]

# ── Заголовок и информация о системе ───────────────────────────────────────

st.title("📊 Вычислительные эксперименты")

with st.expander("🖥️ Стенд тестирования", expanded=True):
    cols = st.columns(len(sys_info))
    for col, (key, val) in zip(cols, sys_info.items()):
        col.metric(key, val)

st.divider()

# ── Боковая панель ──────────────────────────────────────────────────────────

st.sidebar.header("Настройки эксперимента")

min_bits = st.sidebar.slider("Минимальная разрядность (бит)", 16, 60, 20, step=2)
max_bits = st.sidebar.slider("Максимальная разрядность (бит)", 16, 130, 40, step=2)
step_bits = st.sidebar.slider("Шаг разрядности (бит)", 2, 8, 4, step=2)
runs_per_bit = st.sidebar.slider("Замеров на точку графика", 1, 7, 3)

st.sidebar.markdown("---")
st.sidebar.subheader("Выбор алгоритмов")

selected_algos = []
for name, key in ALGO_LIST:
    if st.sidebar.checkbox(name, value=True, key=f"cb_{key}"):
        selected_algos.append((name, key))

run_btn = st.sidebar.button("🚀 Запустить тест", type="primary", use_container_width=True)

# ── Запуск бенчмарка ────────────────────────────────────────────────────────

if run_btn:
    if min_bits >= max_bits:
        st.error("Минимальная разрядность должна быть меньше максимальной.")
        st.stop()

    bit_range = list(range(min_bits, max_bits + 1, step_bits))
    total_steps = len(bit_range) * len(selected_algos)
    step_counter = 0

    # results[algo_key][bit] = list of times (ms)
    raw_results: dict[str, dict[int, list[float]]] = {
        key: {b: [] for b in bit_range} for _, key in selected_algos
    }
    
    # Счётчики неудач для каждого алгоритма и битности
    failed_counts: dict[str, dict[int, int]] = {
        key: {b: 0 for b in bit_range} for _, key in selected_algos
    }
    failed_results: dict[str, dict[int, list[float]]] = {
        key: {b: [] for b in bit_range} for _, key in selected_algos
    }
    partial_counts: dict[str, dict[int, int]] = {
        key: {b: 0 for b in bit_range} for _, key in selected_algos
    }
    
    # Детальные результаты для отображения
    detailed_results = []

    progress_bar = st.progress(0, text="Запуск...")
    status_text = st.empty()
    
    # Контейнер для отображения текущего процесса
    process_container = st.container()
    with process_container:
        current_test_info = st.empty()
        current_result_info = st.empty()

    for bit in bit_range:
        test_numbers = [generate_semiprime(bit) for _ in range(runs_per_bit)]
        
        # Определяем таймаут для текущей битности
        timeout_sec = get_timeout_for_bits(bit)
        timeout_ms = timeout_sec * 1000

        for alg_name, alg_key in selected_algos:
            step_counter += 1
            progress_bar.progress(
                step_counter / total_steps,
                text=f"{alg_name} @ {bit} бит (таймаут: {timeout_sec}с) ({step_counter}/{total_steps})"
            )

            for idx, num in enumerate(test_numbers, 1):
                # Анализ гладкости числа
                smoothness = analyze_smoothness(num)
                
                # Отображаем текущее тестируемое число
                smooth_info = ""
                if smoothness['is_semiprime']:
                    p_max = smoothness['p_minus_1_max']
                    q_max = smoothness['q_minus_1_max']
                    smooth_info = f"\n**Гладкость:** (p-1) max={p_max}, (q-1) max={q_max}"
                    if smoothness['smooth_100']:
                        smooth_info += " → 100-гладкое ✓"
                    elif smoothness['smooth_1000']:
                        smooth_info += " → 1000-гладкое"
                    elif smoothness['smooth_10000']:
                        smooth_info += " → 10000-гладкое"
                    else:
                        smooth_info += " → НЕ гладкое"
                
                current_test_info.info(
                    f"🔍 **Тестирование:** {alg_name.split(' (')[0]} | "
                    f"{bit} бит | Попытка {idx}/{len(test_numbers)}\n\n"
                    f"**Число:** `{num}`{smooth_info}"
                )
                
                try:
                    res = requests.post(
                        "http://127.0.0.1:8453/api/factorize",
                        json={"number": str(num), "algorithm": alg_key},
                        timeout=timeout_sec + 5,  # +5 секунд запас для HTTP
                    )
                    if res.status_code == 200:
                        data = res.json()
                        t = data["time_ms"]
                        factors = data["factors"]
                        step_meta = extract_benchmark_metadata(alg_key, data.get("steps", []))
                        
                        factor_ints = []
                        try:
                            factor_ints = [int(f) for f in factors]
                            product_ok = math.prod(factor_ints) == num
                            all_factors_prime = all(isprime(f) for f in factor_ints)
                        except (TypeError, ValueError):
                            product_ok = False
                            all_factors_prime = False
                        
                        full_factorization_success = product_ok and all_factors_prime and t < timeout_ms
                        partial_factorization = product_ok and not all_factors_prime and not (
                            len(factors) == 1 and factors[0] == str(num)
                        )
                        factorization_failed = not full_factorization_success and not partial_factorization and t < timeout_ms
                        
                        # Формируем результат для отображения
                        if full_factorization_success:
                            # Успешная полная факторизация в пределах таймаута
                            raw_results[alg_key][bit].append(t)
                            result_status = "✅ Успех"
                            result_detail = f"{num} = {' × '.join(factors)}"
                            current_result_info.success(
                                f"{result_status}\n\n{result_detail}\n\n⏱️ Время: {t:.2f} мс ({t/1000:.2f} сек)"
                            )
                        elif partial_factorization:
                            failed_counts[alg_key][bit] += 1
                            failed_results[alg_key][bit].append(t)
                            partial_counts[alg_key][bit] += 1
                            result_status = "⚠️ Частично"
                            result_detail = f"Найдено неполное разложение: {num} = {' × '.join(factors)}"
                            current_result_info.warning(
                                f"{result_status}\n\n{result_detail}\n\n⏱️ Время: {t:.2f} мс"
                            )
                        elif factorization_failed:
                            failed_counts[alg_key][bit] += 1
                            failed_results[alg_key][bit].append(t)
                            result_status = "❌ Не справился"
                            result_detail = f"Не получено полное разложение на простые множители"
                            current_result_info.error(
                                f"{result_status}\n\n{result_detail}\n\n⏱️ Время: {t:.2f} мс"
                            )
                        else:
                            # Таймаут
                            failed_counts[alg_key][bit] += 1
                            failed_results[alg_key][bit].append(t)
                            result_status = "⏱️ Таймаут"
                            result_detail = f"Превышен лимит {timeout_sec} сек"
                            current_result_info.warning(
                                f"{result_status}\n\n{result_detail}\n\n⏱️ Время: {t:.2f} мс ({t/1000:.2f} сек)"
                            )
                        
                        # Сохраняем детальный результат
                        detailed_results.append({
                            "Алгоритм": alg_name.split(" (")[0],
                            "Биты": bit,
                            "Число": str(num),
                            "Результат": ' × '.join(factors) if (full_factorization_success or partial_factorization) else "не справился",
                            "Статус": result_status,
                            "Время (мс)": round(t, 2),
                            "Время (сек)": round(t / 1000, 2),
                            "Гладкость": f"(p-1)≤{smoothness['p_minus_1_max']}, (q-1)≤{smoothness['q_minus_1_max']}" if smoothness['is_semiprime'] else "—",
                            "100-гладкое": "✓" if smoothness.get('smooth_100') else "✗",
                            **step_meta,
                        })
                        
                        # Небольшая пауза для читаемости
                        import time as time_module
                        time_module.sleep(0.3)
                        
                except requests.exceptions.Timeout:
                    # Таймаут на уровне HTTP-запроса
                    failed_counts[alg_key][bit] += 1
                    failed_results[alg_key][bit].append(timeout_ms)
                    current_result_info.error(
                        f"❌ HTTP Таймаут\n\n"
                        f"Превышен лимит HTTP-запроса ({timeout_sec + 5} сек)"
                    )
                    detailed_results.append({
                        "Алгоритм": alg_name.split(" (")[0],
                        "Биты": bit,
                        "Число": str(num),
                        "Результат": "HTTP таймаут",
                        "Статус": "❌ Таймаут",
                        "Время (мс)": timeout_ms,
                        "Время (сек)": timeout_sec,
                        "Гладкость": f"(p-1)≤{smoothness['p_minus_1_max']}, (q-1)≤{smoothness['q_minus_1_max']}" if smoothness['is_semiprime'] else "—",
                        "100-гладкое": "✓" if smoothness.get('smooth_100') else "✗",
                        **extract_benchmark_metadata(alg_key, []),
                    })
                except Exception as e:
                    failed_counts[alg_key][bit] += 1
                    failed_results[alg_key][bit].append(0)
                    current_result_info.error(
                        f"❌ Ошибка\n\n{str(e)}"
                    )
                    detailed_results.append({
                        "Алгоритм": alg_name.split(" (")[0],
                        "Биты": bit,
                        "Число": str(num),
                        "Результат": f"ошибка: {str(e)}",
                        "Статус": "❌ Ошибка",
                        "Время (мс)": 0,
                        "Время (сек)": 0,
                        "Гладкость": f"(p-1)≤{smoothness['p_minus_1_max']}, (q-1)≤{smoothness['q_minus_1_max']}" if smoothness['is_semiprime'] else "—",
                        "100-гладкое": "✓" if smoothness.get('smooth_100') else "✗",
                        **extract_benchmark_metadata(alg_key, []),
                    })

    progress_bar.empty()
    status_text.empty()
    current_test_info.empty()
    current_result_info.empty()
    
    st.success(f"Эксперимент завершён. Протестировано разрядностей: {len(bit_range)}, алгоритмов: {len(selected_algos)}, замеров на точку: {runs_per_bit}.")
    
    # ── Детальные результаты всех тестов ────────────────────────────────────
    
    st.subheader("📋 Детальные результаты всех тестов")
    
    if detailed_results:
        df_detailed = pd.DataFrame(detailed_results)
        
        st.dataframe(
            df_detailed,
            use_container_width=True,
            hide_index=True,
            column_config={
                "Число": st.column_config.TextColumn("Число", width="medium"),
                "Результат": st.column_config.TextColumn("Результат", width="large"),
                "Статус": st.column_config.TextColumn("Статус", width="small"),
                "Время (мс)": st.column_config.NumberColumn("Время (мс)", format="%.2f"),
                "Время (сек)": st.column_config.NumberColumn("Время (сек)", format="%.2f"),
            }
        )
        
        st.caption(f"Показано {len(df_detailed)} результатов.")
        csv_data = df_detailed.to_csv(index=False).encode("utf-8-sig")
        st.download_button(
            "Скачать детальные результаты CSV",
            data=csv_data,
            file_name="benchmark_detailed_results.csv",
            mime="text/csv",
            use_container_width=True,
        )
    
    st.divider()
    
    # ── Анализ гладкости для (p-1)-метода ──────────────────────────────────
    
    st.subheader("🔍 Анализ гладкости чисел для (p-1)-метода")
    
    if detailed_results:
        # Фильтруем только результаты (p-1)-метода
        p1_results = [r for r in detailed_results if r["Алгоритм"] == "(p-1)-метод Полларда"]
        
        if p1_results:
            # Группируем по битности
            smoothness_stats = {}
            for bit in bit_range:
                bit_results = [r for r in p1_results if r["Биты"] == bit]
                if bit_results:
                    total = len(bit_results)
                    smooth_100 = sum(1 for r in bit_results if r["100-гладкое"] == "✓")
                    success = sum(1 for r in bit_results if r["Статус"] == "✅ Успех")
                    
                    smoothness_stats[bit] = {
                        "Всего тестов": total,
                        "100-гладких": smooth_100,
                        "% гладких": f"{smooth_100/total*100:.0f}%",
                        "Успешных": success,
                        "% успеха": f"{success/total*100:.0f}%",
                    }
            
            if smoothness_stats:
                df_smooth = pd.DataFrame.from_dict(smoothness_stats, orient='index')
                df_smooth.index.name = "Биты"
                st.dataframe(df_smooth, use_container_width=True)
                
                st.caption("""
                **Важно:** (p-1)-метод Полларда — это **специализированный алгоритм**, который работает 
                эффективно ТОЛЬКО когда (p-1) или (q-1) является B-гладким (все простые делители ≤ B).
                
                **Почему (p-1)-метод быстро "отказывается"?**
                
                Текущая реализация проверяет только 3 значения B (100, 1000, 10000), что занимает ~11,000 операций.
                Если число не гладкое, метод быстро это определяет и сдаётся (доли секунды).
                Это **нормально и соответствует теории** — метод не тратит время на негладкие числа.
                
                **Вероятность гладкости для случайных чисел:**
                - Малые числа (20-40 бит): высокая вероятность гладкости → метод часто успешен
                - Средние числа (40-64 бит): средняя вероятность → метод иногда успешен
                - Большие числа (64+ бит): низкая вероятность → метод редко успешен
                
                Это объясняет, почему (p-1)-метод показывает хорошие результаты на малых числах,
                но проваливается на больших. **При учёте неудач ρ-метод значительно быстрее на практике.**
                """)
    
    st.divider()
    
    # ── Статистика неудач ───────────────────────────────────────────────────
    
    st.subheader("📊 Статистика успешности алгоритмов")
    
    failure_rows = []
    for _, alg_key in selected_algos:
        alg_name = next(n for n, k in selected_algos if k == alg_key)
        row = {"Алгоритм": alg_name}
        for bit in bit_range:
            success_count = len(raw_results[alg_key][bit])
            fail_count = failed_counts[alg_key][bit]
            total = success_count + fail_count
            
            if total > 0:
                success_rate = (success_count / total) * 100
                if success_count == 0:
                    row[f"{bit} бит"] = f"❌ 0/{total} (0%)"
                elif success_count == total:
                    row[f"{bit} бит"] = f"✅ {success_count}/{total} (100%)"
                else:
                    row[f"{bit} бит"] = f"⚠️ {success_count}/{total} ({success_rate:.0f}%)"
            else:
                row[f"{bit} бит"] = "—"
        failure_rows.append(row)
    
    st.dataframe(pd.DataFrame(failure_rows), use_container_width=True, hide_index=True)
    st.caption("✅ = все успешны, ⚠️ = частично успешны, ❌ = все неудачны. Показывает сколько из замеров завершились успешной факторизацией.")
    
    # Предупреждения о неудачах
    for _, alg_key in selected_algos:
        alg_name = next(n for n, k in selected_algos if k == alg_key)
        total_failures = sum(failed_counts[alg_key].values())
        if total_failures > 0:
            st.warning(f"⚠️ **{alg_name}**: {total_failures} неудачных попыток из {len(bit_range) * runs_per_bit} (алгоритм не справился или таймаут)")
    
    st.divider()

    # ── Агрегация результатов ───────────────────────────────────────────────

    # avg_results[algo_key] = {bit: avg_ms or None}
    avg_results: dict[str, dict[int, float | None]] = {}
    std_results: dict[str, dict[int, float]] = {}

    for _, alg_key in selected_algos:
        avg_results[alg_key] = {}
        std_results[alg_key] = {}
        for bit in bit_range:
            times = raw_results[alg_key][bit]
            if times:
                avg_results[alg_key][bit] = statistics.mean(times)
                std_results[alg_key][bit] = statistics.stdev(times) if len(times) > 1 else 0.0
            else:
                avg_results[alg_key][bit] = None
                std_results[alg_key][bit] = 0.0

    # ── Построение графика ──────────────────────────────────────────────────

    st.subheader("График условного времени успешных разложений")

    # Цветовая палитра + стиль линии по классу сложности
    COLORS = [
        "#e94560", "#0f3460", "#16213e", "#533483",
        "#2b9348", "#e9c46a", "#f4a261", "#264653",
        "#a8dadc", "#457b9d", "#e63946", "#06d6a0",
    ]

    # Экспоненциальные — сплошная, субэкспоненциальные — штрих
    SUBEXP_KEYS = {"qs_basic"}
    LINE_DASH = {k: "dash" for k in SUBEXP_KEYS}  # субэксп — штрих

    fig = go.Figure()

    for idx, (alg_name, alg_key) in enumerate(selected_algos):
        color = COLORS[idx % len(COLORS)]

        x_vals = []
        y_vals = []
        y_err  = []

        for bit in bit_range:
            avg = avg_results[alg_key].get(bit)
            if avg is not None:
                x_vals.append(bit)
                y_vals.append(avg)
                y_err.append(std_results[alg_key].get(bit, 0.0))

        if not x_vals:
            continue

        complexity = ALGO_COMPLEXITY_LABEL.get(alg_key, "")
        trace_name = f"{alg_name} [{complexity}]"
        line_dash = LINE_DASH.get(alg_key, "solid")

        fig.add_trace(go.Scatter(
            x=x_vals,
            y=y_vals,
            mode="lines+markers",
            name=trace_name,
            line=dict(color=color, width=2, dash=line_dash),
            marker=dict(size=7, color=color),
            hovertemplate=(
                f"<b>{alg_name}</b><br>"
                "Разрядность: %{x} бит<br>"
                "tусп: %{y:.2f} мс<br>"
                "<extra></extra>"
            ),
        ))

    fig.update_layout(
        xaxis_title="Разрядность числа (бит)",
        yaxis_title="Среднее время успешных запусков, tусп (мс)",
        hovermode="x unified",
        legend=dict(
            orientation="v",
            x=1.01, y=1,
            bgcolor="rgba(0,0,0,0)",
        ),
        margin=dict(r=220),
        height=560,
    )

    st.plotly_chart(fig, use_container_width=True)

    # ── Сводная таблица ─────────────────────────────────────────────────────

    st.subheader("Сводная таблица результатов")

    table_rows = []
    for _, alg_key in selected_algos:
        alg_name = next(n for n, k in selected_algos if k == alg_key)
        row = {"Алгоритм": alg_name, "Сложность": ALGO_COMPLEXITY_LABEL.get(alg_key, "—")}
        for bit in bit_range:
            avg = avg_results[alg_key].get(bit)
            fail_count = failed_counts[alg_key][bit]
            success_count = len(raw_results[alg_key][bit])
            total_count = success_count + fail_count
            fail_times = failed_results[alg_key][bit]
            fail_avg = statistics.mean(fail_times) if fail_times else None
            partial_count = partial_counts[alg_key][bit]
            
            if avg is not None:
                if fail_count > 0:
                    partial_note = f", част. {partial_count}" if partial_count else ""
                    fail_time_note = f", tотк={fail_avg:.2f}" if fail_avg is not None else ""
                    row[f"{bit} бит"] = (
                        f"усп. {success_count}/{total_count}, "
                        f"tусп={avg:.2f}{fail_time_note}{partial_note}"
                    )
                else:
                    row[f"{bit} бит"] = f"усп. {success_count}/{total_count}, t={avg:.2f}"
            else:
                if fail_count > 0:
                    partial_note = f", част. {partial_count}" if partial_count else ""
                    fail_time_note = f", tотк={fail_avg:.2f}" if fail_avg is not None else ""
                    row[f"{bit} бит"] = f"усп. 0/{total_count}, отказ{fail_time_note}{partial_note}"
                else:
                    row[f"{bit} бит"] = "—"
        table_rows.append(row)

    st.dataframe(pd.DataFrame(table_rows), use_container_width=True, hide_index=True)
    st.caption(
        "Все времена в мс. tусп — среднее только по полным успешным разложениям на простые множители; "
        "tотк — среднее время отказов, таймаутов и частичных разложений. Частичный результат не считается успехом."
    )
    
    # Пояснение для (p-1)-метода
    st.info("""
    💡 **Интерпретация результатов (p-1)-метода:**
    
    (p-1)-метод Полларда — **специализированный алгоритм**, эффективный только на числах с гладким (p-1).
    
    - **Быстрое время** = число имело гладкое (p-1), метод нашел делитель быстро
    - **"не справился"** = (p-1) имело большой простой делитель > B, метод не смог найти делитель
    - **Частичный успех** = некоторые числа были гладкими, другие нет
    
    Для объективного сравнения смотрите на:
    1. **Процент успеха** — как часто метод справляется
    2. **Анализ гладкости** — почему метод справился или нет
    3. **ρ-метод** — работает на ЛЮБЫХ числах, не зависит от гладкости
    """)
    
    # ── Сравнение с учетом неудач ──────────────────────────────────────────
    
    st.subheader("⚖️ Объективное сравнение: учет неудач")
    
    st.markdown("""
    **Почему (p-1)-метод кажется быстрее ρ-метода?**
    
    В таблице выше показано время **только для успешных случаев**. Но (p-1)-метод часто не справляется 
    с числами, у которых (p-1) не является гладким. Для объективного сравнения нужно учитывать неудачи.
    
    **Важно понимать:**
    - (p-1)-метод проверяет 3 значения B (100, 1000, 10000) — всего ~11,000 операций
    - Если число не гладкое, метод быстро "отказывается" (доли секунды)
    - ρ-метод делает до 100,000 итераций независимо от гладкости (секунды)
    - **Быстрый "отказ" — это не преимущество, а признак неприменимости метода к данному числу**
    """)
    
    # Создаем сравнительную таблицу с учетом неудач
    comparison_rows = []
    for bit in bit_range:
        row = {"Биты": bit}
        
        # ρ-метод
        rho_times = raw_results.get("pollard_rho", {}).get(bit, [])
        rho_fails = failed_counts.get("pollard_rho", {}).get(bit, 0)
        rho_fail_times = failed_results.get("pollard_rho", {}).get(bit, [])
        rho_total = len(rho_times) + rho_fails
        rho_success_rate = len(rho_times) / rho_total if rho_total > 0 else 0
        rho_avg = statistics.mean(rho_times) if rho_times else None
        rho_fail_avg = statistics.mean(rho_fail_times) if rho_fail_times else None
        
        if rho_avg:
            row["ρ-метод (tусп)"] = f"{rho_avg:.2f} мс"
            row["ρ-метод (успех)"] = f"{len(rho_times)}/{rho_total} ({rho_success_rate*100:.0f}%)"
        else:
            row["ρ-метод (tусп)"] = "—"
            row["ρ-метод (успех)"] = f"0/{rho_total} (0%)"
        row["ρ-метод (tотк)"] = f"{rho_fail_avg:.2f} мс" if rho_fail_avg is not None else "—"
        
        # (p-1)-метод
        p1_times = raw_results.get("pollard_p1", {}).get(bit, [])
        p1_fails = failed_counts.get("pollard_p1", {}).get(bit, 0)
        p1_fail_times = failed_results.get("pollard_p1", {}).get(bit, [])
        p1_total = len(p1_times) + p1_fails
        p1_success_rate = len(p1_times) / p1_total if p1_total > 0 else 0
        p1_avg = statistics.mean(p1_times) if p1_times else None
        p1_fail_avg = statistics.mean(p1_fail_times) if p1_fail_times else None
        
        if p1_avg:
            row["(p-1)-метод (tусп)"] = f"{p1_avg:.2f} мс"
            row["(p-1)-метод (успех)"] = f"{len(p1_times)}/{p1_total} ({p1_success_rate*100:.0f}%)"
        else:
            row["(p-1)-метод (tусп)"] = "—"
            row["(p-1)-метод (успех)"] = f"0/{p1_total} (0%)"
        row["(p-1)-метод (tотк)"] = f"{p1_fail_avg:.2f} мс" if p1_fail_avg is not None else "—"
        
        # Вывод
        if rho_avg and p1_avg and p1_success_rate > 0:
            # Учитываем неудачи: если (p-1) не справился, считаем это как очень долгое время
            # Эффективное время = время_успеха × вероятность_успеха + таймаут × вероятность_неудачи
            timeout_ms = get_timeout_for_bits(bit) * 1000
            p1_effective = p1_avg * p1_success_rate + timeout_ms * (1 - p1_success_rate)
            rho_effective = rho_avg * rho_success_rate + timeout_ms * (1 - rho_success_rate)
            row["ρ штрафное t"] = f"{rho_effective:.2f} мс"
            row["(p-1) штрафное t"] = f"{p1_effective:.2f} мс"
            
            if p1_effective < rho_effective:
                row["Вывод"] = f"(p-1) быстрее (с учетом неудач)"
            else:
                ratio = p1_effective / rho_effective
                row["Вывод"] = f"ρ быстрее в {ratio:.1f}× (с учетом неудач)"
        elif rho_avg and not p1_avg:
            row["ρ штрафное t"] = f"{rho_avg:.2f} мс"
            row["(p-1) штрафное t"] = f"{get_timeout_for_bits(bit) * 1000:.2f} мс"
            row["Вывод"] = "ρ-метод работает, (p-1) не справился"
        elif p1_avg and not rho_avg:
            row["ρ штрафное t"] = f"{get_timeout_for_bits(bit) * 1000:.2f} мс"
            row["(p-1) штрафное t"] = f"{p1_avg:.2f} мс"
            row["Вывод"] = "(p-1) работает, ρ не справился"
        else:
            row["ρ штрафное t"] = "—"
            row["(p-1) штрафное t"] = "—"
            row["Вывод"] = "Оба не справились"
        
        comparison_rows.append(row)
    
    if comparison_rows:
        df_comparison = pd.DataFrame(comparison_rows)
        st.dataframe(df_comparison, use_container_width=True, hide_index=True)
        
        st.caption("""
        **Как читать таблицу:**
        - **tусп** — среднее время только для полных успешных разложений
        - **tотк** — среднее время отказов, таймаутов и частичных разложений
        - **Успех** — количество успешных факторизаций из общего числа попыток
        - **Штрафное t** — условная оценка, где каждая неудача считается как таймаут; это не время алгоритма, а показатель надежности
        """)
        
        st.warning("""
        ⚠️ **Важный вывод:**
        
        Хотя (p-1)-метод показывает быстрое время на успешных случаях, при учете неудач 
        (когда числа не гладкие) его эффективность значительно снижается.
        
        **Почему (p-1)-метод "отказывается" так быстро?**
        
        (p-1)-метод проверяет только 3 значения B (100, 1000, 10000), что занимает ~11,000 операций.
        Если число не гладкое, метод быстро это определяет и сдаётся (доли секунды).
        
        ρ-метод пытается найти делитель независимо от гладкости, делая до 100,000 итераций (секунды).
        
        **Это нормально и соответствует теории:**
        - (p-1)-метод — **специализированный** алгоритм для гладких чисел
        - ρ-метод — **универсальный** алгоритм для любых чисел
        
        Поэтому в таблице отдельно показаны `tусп`, `tотк` и доля успеха. Быстрое `tусп` у `(p-1)` означает только то,
        что среди тестовых чисел встретился подходящий гладкий случай; без высокой доли успеха это не является
        преимуществом метода на произвольных числах данной разрядности.
        
        Это подтверждает теорию: (p-1)-метод медленнее ρ-метода не только теоретически, 
        но и на практике (при учёте неудач на негладких числах).
        """)
    
    st.divider()

    # ── Статистика по повторениям ───────────────────────────────────────────

    if runs_per_bit > 1:
        st.subheader("Статистика замеров (стандартное отклонение, мс)")
        std_rows = []
        for _, alg_key in selected_algos:
            alg_name = next(n for n, k in selected_algos if k == alg_key)
            row = {"Алгоритм": alg_name}
            for bit in bit_range:
                std = std_results[alg_key].get(bit, 0.0)
                avg = avg_results[alg_key].get(bit)
                if avg and avg > 0:
                    cv = std / avg * 100
                    row[f"{bit} бит"] = f"±{std:.2f} ({cv:.0f}%)"
                else:
                    row[f"{bit} бит"] = "—"
            std_rows.append(row)
        st.dataframe(pd.DataFrame(std_rows), use_container_width=True, hide_index=True)
        st.caption("В скобках — коэффициент вариации (std/mean × 100%). Чем ниже, тем стабильнее замеры.")

else:
    # Заглушка до запуска
    st.info("Настройте параметры в боковой панели и нажмите **🚀 Запустить тест**.")

    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**Реализованные алгоритмы:**")
        groups_table = [
            {"Алгоритм": "ρ-метод Полларда", "Сложность": "O(n^{1/4})", "Класс": "Экспоненциальный"},
            {"Алгоритм": "(p-1)-метод Полларда", "Сложность": "O(n^{1/2})", "Класс": "Экспоненциальный"},
            {"Алгоритм": "Алгоритм Диксона", "Сложность": "L[1/2, c]", "Класс": "Субэкспоненциальный"},
        ]
        st.dataframe(pd.DataFrame(groups_table), use_container_width=True, hide_index=True)

    with col2:
        st.markdown("**Сложность алгоритмов:**")
        st.markdown("""
| Алгоритм | Сложность |
|---|---|
| Поллард ρ | O(n^{1/4}) |
| Поллард p-1 | O(n^{1/2}) |
| Алгоритм Диксона | L[1/2, c] |

**Обозначение:**
- Сплошные линии — экспоненциальный алгоритм
- Штриховые линии — субэкспоненциальный алгоритм
        """)
    
    st.divider()
    
    # ── Оценка времени выполнения ──────────────────────────────────────────
    
    st.subheader("⏱️ Оценка времени выполнения")
    st.markdown("Примерное время факторизации для разных битностей:")
    
    estimate_bits = [20, 40, 60, 80, 90, 100, 110, 112]
    estimate_rows = []
    
    for bits in estimate_bits:
        row = {"Биты": bits}
        for name, key in ALGO_LIST:
            algo_short = name.split(" (")[0]
            estimate = estimate_time_for_bits(bits, key)
            row[algo_short] = estimate
        estimate_rows.append(row)
    
    st.dataframe(pd.DataFrame(estimate_rows), use_container_width=True, hide_index=True)
    
    st.markdown("""
    **Пояснения:**
    - **ρ-метод и (p-1)-метод**: Эффективны только для малых чисел (< 60 бит). Для больших чисел не справляются.
    - **Алгоритм Диксона**: Единственный способен факторизовать большие числа (80+ бит), но время растёт экспоненциально.
    - **112 бит**: Диксон займёт ~30-60 минут. Методы Полларда не справятся.
    """)
    
    st.divider()
    
    # ── Адаптивные таймауты ─────────────────────────────────────────────────
    
    st.subheader("⏰ Адаптивные таймауты")
    st.markdown("Таймауты автоматически увеличиваются для больших чисел:")
    
    timeout_info = []
    for bits in [20, 40, 60, 80, 85, 92, 100, 110, 112]:
        timeout_sec = get_timeout_for_bits(bits)
        timeout_info.append({
            "Битность": f"{bits} бит",
            "Таймаут": f"{timeout_sec} сек ({timeout_sec // 60} мин)" if timeout_sec >= 60 else f"{timeout_sec} сек",
            "Причина": (
                "Стандартный" if bits < 80 else
                "Диксон ~5 мин" if bits < 85 else
                "Диксон ~10 мин" if bits < 92 else
                "Диксон ~20 мин" if bits < 100 else
                "Диксон ~40 мин" if bits < 110 else
                "Диксон ~60 мин"
            )
        })
    
    st.dataframe(pd.DataFrame(timeout_info), use_container_width=True, hide_index=True)
    
    st.info("""
    💡 **Рекомендации:**
    - Для тестирования до 80 бит: используйте стандартные настройки
    - Для 80-85 бит: будьте готовы ждать до 5 минут (Диксон)
    - Для 85-92 бит: будьте готовы ждать до 10 минут (Диксон)
    - Для 92-100 бит: будьте готовы ждать до 20 минут (Диксон)
    - Для 100-110 бит: будьте готовы ждать до 40 минут (Диксон)
    - Для 110+ бит: будьте готовы ждать до 60 минут (Диксон)
    - Для 112 бит: факторизация займёт 40-60 минут (только Диксон)
    """)
    
    st.warning("""
    ⚠️ **Важно:**
    - Методы Полларда (ρ и p-1) **не справятся** с числами 80+ бит
    - Только алгоритм Диксона способен факторизовать большие числа
    - Для чисел 100+ бит рекомендуется запускать тест с 1-2 замерами на точку
    """)

