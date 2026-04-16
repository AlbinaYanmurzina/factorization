import streamlit as st
import requests
import random
import math
import platform
import statistics
import plotly.graph_objects as go
import pandas as pd
from sympy import nextprime

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

# ── Теоретические кривые сложности ─────────────────────────────────────────

def theoretical_curve(algo_key: str, bits_list: list) -> list | None:
    """
    Возвращает список относительных значений теоретической сложности
    (нормированных к первой точке), или None если кривая не определена.
    """
    curves = {
        # O(n^(1/4)) = O(2^(bits/4))
        "pollard_rho":  lambda b: 2 ** (b / 4),
        # O(n^(1/4)) аналогично
        "squfof":       lambda b: 2 ** (b / 4),
        # O(n^(1/2)) = O(2^(bits/2))
        "fermat":       lambda b: 2 ** (b / 2),
        "pollard_p1":   lambda b: 2 ** (b / 2),
        "williams_p1":  lambda b: 2 ** (b / 2),
        # L-нотация: exp(c * sqrt(bits * ln2 * ln(bits * ln2)))
        "cfrac":        lambda b: math.exp(math.sqrt(b * math.log(2) * math.log(b * math.log(2) + 1))),
        "qs_basic":     lambda b: math.exp(math.sqrt(b * math.log(2) * math.log(b * math.log(2) + 1))),
        "qs_optimized": lambda b: math.exp(math.sqrt(b * math.log(2) * math.log(b * math.log(2) + 1))),
        "qs_auto":      lambda b: math.exp(math.sqrt(b * math.log(2) * math.log(b * math.log(2) + 1))),
        "qs_lpv":       lambda b: math.exp(math.sqrt(b * math.log(2) * math.log(b * math.log(2) + 1))),
        "qs_mpqs":      lambda b: math.exp(math.sqrt(b * math.log(2) * math.log(b * math.log(2) + 1))),
        "qs_mpqs_parallel": lambda b: math.exp(math.sqrt(b * math.log(2) * math.log(b * math.log(2) + 1))),
    }
    fn = curves.get(algo_key)
    if fn is None:
        return None
    raw = [fn(b) for b in bits_list]
    if raw[0] == 0:
        return None
    # Нормируем: первая точка = первое реальное ненулевое значение
    return raw

ALGO_COMPLEXITY_LABEL = {
    "pollard_rho":      "O(n^{1/4})",
    "squfof":           "O(n^{1/4})",
    "fermat":           "O(n^{1/2})",
    "pollard_p1":       "O(n^{1/2})",
    "williams_p1":      "O(n^{1/2})",
    "cfrac":            "L[1/2, c]",
    "qs_basic":         "L[1/2, c]",
    "qs_optimized":     "L[1/2, c]",
    "qs_auto":          "L[1/2, c]",
    "qs_lpv":           "L[1/2, c]",
    "qs_mpqs":          "L[1/2, c]",
    "qs_mpqs_parallel": "L[1/2, c]",
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

ALGO_LIST = [
    ("ρ-метод Полларда (разд. 3.4)",                        "pollard_rho"),
    ("(p-1)-метод Полларда (разд. 3.2)",                    "pollard_p1"),
    ("Алгоритм Диксона (Basic QS, разд. 6.1)",              "qs_basic"),
]

# Алгоритмы, которые слишком медленны на больших числах
SLOW_ABOVE_BITS = {
    "qs_basic": 32,
}

TIMEOUT_MS = 30_000  # считаем тайм-аут если время > 30 с


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
max_bits = st.sidebar.slider("Максимальная разрядность (бит)", 16, 60, 40, step=2)
step_bits = st.sidebar.slider("Шаг разрядности (бит)", 2, 8, 4, step=2)
runs_per_bit = st.sidebar.slider("Замеров на точку графика", 1, 7, 3)

st.sidebar.markdown("---")
st.sidebar.subheader("Выбор алгоритмов")

selected_algos = []
for name, key in ALGO_LIST:
    if st.sidebar.checkbox(name, value=True, key=f"cb_{key}"):
        selected_algos.append((name, key))

st.sidebar.markdown("---")
st.sidebar.subheader("Отображение")
show_theory = st.sidebar.toggle("Теоретические кривые O(f(n))", value=True)

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

    progress_bar = st.progress(0, text="Запуск...")
    status_text = st.empty()

    for bit in bit_range:
        test_numbers = [generate_semiprime(bit) for _ in range(runs_per_bit)]

        for alg_name, alg_key in selected_algos:
            step_counter += 1
            progress_bar.progress(
                step_counter / total_steps,
                text=f"{alg_name} @ {bit} бит ({step_counter}/{total_steps})"
            )

            if alg_key in SLOW_ABOVE_BITS and bit > SLOW_ABOVE_BITS[alg_key]:
                # Пропускаем — слишком медленно
                continue

            for num in test_numbers:
                try:
                    res = requests.post(
                        "http://127.0.0.1:8453/api/factorize",
                        json={"number": str(num), "algorithm": alg_key},
                        timeout=35,
                    )
                    if res.status_code == 200:
                        t = res.json()["time_ms"]
                        # Тайм-аут на бэкенде возвращает ~30000 мс — помечаем как None
                        if t < TIMEOUT_MS:
                            raw_results[alg_key][bit].append(t)
                except Exception:
                    pass

    progress_bar.empty()
    status_text.empty()
    st.success(f"Эксперимент завершён. Протестировано разрядностей: {len(bit_range)}, алгоритмов: {len(selected_algos)}, замеров на точку: {runs_per_bit}.")

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

    st.subheader("График зависимости времени от разрядности")

    # Цветовая палитра + стиль линии по классу сложности
    COLORS = [
        "#e94560", "#0f3460", "#16213e", "#533483",
        "#2b9348", "#e9c46a", "#f4a261", "#264653",
        "#a8dadc", "#457b9d", "#e63946", "#06d6a0",
    ]

    # Экспоненциальные — сплошная, субэкспоненциальные — штрих
    SUBEXP_KEYS = {
        "cfrac", "qs_basic", "qs_optimized", "qs_auto",
        "qs_lpv", "qs_mpqs", "qs_mpqs_parallel",
    }
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
                "Время: %{y:.2f} мс<br>"
                "<extra></extra>"
            ),
        ))

        # Теоретическая кривая
        if show_theory:
            theory_raw = theoretical_curve(alg_key, x_vals)
            if theory_raw and y_vals:
                # Нормируем: масштабируем теорию к первой реальной точке
                scale = y_vals[0] / theory_raw[0] if theory_raw[0] != 0 else 1
                theory_scaled = [v * scale for v in theory_raw]

                fig.add_trace(go.Scatter(
                    x=x_vals,
                    y=theory_scaled,
                    mode="lines",
                    name=f"{alg_name} (теория)",
                    line=dict(color=color, width=1.5, dash="dot"),
                    opacity=0.5,
                    showlegend=False,
                    hovertemplate=(
                        f"<b>{alg_name} — теория {complexity}</b><br>"
                        "Разрядность: %{x} бит<br>"
                        "Норм. значение: %{y:.2f}<br>"
                        "<extra></extra>"
                    ),
                ))

    fig.update_layout(
        xaxis_title="Разрядность числа (бит)",
        yaxis_title="Среднее время (мс)",
        hovermode="x unified",
        legend=dict(
            orientation="v",
            x=1.01, y=1,
            bgcolor="rgba(0,0,0,0)",
        ),
        margin=dict(r=220),
        height=560,
    )

    # Аннотация: пунктир = теория, стиль линий = класс сложности
    if show_theory:
        fig.add_annotation(
            text="— — пунктир: теоретическая O(f(n))",
            xref="paper", yref="paper",
            x=0, y=-0.12,
            showarrow=False,
            font=dict(size=11, color="gray"),
        )
    fig.add_annotation(
        text="сплошная = экспоненциальные · штрих = субэкспоненциальные",
        xref="paper", yref="paper",
        x=0, y=-0.17,
        showarrow=False,
        font=dict(size=11, color="gray"),
    )

    st.plotly_chart(fig, use_container_width=True)

    # ── Сводная таблица ─────────────────────────────────────────────────────

    st.subheader("Сводная таблица (среднее время, мс)")

    table_rows = []
    for _, alg_key in selected_algos:
        alg_name = next(n for n, k in selected_algos if k == alg_key)
        row = {"Алгоритм": alg_name, "Сложность": ALGO_COMPLEXITY_LABEL.get(alg_key, "—")}
        for bit in bit_range:
            avg = avg_results[alg_key].get(bit)
            row[f"{bit} бит"] = f"{avg:.2f}" if avg is not None else "—"
        table_rows.append(row)

    st.dataframe(pd.DataFrame(table_rows), use_container_width=True, hide_index=True)

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
            {"Алгоритм": "Алгоритм Диксона (QS Basic)", "Сложность": "L[1/2, c]", "Класс": "Субэкспоненциальный"},
        ]
        st.dataframe(pd.DataFrame(groups_table), use_container_width=True, hide_index=True)

    with col2:
        st.markdown("**Теоретические кривые:**")
        st.markdown("""
| Алгоритм | Сложность |
|---|---|
| Поллард ρ | O(n^{1/4}) |
| Поллард p-1 | O(n^{1/2}) |
| QS Basic | L[1/2, c] |

**Как читать график:**
- Сплошные линии — экспоненциальные алгоритмы
- Штриховые линии — субэкспоненциальные
- Пунктирные линии — теоретические кривые O(f(n))
        """)

