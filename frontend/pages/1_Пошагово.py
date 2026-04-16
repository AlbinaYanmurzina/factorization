import streamlit as st
import requests
import pandas as pd
import plotly.express as px

st.set_page_config(page_title="Пошаговая визуализация", page_icon="🔍")

st.title("Пошаговая визуализация алгоритма")

# Формулы для каждого типа шага
STEP_LATEX = {
    "факторная база": r"\left(\frac{n}{p}\right) = 1",
    "гладк": r"Q(x) = x^2 - n",
    "гаусс": r"Ax = 0 \pmod{2}",
    "зависимост": r"\gcd(X - Y,\; n)",
    "cfrac": r"A_k^2 \equiv (-1)^k \cdot d_k \pmod{n}",
    "цепная дробь": r"A_k^2 \equiv Q_k \pmod{n}",
    "просеивани": r"\text{sieve}[i] \mathrel{+}= \log_2 p \quad \text{если } p \mid Q(x_i)",
    "полином": r"Q(x) = ax^2 + 2bx + c, \quad a = t^2",
    "инициализац": r"f(x) = (x^2 + c) \bmod n",
    "итерац": r"\gcd(|x - y|,\; n)",
    "лукас": r"V_k = P \cdot V_{k-1} - V_{k-2} \pmod{n}",
    "вильямс": r"\gcd(V_M - 2,\; n)",
    "p+1": r"\gcd(V_M - 2,\; n)",
    "p−1": r"a^M \equiv 1 \pmod{p} \;\Rightarrow\; p \mid \gcd(a^M - 1,\; n)",
    "p-1": r"a^M \equiv 1 \pmod{p} \;\Rightarrow\; p \mid \gcd(a^M - 1,\; n)",
    "золотая": r"Q_i = s^2 \;\Rightarrow\; \gcd(Q_{\text{prev}},\; n)",
    "squfof": r"D = kn,\quad q = \left\lfloor\frac{\lfloor\sqrt{D}\rfloor + P}{Q}\right\rfloor",
    "квадратичн": r"n = x^2 - y^2 = (x-y)(x+y)",
    "ферма": r"w = x^2 - n = y^2",
    "представлен": r"n = x^2 - y^2,\quad p = x - y,\quad q = x + y",
    "l-нотац": r"L(n,\tfrac{1}{2}) = \exp\!\left(\sqrt{\ln n \cdot \ln\ln n}\right)",
    "параметр b": r"B = \exp\!\left(0.5\sqrt{\ln n \cdot \ln\ln n}\right)",
    "выбор параметра": r"B = \exp\!\left(0.5\sqrt{\ln n \cdot \ln\ln n}\right)",
    "расчёт параметр": r"B = L(n,\tfrac{1}{2})^{\alpha}",
}

def get_step_latex(step_name: str) -> str | None:
    lower = step_name.lower()
    for key, formula in STEP_LATEX.items():
        if key in lower:
            return formula
    return None

def render_step(idx: int, step_data: dict, is_qs_or_cfrac: bool):
    step_name = step_data.get("step", f"Шаг {idx + 1}")
    details = step_data.get("details", {})

    # Определяем иконку по содержимому шага
    name_lower = step_name.lower()
    if "завершена" in name_lower or "успех" in name_lower:
        icon = "✅"
    elif "провал" in name_lower or "ошибка" in name_lower:
        icon = "❌"
    elif "этап" in name_lower:
        icon = "🔷"
    elif "зависимост" in name_lower:
        icon = "🔗"
    elif "результат" in name_lower:
        icon = "📊"
    else:
        icon = "🔹"

    with st.expander(f"{icon} {step_name}", expanded=(idx == 0)):
        # Формула рядом с шагом
        formula = get_step_latex(step_name)
        if formula:
            st.latex(formula)

        if "message" in details:
            st.code(details["message"], language=None)

        # Heatmap матрицы GF(2)
        if "matrix_data" in details and details["matrix_data"]:
            matrix = details["matrix_data"]
            st.caption(f"Матрица GF(2): {len(matrix)} × {len(matrix[0])} (показаны первые 30×30)")
            fig = px.imshow(
                matrix,
                color_continuous_scale=[[0, "#1a1a2e"], [1, "#e94560"]],
                title="Матрица показателей степеней mod 2",
                labels={"x": "Простые из FB", "y": "Гладкие числа"},
                aspect="auto",
            )
            fig.update_layout(
                coloraxis_showscale=False,
                margin=dict(l=0, r=0, t=40, b=0),
                height=max(200, min(len(matrix) * 12, 400)),
            )
            fig.update_traces(hovertemplate="строка %{y}, столбец %{x}: %{z}<extra></extra>")
            st.plotly_chart(fig, use_container_width=True)

        # Интерактивная таблица
        if "table" in details and details["table"]:
            df = pd.DataFrame(details["table"])
            st.dataframe(df, use_container_width=True, hide_index=True)

        # Факторная база
        if "FB" in details:
            fb = details["FB"]
            st.caption(f"Факторная база ({len(fb)} элементов): {fb[:20]}{'...' if len(fb) > 20 else ''}")

        # Примеры полиномов MPQS
        if "Примеры полиномов" in details and details["Примеры полиномов"]:
            st.caption("Примеры полиномов:")
            st.dataframe(pd.DataFrame(details["Примеры полиномов"]), use_container_width=True, hide_index=True)


# ── UI ──────────────────────────────────────────────────────────────────────

with st.container():
    col1, col2 = st.columns([2, 1])
    with col1:
        number_input = st.text_input("Введите составное число (n):", value="8051")
    with col2:
        algo_choice = st.selectbox(
            "Выберите алгоритм:",
            [
                "ρ-метод Полларда (разд. 3.4)",
                "(p-1)-метод Полларда (разд. 3.2)",
                "Алгоритм Диксона (Basic QS, разд. 6.1)",
            ]
        )

# Словарь для связи названий с ключами API бэкенда
algo_map = {
    "ρ-метод Полларда (разд. 3.4)": "pollard_rho",
    "(p-1)-метод Полларда (разд. 3.2)": "pollard_p1",
    "Алгоритм Диксона (Basic QS, разд. 6.1)": "qs_basic",
}

IS_QS_OR_CFRAC = {
    "qs_basic",
}

if st.button("Факторизовать", type="primary"):
    algo_key = algo_map[algo_choice]
    is_qs = algo_key in IS_QS_OR_CFRAC

    with st.status(f"Выполнение: {algo_choice}...", expanded=True) as status:
        st.write("Отправка запроса на сервер...")
        try:
            response = requests.post(
                "http://127.0.0.1:8453/api/factorize",
                json={"number": number_input, "algorithm": algo_key}
            )

            if response.status_code == 200:
                data = response.json()
                st.write(f"Получен ответ. Шагов: {len(data['steps'])}")
                status.update(
                    label=f"Готово за {data['time_ms']:.3f} мс",
                    state="complete",
                    expanded=False,
                )
            else:
                status.update(label="Ошибка сервера", state="error", expanded=True)
                st.error(response.json().get("detail", "Неизвестная ошибка"))
                st.stop()

        except requests.exceptions.ConnectionError:
            status.update(label="Нет соединения с сервером", state="error", expanded=True)
            st.error("Не удалось подключиться к серверу. Убедитесь, что FastAPI (backend) запущен.")
            st.stop()

    # Результат
    factors_str = " × ".join(data["factors"])
    st.success(f"**{number_input} = {factors_str}**  |  время: {data['time_ms']:.3f} мс")

    st.divider()
    st.subheader("Шаги работы алгоритма")

    for idx, step_data in enumerate(data["steps"]):
        render_step(idx, step_data, is_qs)
