import streamlit as st
import requests
import pandas as pd
import plotly.express as px

st.set_page_config(page_title="Пошаговая визуализация", page_icon="🔍")

st.title("Пошаговая визуализация алгоритма")

# ── Формулы для каждого типа шага ───────────────────────────────────────────

STEP_LATEX = {
    # ── ρ-метод Полларда ────────────────────────────────────────────────────
    "инициализац": r"f(x) = (x^2 + c) \bmod n, \quad x_0 = y_0 = 2",
    "итерац": (
        r"x_{i+1} = f(x_i),\quad y_{i+1} = f(f(y_i)),"
        r"\quad d_i = \gcd\!\bigl(|x_i - y_i|,\; n\bigr)"
    ),
    "вырожден": (
        r"d = n \;\Rightarrow\; x \equiv y \pmod{n}"
        r"\;\Rightarrow\;\text{цикл замкнулся, меняем } c"
    ),
    "найден нетривиальный": (
        r"1 < d = \gcd(|x - y|,\, n) < n"
        r"\;\Rightarrow\; d \mid n"
    ),
    "разбиение числа": r"n = d \times \frac{n}{d}, \quad d = \gcd(|x-y|,\,n)",
    "начало факторизации": (
        r"n = p_1^{a_1} \cdot p_2^{a_2} \cdots p_k^{a_k},"
        r"\quad p_i \text{ — простые}"
    ),
    "простое число найдено": r"p \text{ — простое} \;\Leftrightarrow\; \nexists\, d : 1 < d < p,\; d \mid p",
    "составное число": r"\exists\, d : 1 < d < n,\; d \mid n",

    # ── (p−1)-метод Полларда ────────────────────────────────────────────────
    "p−1": (
        r"a^{p-1} \equiv 1 \pmod{p} \;\Rightarrow\;"
        r"p \mid \gcd(a^M - 1,\; n),"
        r"\quad M = \prod_{q \le B} q^{\lfloor \log_q B \rfloor}"
    ),
    "p-1": (
        r"a^{p-1} \equiv 1 \pmod{p} \;\Rightarrow\;"
        r"p \mid \gcd(a^M - 1,\; n),"
        r"\quad M = \prod_{q \le B} q^{\lfloor \log_q B \rfloor}"
    ),
    "стратегия факторизации": (
        r"(p-1)\text{-гладкое: все простые } q \mid (p-1)"
        r"\text{ удовлетворяют } q \le B"
    ),
    "попытка с b": (
        r"M_B = \prod_{\substack{q \le B \\ q\text{ — простое}}} q^{\lfloor \log_q B \rfloor},"
        r"\quad a \leftarrow a^{q^k} \bmod n"
    ),
    "промежуточные шаги": (
        r"a \leftarrow a^{p^k} \bmod n,"
        r"\quad d = \gcd(a - 1,\; n)"
    ),
    "найден делитель": (
        r"(p-1) \text{ является } B\text{-гладким}"
        r"\;\Rightarrow\; p \mid \gcd(a^{M_B}-1,\,n)"
    ),
    "неудача при": (
        r"\gcd(a^{M_B}-1,\,n) \in \{1, n\}"
        r"\;\Rightarrow\; \exists\, q \mid (p-1) : q > B"
    ),

    # ── Квадратичное решето (Диксон) ────────────────────────────────────────
    "запуск": (
        r"x^2 \equiv y^2 \pmod{n}"
        r"\;\Rightarrow\; n \mid (x-y)(x+y)"
        r"\;\Rightarrow\; \gcd(x-y,\,n) \text{ — делитель}"
    ),
    "выбор параметра b": (
        r"B_{\mathrm{opt}} = \exp\!\Bigl({\tfrac{1}{2}}\sqrt{\ln n \cdot \ln\ln n}\Bigr),"
        r"\quad \pi(B) \approx \frac{B}{\ln B}"
    ),
    "факторная база": (
        r"\mathrm{FB} = \Bigl\{p \le B \;\Big|\; \left(\frac{n}{p}\right) = 1\Bigr\},"
        r"\quad \left(\frac{n}{p}\right) = n^{\frac{p-1}{2}} \bmod p"
    ),
    "гладк": (
        r"Q(x) = x^2 - n \equiv 0 \pmod{p_i}"
        r"\;\Leftrightarrow\; Q(x) = \prod_{p_i \in \mathrm{FB}} p_i^{e_i}"
    ),
    "построение матрицы": (
        r"A \in \mathrm{GF}(2)^{m \times k},"
        r"\quad a_{ij} = e_{ij} \bmod 2,"
        r"\quad Q(x_i) = \prod_j p_j^{e_{ij}}"
    ),
    "расширенная матрица": (
        r"\bigl[\,A \;\big|\; I_m\,\bigr] \in \mathrm{GF}(2)^{m \times (k+m)},"
        r"\quad I_m = \mathrm{diag}(1,\ldots,1)"
    ),
    "ход метода гаусса": (
        r"\mathbf{r}_j \leftarrow \mathbf{r}_j \oplus \mathbf{r}_{\mathrm{pivot}}"
        r"\text{ если } a_{j,c}=1;"
        r"\quad 0 \oplus 0 = 0,\; 1 \oplus 1 = 0,\; 0 \oplus 1 = 1"
    ),
    "результат гаусса": (
        r"\exists\, S \subseteq [m]:\; \sum_{i \in S} \mathbf{a}_i \equiv \mathbf{0} \pmod{2}"
        r"\;\Rightarrow\; \prod_{i \in S} Q(x_i) = Y^2"
    ),
    "зависимост": (
        r"X = \prod_{i \in S} x_i \bmod n,"
        r"\quad Y = \prod_j p_j^{\,\sum_{i\in S} e_{ij}/2} \bmod n,"
        r"\quad d = \gcd(X - Y,\, n)"
    ),
    "этап 4": (
        r"X^2 \equiv Y^2 \pmod{n},"
        r"\quad d_1 = \gcd(X-Y,\,n),"
        r"\quad d_2 = \gcd(X+Y,\,n)"
    ),
    "факторизация завершена": (
        r"n = p \times q, \quad p = \gcd(X - Y,\, n), \quad q = \frac{n}{p}"
    ),
    "провал": (
        r"\forall\, S:\; \gcd(X_S - Y_S,\,n) \in \{1,\,n\}"
        r"\;\Rightarrow\; \text{увеличить } B"
    ),
    "профилирован": (
        r"T = T_{\mathrm{FB}} + T_{\mathrm{sieve}} + T_{\mathrm{GF2}},"
        r"\quad T_{\mathrm{sieve}} = O\!\left(B \cdot \frac{\sqrt{n}}{B}\right)"
    ),
}

def get_step_latex(step_name: str) -> str | None:
    lower = step_name.lower()
    for key, formula in STEP_LATEX.items():
        if key in lower:
            return formula
    return None

# ── Тестовая база чисел ──────────────────────────────────────────────────────

TEST_NUMBERS = [
    # (number, bits, factors_hint)
    ("7904670551851",                                    45,  "2623991 × 3012461"),
    ("334384951600421767",                               60,  "314902243 × 1061869069"),
    ("9430176438789578115359",                           75,  "106294879963 × 88717127693"),
    ("292097409001938271188895309",                      90,  "10100880729463 × 28918013866843"),
    ("15879781158052216803368704409591",                 105, "3712658127928111 × 4277199949706681"),
]

# Числа специально подобранные для демонстрации (p-1)-метода Полларда.
# У каждого числа один из простых множителей p имеет (p-1) = B-гладкое число
# (все простые делители (p-1) малы), что гарантирует успех метода.
TEST_NUMBERS_P1 = [
    # (number, bits, factors_hint, max_prime_of_p_minus_1)
    ("838867",              20, "751 × 1117",                5),
    ("22223627",            25, "6301 × 3527",               7),
    ("1029256801",          30, "52489 × 19609",             3),
    ("1270723483397",       41, "1037233 × 1225109",         7),
    ("407009128958713",     49, "24000001 × 16958713",       5),
    ("1909115023580018711", 61, "560010907 × 3409067573",   17),
]

ALGO_OPTIONS = [
    "ρ-метод Полларда (разд. 3.4)",
    "(p-1)-метод Полларда (разд. 3.2)",
    "Алгоритм Диксона (разд. 6.1)",
]

ALGO_MAP = {
    "ρ-метод Полларда (разд. 3.4)":    "pollard_rho",
    "(p-1)-метод Полларда (разд. 3.2)": "pollard_p1",
    "Алгоритм Диксона (разд. 6.1)":    "qs_basic",
}

IS_QS = {"qs_basic"}

# ── Session state ────────────────────────────────────────────────────────────

if "history" not in st.session_state:
    st.session_state.history = []   # список dict с результатами запусков
if "last_steps" not in st.session_state:
    st.session_state.last_steps = None
if "last_meta" not in st.session_state:
    st.session_state.last_meta = None

# ── Рендер шага ─────────────────────────────────────────────────────────────

def render_step(idx: int, step_data: dict, is_qs_or_cfrac: bool):
    step_name = step_data.get("step", f"Шаг {idx + 1}")
    details = step_data.get("details", {})

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

    # ── Разделитель перед блоком СЛАУ ────────────────────────────────────────
    if "этап 3: построение матрицы" in name_lower:
        st.markdown("---")
        st.markdown(
            "<h4 style='margin-bottom:4px'>📐 ЭТАПЫ РЕШЕНИЯ СЛАУ над GF(2)</h4>"
            "<p style='color:gray;font-size:0.85em;margin-top:0'>Система линейных уравнений над полем из двух элементов — "
            "ключевой шаг алгоритма Диксона. Решается методом Гаусса с XOR-операциями.</p>",
            unsafe_allow_html=True,
        )

    with st.expander(f"{icon} {step_name}", expanded=(idx == 0)):
        formula = get_step_latex(step_name)
        if formula:
            st.latex(formula)

        # ── Дополнительные контекстные формулы по типу шага ─────────────────
        name_lower_r = step_name.lower()

        if "итерац" in name_lower_r and "алгоритма" in name_lower_r:
            st.latex(
                r"x_i \equiv y_i \pmod{p}"
                r"\;\Rightarrow\; p \mid (x_i - y_i)"
                r"\;\Rightarrow\; p \mid \gcd(x_i - y_i,\, n)"
            )

        if "построение матрицы" in name_lower_r:
            st.latex(
                r"Q(x_i) = \prod_{j=1}^{k} p_j^{\,e_{ij}}"
                r"\;\Longrightarrow\;"
                r"a_{ij} = e_{ij} \bmod 2 \in \{0, 1\}"
            )
            st.latex(
                r"\text{Цель: найти } \mathbf{v} \in \mathrm{GF}(2)^m \setminus \{\mathbf{0}\}"
                r"\text{ такой, что } \mathbf{v} \cdot A = \mathbf{0}"
            )

        if "расширенная матрица" in name_lower_r:
            st.latex(
                r"\begin{pmatrix} A & I_m \end{pmatrix}"
                r"\xrightarrow{\text{Гаусс над GF(2)}}"
                r"\begin{pmatrix} R & T \end{pmatrix},"
                r"\quad R\text{ — ступенчатый вид}"
            )

        if "ход метода гаусса" in name_lower_r:
            st.latex(
                r"\text{Для столбца } c: \quad "
                r"\mathbf{r}_j \leftarrow \mathbf{r}_j \oplus \mathbf{r}_p"
                r"\text{ при } a_{j,c} = 1,\; j \neq p"
            )
            st.latex(
                r"\mathrm{rank}(A) = r"
                r"\;\Rightarrow\; \dim\ker(A) = m - r"
                r"\;\Rightarrow\; m - r \text{ зависимостей}"
            )

        if "результат гаусса" in name_lower_r:
            st.latex(
                r"\sum_{i \in S} \mathbf{a}_i \equiv \mathbf{0} \pmod{2}"
                r"\;\Leftrightarrow\;"
                r"\prod_{i \in S} Q(x_i) = \prod_j p_j^{\,2k_j} = Y^2"
            )
            st.latex(
                r"\Bigl(\prod_{i \in S} x_i\Bigr)^2 \equiv Y^2 \pmod{n}"
                r"\;\Rightarrow\; n \mid (X - Y)(X + Y)"
            )

        if "зависимость #" in name_lower_r:
            st.latex(
                r"X = \prod_{i \in S} x_i \bmod n, \quad "
                r"Y = \prod_{j} p_j^{\,e_j / 2} \bmod n"
            )
            st.latex(
                r"d_1 = \gcd(X - Y,\, n), \quad d_2 = \gcd(X + Y,\, n)"
            )

        if "этап 1" in name_lower_r and "факторн" in name_lower_r:
            st.latex(
                r"\left(\frac{n}{p}\right) = n^{\frac{p-1}{2}} \bmod p"
                r"= \begin{cases} 1 & x^2 \equiv n \pmod{p}\text{ имеет решение}\\"
                r"-1 & \text{иначе} \end{cases}"
            )

        if "этап 2" in name_lower_r or ("поиск" in name_lower_r and "гладк" in name_lower_r):
            st.latex(
                r"Q(x) = x^2 - n, \quad x = \lfloor\sqrt{n}\rfloor + 1,\;"
                r"\lfloor\sqrt{n}\rfloor + 2, \ldots"
            )
            st.latex(
                r"Q(x) \text{ является } B\text{-гладким} \;\Leftrightarrow\;"
                r"Q(x) = \prod_{\substack{p \in \mathrm{FB}}} p^{e_p},\;"
                r"\text{остаток} = 1"
            )

        if "выбор параметра b" in name_lower_r:
            st.latex(
                r"L_n[\alpha, c] = \exp\!\bigl(c \cdot (\ln n)^\alpha (\ln\ln n)^{1-\alpha}\bigr)"
            )
            st.latex(
                r"B_{\mathrm{opt}} = L_n\!\left[\tfrac{1}{2},\, \tfrac{1}{2}\right]"
                r"= \exp\!\Bigl(\tfrac{1}{2}\sqrt{\ln n \cdot \ln\ln n}\Bigr)"
            )

        if "инициализац" in name_lower_r and "p" in name_lower_r and "метод" in name_lower_r:
            st.latex(
                r"M_B = \mathrm{lcm}(1, 2, \ldots, B)"
                r"= \prod_{\substack{q \le B \\ q\text{ — простое}}} q^{\lfloor \log_q B \rfloor}"
            )
            st.latex(
                r"a^{M_B} \bmod n: \quad "
                r"a \leftarrow a^{q^k} \bmod n \text{ для каждого } q^k \le B"
            )

        if "message" in details:
            st.code(details["message"], language=None)

        if "sieve_time_ms" in details:
            st.metric("⏱ Время поиска гладких чисел", f"{details['sieve_time_ms']:.3f} мс")

        if "matrix_data" in details and details["matrix_data"]:
            matrix = details["matrix_data"]

            # ── Система уравнений СЛАУ ───────────────────────────────────────
            if "equations_examples" in details and "factor_base_first" in details:
                fb = details["factor_base_first"]
                eqs = details["equations_examples"]
                n_cols = min(len(fb), 8)
                has_more_cols = len(fb) > n_cols

                # Разложения Q(xᵢ)
                decomp_latex = r"\begin{aligned}"
                for eq in eqs:
                    i, xi, qx = eq["i"], eq["x_i"], eq["q_x"]
                    fac = eq["factorization_latex"]
                    decomp_latex += rf"Q(x_{{{i}}}) &= {xi}^2 - n = {qx} = {fac} \\"
                decomp_latex += r"\end{aligned}"
                st.latex(decomp_latex)

                # Матричное уравнение — компоненты в колонках с подписями
                total_rows = len(matrix)
                rows_latex = []
                for eq in eqs:
                    coeffs = eq["row_mod2"][:n_cols]
                    row_str = " & ".join(map(str, coeffs))
                    if has_more_cols:
                        row_str += r" & \cdots"
                    rows_latex.append(row_str)
                if len(eqs) < total_rows:
                    rows_latex.append(r"\vdots & " * (n_cols - 1) + r"\vdots" + (r" & \ddots" if has_more_cols else ""))

                mat_body  = r" \\ ".join(rows_latex)
                v_body    = r" \\ ".join([rf"v_{{{eq['i']}}}" for eq in eqs])
                zero_body = r" \\ ".join(["0"] * len(eqs))
                if len(eqs) < total_rows:
                    v_body    += r" \\ \vdots"
                    zero_body += r" \\ \vdots"

                col_A, col_eq, col_v, col_eq2, col_0 = st.columns([4, 1, 1, 1, 1])
                with col_A:
                    with st.container(border=True):
                        st.caption("Матрица A — коэффициенты mod 2")
                        st.latex(r"\begin{pmatrix}" + mat_body + r"\end{pmatrix}")
                with col_eq:
                    st.write("")
                    st.write("")
                    st.latex(r"\cdot")
                with col_v:
                    with st.container(border=True):
                        st.caption("Вектор v — неизвестные")
                        st.latex(r"\begin{pmatrix}" + v_body + r"\end{pmatrix}")
                with col_eq2:
                    st.write("")
                    st.write("")
                    st.latex(r"\equiv")
                with col_0:
                    with st.container(border=True):
                        st.caption("Нулевой вектор")
                        st.latex(r"\begin{pmatrix}" + zero_body + r"\end{pmatrix}")

                st.caption("(mod 2) — сложение без переноса: 1+1=0, 0+1=1")

                # Уравнения в развёрнутом виде
                st.caption("Уравнения в развёрнутом виде (первые строки):")
                for eq in eqs[:3]:
                    i, xi = eq["i"], eq["x_i"]
                    coeffs = eq["row_mod2"]
                    nonzero_j = [j for j, c in enumerate(coeffs) if c == 1 and j < len(fb)]
                    if nonzero_j:
                        terms_eq = " + ".join([rf"v_{{{j+1}}}" for j in nonzero_j])
                        if len(coeffs) > len(fb):
                            terms_eq += r" + \cdots"
                    else:
                        terms_eq = "0"
                    st.latex(
                        rf"x_{{{i}}} = {xi}: \quad "
                        + terms_eq
                        + r" \equiv 0 \pmod{2}"
                    )

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

        if "table" in details and details["table"]:
            df = pd.DataFrame(details["table"])
            # PyArrow не поддерживает Python int произвольной точности —
            # конвертируем все целочисленные колонки в строки
            for col in df.columns:
                if df[col].dtype == object:
                    df[col] = df[col].apply(
                        lambda v: str(v) if isinstance(v, int) and (v > 2**63 - 1 or v < -(2**63)) else v
                    )
                elif df[col].dtype == "int64":
                    pass  # int64 PyArrow поддерживает
                else:
                    try:
                        df[col] = df[col].apply(
                            lambda v: str(v) if isinstance(v, int) and (v > 2**63 - 1 or v < -(2**63)) else v
                        )
                    except Exception:
                        df[col] = df[col].astype(str)
            st.dataframe(df, use_container_width=True, hide_index=True)

        if "FB" in details:
            fb = details["FB"]
            fb_size = details.get("FB_size", len(fb))
            pi_b = details.get("pi_B")
            b_val = details.get("B")

            # ── Заголовок с размером базы ────────────────────────────────────
            m1, m2, m3 = st.columns(3)
            m1.metric("|FB| — размер факторной базы", fb_size,
                help="Количество простых чисел p ≤ B, для которых символ Лежандра (n/p) = 1. "
                     "Именно столько столбцов в матрице СЛАУ.")
            m2.metric("Наименьший элемент", fb[0] if fb else "—",
                help="Первый элемент факторной базы — всегда 2.")
            m3.metric("Наибольший элемент", fb[-1] if fb else "—",
                help=f"Наибольшее простое число в базе — не превышает B = {b_val}.")

            if pi_b and b_val:
                st.latex(
                    rf"\mathrm{{FB}} = \{{p \le {b_val} \mid (n/p) = 1\}},"
                    rf"\quad |\mathrm{{FB}}| = {fb_size} < \pi({b_val}) = {pi_b}"
                )

            # ── Полная таблица факторной базы ────────────────────────────────
            with st.expander(f"📋 Все {fb_size} элементов факторной базы", expanded=False):
                # Разбиваем на страницы по 100 элементов
                page_size = 100
                n_pages = (fb_size + page_size - 1) // page_size

                if n_pages > 1:
                    page = st.selectbox(
                        "Страница:",
                        options=list(range(1, n_pages + 1)),
                        format_func=lambda p: f"Стр. {p}  ({(p-1)*page_size + 1}–{min(p*page_size, fb_size)})",
                        key=f"fb_page_{idx}",
                    )
                else:
                    page = 1

                chunk = fb[(page - 1) * page_size : page * page_size]

                # Отображаем в виде сетки: 10 чисел в строке
                cols_per_row = 10
                rows = [chunk[i:i+cols_per_row] for i in range(0, len(chunk), cols_per_row)]
                df_grid = pd.DataFrame(
                    rows,
                    columns=[f"p{(page-1)*page_size + i + 1}" for i in range(cols_per_row)],
                ).fillna("")
                # Убираем заголовки столбцов — они не нужны
                st.dataframe(
                    df_grid,
                    use_container_width=True,
                    hide_index=True,
                    column_config={c: st.column_config.TextColumn(c, width="small") for c in df_grid.columns},
                )
                st.caption(
                    f"Показаны элементы {(page-1)*page_size + 1}–{min(page*page_size, fb_size)} "
                    f"из {fb_size}. Все числа — простые p ≤ {b_val}, для которых x² ≡ n (mod p) разрешимо."
                )

        if "pi_B" in details and "FB" not in details:
            pi_b = details["pi_B"]
            b_val = details.get("B")
            st.info(f"π({b_val}) = {pi_b} — столько простых чисел ≤ {b_val} (размер факторной базы)")

        if "pi_B_data" in details:
            d = details["pi_B_data"]
            B      = d["B"]
            B_auto = d["B_auto"]
            pi_exact  = d["pi_B_exact"]
            pi_approx = d["pi_B_approx"]
            ln_B   = d["ln_B"]

            st.markdown("**Функция π(x) — количество простых чисел, не превышающих x:**")

            # Теорема о простых числах
            st.latex(
                r"\pi(x) \sim \frac{x}{\ln x} \quad (x \to \infty)"
                r"\qquad \text{(теорема о простых числах)}"
            )

            # Подстановка конкретного B
            st.latex(
                rf"\pi({B}) \approx \frac{{{B}}}{{\ln {B}}}"
                rf"= \frac{{{B}}}{{{ln_B}}}"
                rf"\approx {pi_approx}"
            )

            # Точное значение
            st.latex(
                rf"\pi({B}) = {pi_exact} \quad \text{{(точное значение, решето Эратосфена)}}"
            )

            # Оптимальный B из L-нотации
            st.markdown("**Оптимальная граница B по L-нотации:**")
            st.latex(
                r"B_{\mathrm{opt}} = L_n\!\left[\tfrac{1}{2},\,\tfrac{1}{2}\right]"
                r"= \exp\!\Bigl(\tfrac{1}{2}\sqrt{\ln n \cdot \ln\ln n}\Bigr)"
                rf"\approx {B_auto}"
            )

            # Итоговый размер базы
            col1, col2, col3 = st.columns(3)
            col1.metric("Граница B", B,
                help=(
                    "Граница гладкости — максимальный простой делитель, "
                    "который мы допускаем в разложении Q(x) = x² − n. "
                    "Чем больше B, тем больше чисел считаются гладкими, "
                    "но тем больше размер матрицы СЛАУ."
                ))
            col2.metric("π(B) точное", pi_exact,
                help=(
                    f"Точное количество простых чисел ≤ {B}, "
                    f"вычисленное решетом Эратосфена. "
                    f"Именно столько столбцов будет в матрице A над GF(2). "
                    f"Чтобы система имела нетривиальное решение, нужно найти "
                    f"хотя бы π(B)+1 гладких чисел."
                ))
            col3.metric("π(B) по ТПЧ", f"≈ {pi_approx}",
                help=(
                    f"Приближение по теореме о простых числах: π(B) ≈ B / ln(B). "
                    f"При B = {B}: {B} / ln({B}) = {B} / {ln_B} ≈ {pi_approx}. "
                    f"Теорема даёт асимптотическую оценку — точность растёт с ростом B."
                ))

        if "Примеры полиномов" in details and details["Примеры полиномов"]:
            st.caption("Примеры полиномов:")
            st.dataframe(pd.DataFrame(details["Примеры полиномов"]), use_container_width=True, hide_index=True)

        # Профилирование времени (круговая диаграмма)
        if "profiling" in details:
            prof = details["profiling"]
            labels = list(prof.keys())
            values = list(prof.values())
            
            fig = px.pie(
                names=labels,
                values=values,
                title="Распределение времени выполнения алгоритма Диксона",
                hole=0.3,  # donut chart
                color_discrete_sequence=px.colors.qualitative.Set3,
            )
            fig.update_traces(
                textposition='inside',
                textinfo='label+percent',
                hovertemplate='<b>%{label}</b><br>%{value:.3f} мс<br>%{percent}<extra></extra>',
            )
            fig.update_layout(
                showlegend=True,
                height=450,
                margin=dict(l=20, r=20, t=60, b=20),
            )
            st.plotly_chart(fig, use_container_width=True)
            
            # Текстовая таблица для точных значений
            prof_df = pd.DataFrame({
                "Этап": labels,
                "Время, мс": values,
                "Доля, %": [f"{v/sum(values)*100:.1f}%" for v in values]
            })
            st.dataframe(prof_df, use_container_width=True, hide_index=True)


# ── UI: генератор числа по битности ─────────────────────────────────────────

with st.container(border=True):
    st.caption("🎲 Генератор составного числа")
    gen_col_slider, gen_col_btn = st.columns([4, 1])
    with gen_col_slider:
        gen_bits = st.slider(
            "Битность числа:",
            min_value=10,
            max_value=200,
            value=16,
            step=1,
            key="gen_bits_slider",
            help="Генерирует число вида p × q, где p и q — случайные простые заданной битности",
        )
    with gen_col_btn:
        st.write("")
        gen_clicked = st.button("🎲 Сгенерировать", use_container_width=True, key="gen_btn")

    if gen_clicked:
        import random as _random
        from backend.algorithms.math_utils import is_prime as _is_prime

        def _gen_prime(bits: int) -> int:
            """Генерирует случайное простое число заданной битности."""
            while True:
                # Старший бит всегда 1, чтобы гарантировать битность
                n = _random.getrandbits(bits) | (1 << (bits - 1)) | 1
                if _is_prime(n, k=10):
                    return n

        half = gen_bits // 2
        # p и q примерно одинаковой битности, чтобы число было составным и интересным
        p_bits = half
        q_bits = gen_bits - half
        p = _gen_prime(p_bits)
        q = _gen_prime(q_bits)
        # Гарантируем p ≠ q
        attempts = 0
        while q == p and attempts < 20:
            q = _gen_prime(q_bits)
            attempts += 1
        generated = p * q
        st.session_state["generated_number"] = str(generated)
        st.session_state["generated_info"] = (p, q, gen_bits)

    if "generated_info" in st.session_state:
        p_g, q_g, bits_g = st.session_state["generated_info"]
        num_g = st.session_state["generated_number"]
        actual_bits = len(bin(int(num_g))) - 2
        st.info(
            f"n = {num_g}  |  {actual_bits} бит  |  {p_g} × {q_g}  "
            f"— нажмите «Использовать» чтобы подставить в поле ввода"
        )
        if st.button("✅ Использовать это число", key="use_generated_btn"):
            st.session_state["manual_input"] = st.session_state["generated_number"]
            st.rerun()

# ── UI: ввод числа ───────────────────────────────────────────────────────────

# Инициализируем значение по умолчанию если виджет ещё не создан
if "manual_input" not in st.session_state:
    st.session_state["manual_input"] = "8051"

col_input, col_algo, col_btn = st.columns([3, 2, 1])
with col_input:
    number_input = st.text_input("Составное число (n):", key="manual_input")
with col_algo:
    algo_choice = st.selectbox("Алгоритм:", ALGO_OPTIONS, key="algo_select")
with col_btn:
    st.write("")  # выравнивание по высоте
    st.write("")
    run_manual = st.button("▶ Запустить", type="primary", use_container_width=True)

# ── UI: параметр B (только для Диксона) ─────────────────────────────────────

b_override = None  # None = авто

if algo_choice == "Алгоритм Диксона (разд. 6.1)":
    with st.container(border=True):
        st.caption("⚙️ Параметр B — граница факторной базы (только для алгоритма Диксона)")
        b_mode = st.radio(
            "Режим выбора B:",
            ["Автоматический расчёт (L-нотация)", "Ручной ввод"],
            horizontal=True,
            key="b_mode_radio",
        )
        if b_mode == "Ручной ввод":
            # Вычисляем авто-значение для подсказки
            try:
                n_preview = int(number_input)
                import math as _math
                if n_preview > 1:
                    b_auto_preview = int(_math.exp(0.5 * _math.sqrt(_math.log(n_preview) * _math.log(_math.log(n_preview)))))
                    b_auto_preview = max(b_auto_preview, 100)
                    b_auto_preview = min(b_auto_preview, 50000)
                else:
                    b_auto_preview = 200
            except Exception:
                b_auto_preview = 200

            col_sl, col_hint = st.columns([3, 2])
            with col_sl:
                b_override = st.slider(
                    "Граница B:",
                    min_value=10,
                    max_value=50000,
                    value=b_auto_preview,
                    step=50,
                    key="b_slider",
                    help="Авто-значение выделено. Двигайте влево — меньше гладких чисел, вправо — тяжелее Гаусс.",
                )
            with col_hint:
                if b_override < b_auto_preview // 2:
                    st.warning(f"B = {b_override} сильно меньше оптимального ({b_auto_preview}).\nВероятно, гладких чисел не хватит.")
                elif b_override > b_auto_preview * 3:
                    st.warning(f"B = {b_override} сильно больше оптимального ({b_auto_preview}).\nМатрица Гаусса будет очень большой.")
                else:
                    st.info(f"Оптимальное B ≈ {b_auto_preview}\nВыбрано: B = {b_override}")


# ── UI: таблица тестовых чисел ───────────────────────────────────────────────

with st.expander("📋 Тестовая база криптографических чисел", expanded=False):
    st.caption("Нажмите «▶» для запуска одного алгоритма или «▶▶▶» для запуска всех трёх сразу.")

    db_cols = st.columns([4, 1, 2, 1, 1, 1, 1])
    headers = ["n", "Бит", "Множители", "ρ-Полларда", "(p-1)-Полларда", "Диксон", "Все"]
    for col, h in zip(db_cols, headers):
        col.markdown(f"**{h}**")

    for entry in TEST_NUMBERS:
        number, bits, factors_hint = entry
        row = st.columns([4, 1, 2, 1, 1, 1, 1])
        row[0].code(number, language=None)
        row[1].write(str(bits))
        row[2].write(factors_hint)
        if row[3].button("▶", key=f"rho_{number}", help="ρ-метод Полларда"):
            st.session_state["queued_number"] = number
            st.session_state["queued_algo"] = ALGO_OPTIONS[0]
            st.rerun()
        if row[4].button("▶", key=f"p1_{number}", help="(p-1)-метод Полларда"):
            st.session_state["queued_number"] = number
            st.session_state["queued_algo"] = ALGO_OPTIONS[1]
            st.rerun()
        if row[5].button("▶", key=f"qs_{number}", help="Алгоритм Диксона"):
            st.session_state["queued_number"] = number
            st.session_state["queued_algo"] = ALGO_OPTIONS[2]
            st.rerun()
        if row[6].button("▶▶▶", key=f"all_{number}", help="Запустить все три алгоритма"):
            st.session_state["queued_all_number"] = number
            st.rerun()

    # ── Специальные числа для (p-1)-метода ──────────────────────────────────
    st.divider()
    st.markdown(
        "**🎯 Числа специально подобранные для демонстрации (p-1)-метода Полларда**"
    )
    st.caption(
        "У каждого числа один из простых множителей p имеет (p−1) = B-гладкое число "
        "(все простые делители (p−1) не превышают B=100). "
        "Это гарантирует что (p−1)-метод найдёт делитель, тогда как для случайных чисел это не выполняется."
    )

    p1_cols = st.columns([4, 1, 2, 1, 1])
    p1_headers = ["n", "Бит", "Множители", "Наиб. дел. (p−1)", "(p-1)-Полларда"]
    for col, h in zip(p1_cols, p1_headers):
        col.markdown(f"**{h}**")

    for entry in TEST_NUMBERS_P1:
        number, bits, factors_hint, max_p = entry
        row = st.columns([4, 1, 2, 1, 1])
        row[0].code(number, language=None)
        row[1].write(str(bits))
        row[2].write(factors_hint)
        row[3].markdown(f"`≤ {max_p}`")
        if row[4].button("▶", key=f"p1_special_{number}", help="(p-1)-метод Полларда"):
            st.session_state["queued_number"] = number
            st.session_state["queued_algo"] = ALGO_OPTIONS[1]
            st.rerun()

# Подхватываем число из очереди (нажатие кнопки в таблице)
queued = "queued_number" in st.session_state
queued_all = "queued_all_number" in st.session_state

if queued:
    run_number     = st.session_state.pop("queued_number")
    run_algo       = st.session_state.pop("queued_algo")
    run_b_override = None  # из таблицы всегда авто
elif run_manual:
    run_number     = number_input
    run_algo       = algo_choice
    run_b_override = b_override  # None если авто, int если ручной
else:
    run_number     = None
    run_algo       = None
    run_b_override = None

# ── Выполнение всех алгоритмов сразу ────────────────────────────────────────

if queued_all:
    all_number = st.session_state.pop("queued_all_number")
    all_results = []

    with st.status(f"Запуск всех алгоритмов для n = {all_number}...", expanded=True) as all_status:
        for algo_label in ALGO_OPTIONS:
            algo_key = ALGO_MAP[algo_label]
            st.write(f"▶ {algo_label.split(' (')[0]}...")
            try:
                payload = {"number": all_number, "algorithm": algo_key}
                resp = requests.post("http://127.0.0.1:8453/api/factorize", json=payload)
                if resp.status_code == 200:
                    d = resp.json()
                    factors_str = " × ".join(d["factors"])
                    all_results.append({
                        "Алгоритм": algo_label.split(" (")[0],
                        "Результат": factors_str,
                        "Время, мс": round(d["time_ms"], 3),
                        "Шагов": len(d["steps"]),
                    })
                    # Сохраняем в историю
                    st.session_state.history.append({
                        "n": all_number,
                        "Бит": len(bin(int(all_number))) - 2,
                        "Алгоритм": algo_label.split(" (")[0],
                        "B": "авто",
                        "Результат": factors_str,
                        "Время, мс": round(d["time_ms"], 3),
                        "Шагов": len(d["steps"]),
                    })
                    # Последний алгоритм — сохраняем шаги для пошагового разбора
                    st.session_state.last_steps = d["steps"]
                    st.session_state.last_meta = {
                        "number": all_number,
                        "algo": algo_label,
                        "is_qs": algo_key in IS_QS,
                    }
                else:
                    all_results.append({
                        "Алгоритм": algo_label.split(" (")[0],
                        "Результат": "Ошибка",
                        "Время, мс": "-",
                        "Шагов": "-",
                    })
            except requests.exceptions.ConnectionError:
                all_status.update(label="Нет соединения с сервером", state="error", expanded=True)
                st.error("Не удалось подключиться к серверу.")
                st.stop()

        all_status.update(label=f"Все алгоритмы выполнены для n = {all_number}", state="complete", expanded=False)

    st.subheader(f"Сравнение алгоритмов  |  n = {all_number}")
    df_all = pd.DataFrame(all_results)
    st.dataframe(df_all, use_container_width=True, hide_index=True)

    # Мини-диаграмма сравнения времени
    numeric_rows = [r for r in all_results if isinstance(r["Время, мс"], (int, float))]
    if len(numeric_rows) > 1:
        fig_cmp = px.bar(
            pd.DataFrame(numeric_rows),
            x="Алгоритм", y="Время, мс",
            color="Алгоритм",
            text="Время, мс",
            title="Время выполнения по алгоритмам",
            color_discrete_sequence=px.colors.qualitative.Set2,
        )
        fig_cmp.update_traces(texttemplate="%{text:.3f} мс", textposition="outside")
        fig_cmp.update_layout(showlegend=False, height=320, margin=dict(t=50, b=20))
        st.plotly_chart(fig_cmp, use_container_width=True)

# ── Выполнение запроса ───────────────────────────────────────────────────────

if run_number and run_algo:
    algo_key = ALGO_MAP[run_algo]
    is_qs = algo_key in IS_QS

    # b_override актуален только для Диксона и только если задан вручную
    effective_b = run_b_override if (algo_key == "qs_basic" and run_b_override is not None) else None

    with st.status(f"Выполнение: {run_algo}...", expanded=True) as status:
        st.write("Отправка запроса на сервер...")
        try:
            payload = {"number": run_number, "algorithm": algo_key}
            if effective_b is not None:
                payload["b_override"] = effective_b
            response = requests.post("http://127.0.0.1:8453/api/factorize", json=payload)
            if response.status_code == 200:
                data = response.json()
                st.write(f"Получен ответ. Шагов: {len(data['steps'])}")
                status.update(label=f"Готово за {data['time_ms']:.3f} мс", state="complete", expanded=False)
            else:
                status.update(label="Ошибка сервера", state="error", expanded=True)
                st.error(response.json().get("detail", "Неизвестная ошибка"))
                st.stop()
        except requests.exceptions.ConnectionError:
            status.update(label="Нет соединения с сервером", state="error", expanded=True)
            st.error("Не удалось подключиться к серверу. Убедитесь, что FastAPI (backend) запущен.")
            st.stop()

    factors_str = " × ".join(data["factors"])
    st.success(f"**{run_number} = {factors_str}**  |  время: {data['time_ms']:.3f} мс")

    # Сохраняем в историю
    b_label = str(effective_b) if effective_b is not None else "авто"
    st.session_state.history.append({
        "n": run_number,
        "Бит": len(bin(int(run_number))) - 2,
        "Алгоритм": run_algo.split(" (")[0],
        "B": b_label,
        "Результат": factors_str,
        "Время, мс": round(data["time_ms"], 3),
        "Шагов": len(data["steps"]),
    })
    st.session_state.last_steps = data["steps"]
    st.session_state.last_meta  = {"number": run_number, "algo": run_algo, "is_qs": is_qs}

# ── История запусков ─────────────────────────────────────────────────────────

if st.session_state.history:
    st.divider()
    st.subheader("📊 История запусков")

    hist_col, clear_col = st.columns([5, 1])
    with clear_col:
        if st.button("🗑 Очистить", use_container_width=True):
            st.session_state.history = []
            st.session_state.last_steps = None
            st.session_state.last_meta  = None
            st.rerun()

    df_hist = pd.DataFrame(st.session_state.history)
    st.dataframe(df_hist, use_container_width=True, hide_index=True)

# ── Пошаговый разбор последнего запуска ─────────────────────────────────────

if st.session_state.last_steps:
    meta = st.session_state.last_meta
    st.divider()
    st.subheader(f"Шаги работы: {meta['algo'].split(' (')[0]}  |  n = {meta['number']}")

    for idx, step_data in enumerate(st.session_state.last_steps):
        render_step(idx, step_data, meta["is_qs"])
