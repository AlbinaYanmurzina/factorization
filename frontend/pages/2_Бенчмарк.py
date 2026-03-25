import streamlit as st
import requests
import time
import random
import plotly.graph_objects as go
import pandas as pd
from sympy import isprime, nextprime # Понадобится для генерации тестов

st.set_page_config(page_title="Сравнение алгоритмов", page_icon="📊", layout="wide")

st.title("📊 Вычислительные эксперименты и бенчмарки")

def generate_semiprime(bits: int) -> int:
    """Генерирует число N = p * q заданной разрядности в битах"""
    # Ищем два простых числа примерно одинаковой длины (bits // 2)
    start_p = random.getrandbits(bits // 2)
    start_q = random.getrandbits(bits // 2)
    
    p = nextprime(start_p)
    q = nextprime(start_q)
    
    return p * q

# Боковая панель настроек
st.sidebar.header("Настройки эксперимента")
min_bits = st.sidebar.slider("Минимальная разрядность (бит)", min_value=16, max_value=60, value=20, step=2)
max_bits = st.sidebar.slider("Максимальная разрядность (бит)", min_value=16, max_value=60, value=40, step=2)
runs_per_bit = st.sidebar.number_input("Кол-во тестов на одну разрядность", value=2, min_value=1)

if st.sidebar.button("🚀 Запустить бенчмарк", type="primary"):
    
    bit_range = list(range(min_bits, max_bits + 1, 4))
    
    results = []
    
    progress_bar = st.progress(0)
    total_steps = len(bit_range)
    
    for i, bit in enumerate(bit_range):
        st.write(f"Тестирование разрядности: **{bit} бит**...")
        
        # Генерируем тестовые числа
        test_numbers = [generate_semiprime(bit) for _ in range(runs_per_bit)]
        
        for alg_name, alg_key in [
            ("Поллард (Rho)", "pollard_rho"),
            ("Поллард (p-1)", "pollard_p1"),
            ("КР базовый", "qs_basic"),
            ("КР оптимизированный", "qs_optimized"),
            ("КР AUTO", "qs_auto"),
            ("КР LPV", "qs_lpv"),
            ("КР MPQS", "qs_mpqs"),
            ("КР MPQS Параллельный", "qs_mpqs_parallel"),
        ]:
            # КР базовый неэффективен на числах > 32 бит — пропускаем
            if alg_key == "qs_basic" and bit > 32:
                results.append({"Разрядность (бит)": bit, "Алгоритм": alg_name, "Среднее время (мс)": None})
                continue
            total_time = 0
            
            for num in test_numbers:
                try:
                    res = requests.post(
                        "http://127.0.0.1:8000/api/factorize",
                        json={"number": str(num), "algorithm": alg_key}
                    )
                    if res.status_code == 200:
                        total_time += res.json()["time_ms"]
                except:
                    pass
            
            avg_time = total_time / runs_per_bit
            results.append({"Разрядность (бит)": bit, "Алгоритм": alg_name, "Среднее время (мс)": avg_time})
            
        progress_bar.progress((i + 1) / total_steps)

    # Строим график
    df = pd.DataFrame(results)
    
    st.success("Эксперимент завершен!")
    
    st.subheader("Линейный график зависимости времени от битности")
    
    fig = go.Figure()
    for alg in df["Алгоритм"].unique():
        alg_data = df[df["Алгоритм"] == alg]
        fig.add_trace(go.Scatter(
            x=alg_data["Разрядность (бит)"], 
            y=alg_data["Среднее время (мс)"],
            mode='lines+markers',
            name=alg
        ))
        
    fig.update_layout(
        xaxis_title="Разрядность числа (бит)",
        yaxis_title="Среднее время выполнения (мс)",
        hovermode="x unified"
    )
    
    st.plotly_chart(fig, use_container_width=True)
    
    # Таблица результатов
    st.subheader("Сводная таблица")
    pivot_df = df.pivot(index="Разрядность (бит)", columns="Алгоритм", values="Среднее время (мс)")
    st.dataframe(pivot_df)