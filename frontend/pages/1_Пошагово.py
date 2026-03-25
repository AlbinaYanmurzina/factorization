import streamlit as st
import requests
import pandas as pd

st.set_page_config(page_title="Пошаговая визуализация", page_icon="🔍")

st.title("Пошаговая визуализация алгоритма")

# Блок ввода данных
with st.container():
    col1, col2 = st.columns([2, 1])
    with col1:
        # Для начала пробуем числа 8051 (маленькое) или 10403
        number_input = st.text_input("Введите составное число (n):", value="8051")
    with col2:
        algo_choice = st.selectbox(
            "Выберите алгоритм:",
            [
                "Поллард (Rho)",
                "Поллард (p-1)",
                "Квадратичное решето (базовый)",
                "Квадратичное решето (оптимизированный)",
                "Квадратичное решето (AUTO)",
                "Квадратичное решето (LPV)",
                "Квадратичное решето (MPQS)",
                "Квадратичное решето (MPQS Параллельный)",
            ]
        )

algo_map = {
    "Поллард (Rho)": "pollard_rho",
    "Поллард (p-1)": "pollard_p1",
    "Квадратичное решето (базовый)": "qs_basic",
    "Квадратичное решето (оптимизированный)": "qs_optimized",
    "Квадратичное решето (AUTO)": "qs_auto",
    "Квадратичное решето (LPV)": "qs_lpv",
    "Квадратичное решето (MPQS)": "qs_mpqs",
    "Квадратичное решето (MPQS Параллельный)": "qs_mpqs_parallel",
}

if st.button("Факторизовать", type="primary"):
    with st.spinner("Выполняются вычисления на сервере..."):
        try:
            # Отправляем запрос на наш FastAPI сервер
            response = requests.post(
                "http://127.0.0.1:8000/api/factorize",
                json={"number": number_input, "algorithm": algo_map[algo_choice]}
            )
            
            if response.status_code == 200:
                data = response.json()
                
                # Вывод результатов
                st.success(f"Успешно! Время работы: {data['time_ms']:.3f} мс")
                
                # Красиво выводим множители
                factors_str = " × ".join(data['factors'])
                st.markdown(f"### Результат: `{number_input} = {factors_str}`")
                
                st.divider()
                st.subheader("Шаги работы алгоритма")
                
                # Парсим логи, которые прислал бэкенд
                for idx, step_data in enumerate(data['steps']):
                    step_name = step_data.get('step', f'Шаг {idx+1}')
                    details = step_data.get('details', {})
                    
                    # Используем expander (раскрывающиеся блоки) для красоты
                    with st.expander(f"🔹 {step_name}", expanded=True if idx == 0 else False):
                        if 'message' in details:
                            st.write(details['message'])
                        
                        if 'table' in details and len(details['table']) > 0:
                            # Отрисовываем таблицу шагов (x, y, НОД) с помощью Pandas
                            df = pd.DataFrame(details['table'])
                            st.dataframe(df, use_container_width=True, hide_index=True)
                            
            else:
                st.error(f"Ошибка сервера: {response.json().get('detail', 'Неизвестная ошибка')}")
                
        except requests.exceptions.ConnectionError:
            st.error("Не удалось подключиться к серверу. Убедитесь, что FastAPI (backend) запущен.")