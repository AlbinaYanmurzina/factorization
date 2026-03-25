factorization/
│
├── backend/                        # Бэкенд на FastAPI (Математика и API)
│   ├── algorithms/                 # Ядро: реализация алгоритмов
│   │   ├── __init__.py
│   │   ├── base.py                 # Базовый класс для алгоритмов (шаблоны)
│   │   ├── pollard.py              # Методы Полларда (Rho, p-1)
│   │   ├── quadratic_sieve.py      # Квадратичное решето (QS / MPQS)
│   │   └── math_utils.py           # Вспомогательная математика (НОД, простые числа, Лежандр)
│   │
│   ├── api/                        # Роутеры FastAPI
│   │   ├── __init__.py
│   │   ├── endpoints.py            # Эндпоинты для вызова факторизации
│   │   └── benchmark.py            # Эндпоинты для запуска массовых тестов
│   │
│   ├── schemas/                    # Pydantic модели (структуры данных API)
│   │   ├── __init__.py
│   │   └── models.py               # Запросы (Request) и ответы (Response)
│   │
│   └── main.py                     # Точка входа FastAPI
│
├── frontend/                       # Фронтенд на Streamlit (Веб-интерфейс)
│   ├── pages/                      # Многостраничность Streamlit
│   │   ├── 1_Step_by_Step.py       # Страница: Пошаговая визуализация
│   │   └── 2_Benchmark.py          # Страница: Графики и сравнение алгоритмов
│   │
│   ├── utils/                      # Вспомогательные функции для UI
│   │   └── api_client.py           # Клиент для запросов к FastAPI
│   │
│   └── app.py                      # Главная страница (Описание работы, вводные)
│
├── tests/                          # Unit-тесты (очень важны для математики!)
│   ├── test_pollard.py
│   └── test_qs.py
│
├── requirements.txt                # Зависимости проекта
├── .gitignore
└── README.md