"""
backend/main.py

Точка входа FastAPI приложения для факторизации целых чисел.

Предоставляет REST API endpoint /api/factorize для выполнения факторизации
с использованием различных алгоритмов из учебного пособия (главы 3 и 6).

Архитектура:
- FastAPI для создания REST API
- Асинхронное выполнение с таймаутом через asyncio
- Маппинг строковых идентификаторов на классы алгоритмов
- Детальное логирование каждого шага для образовательных целей
"""

from fastapi import FastAPI, HTTPException
from schemas.models import FactorizeRequest, FactorizeResponse

# Импорт алгоритмов из главы 3 учебного пособия (простые алгоритмы факторизации)
from algorithms.pollard import PollardRho              # 3.4. ρ-метод Полларда
from algorithms.pollard_p1 import PollardP1            # 3.2. (p-1)-метод Полларда
from algorithms.fermat import FermatFactorization      # 3.1. Метод Ферма
from algorithms.williams_p1 import WilliamsP1          # 3.3. (p+1)-метод Вильямса
from algorithms.cfrac import CFRAC                     # 3.6. Факторизация непрерывными дробями
from algorithms.squfof import SQUFOF                   # 3.8. Факторизация квадратичными формами

# Импорт алгоритмов из главы 6 учебного пособия (метод квадратичного решета)
from algorithms.quadratic_sieve_basic import QuadraticSieveBasic          # 6.1. Алгоритм Диксона
from algorithms.quadratic_sieve_optimized import QuadraticSieveOptimized  # 6.2. Метод Померанца
from algorithms.quadratic_sieve_auto import QuadraticSieveAuto            # 6.5. С оценкой сложности
from algorithms.quadratic_sieve_lpv import QuadraticSieveLPV              # 6.8. Вариация множителя
from algorithms.quadratic_sieve_mpqs import QuadraticSieveMPQS            # 6.9. С множеством полиномов
from algorithms.quadratic_sieve_mpqs_parallel import QuadraticSieveMPQSParallel  # 6.9. Параллельная версия

import time
import asyncio

# Инициализация FastAPI приложения
app = FastAPI(title="ВКР: API Факторизации")

# Словарь маппинга строковых идентификаторов на классы алгоритмов
# Ключи используются frontend для выбора алгоритма
# Значения - классы, реализующие интерфейс FactorizationAlgorithm
ALGO_MAP = {
    "pollard_rho":       PollardRho,           # ρ-метод Полларда (раздел 3.4)
    "pollard_p1":        PollardP1,            # (p-1)-метод Полларда (раздел 3.2)
    "fermat":            FermatFactorization,  # Метод Ферма (раздел 3.1)
    "williams_p1":       WilliamsP1,           # (p+1)-метод Вильямса (раздел 3.3)
    "cfrac":             CFRAC,                # Факторизация непрерывными дробями (раздел 3.6)
    "squfof":            SQUFOF,               # Факторизация квадратичными формами (раздел 3.8)
    "qs_basic":          QuadraticSieveBasic,  # Алгоритм Диксона (раздел 6.1)
    "qs_optimized":      QuadraticSieveOptimized,  # Метод Померанца (раздел 6.2)
    "qs_auto":           QuadraticSieveAuto,   # С оценкой сложности (раздел 6.5)
    "qs_lpv":            QuadraticSieveLPV,    # Вариация множителя (раздел 6.8)
    "qs_mpqs":           QuadraticSieveMPQS,   # С множеством полиномов (раздел 6.9)
    "qs_mpqs_parallel":  QuadraticSieveMPQSParallel,  # Параллельная версия (раздел 6.9)
}

# Таймаут выполнения факторизации в секундах
# Предотвращает зависание на слишком больших числах
TIMEOUT_SECONDS = 30.0

@app.post("/api/factorize", response_model=FactorizeResponse)
async def factorize(request: FactorizeRequest):
    """
    API endpoint для факторизации целых чисел.
    
    Принимает:
        request (FactorizeRequest): Запрос с полями:
            - number (str): Число для факторизации (строка для поддержки больших чисел)
            - algorithm (str): Идентификатор алгоритма из ALGO_MAP
    
    Возвращает:
        FactorizeResponse: Ответ с полями:
            - factors (List[str]): Список простых множителей
            - time_ms (float): Время выполнения в миллисекундах
            - steps (List[Dict]): Детальный лог работы алгоритма
    
    Исключения:
        HTTPException(400): Некорректное число или неизвестный алгоритм
        HTTPException(500): Внутренняя ошибка при выполнении
    
    Особенности:
        - Выполнение в отдельном потоке через asyncio.to_thread()
        - Таймаут TIMEOUT_SECONDS для предотвращения зависания
        - Числа передаются как строки для поддержки произвольной длины
    """
    # Валидация входного числа
    try:
        n = int(request.number)
        if n < 2:
            raise ValueError("Число должно быть больше 1")
    except ValueError:
        raise HTTPException(status_code=400, detail="Некорректный ввод числа")

    # Получение класса алгоритма по идентификатору
    algo_cls = ALGO_MAP.get(request.algorithm)
    if algo_cls is None:
        raise HTTPException(status_code=400, detail="Неизвестный алгоритм")

    # Создание экземпляра алгоритма
    algo = algo_cls()
    start_time = time.perf_counter()

    try:
        # Асинхронное выполнение факторизации в отдельном потоке с таймаутом
        # asyncio.to_thread() позволяет не блокировать event loop FastAPI
        raw_factors = await asyncio.wait_for(
            asyncio.to_thread(algo.factorize, n),
            timeout=TIMEOUT_SECONDS
        )
        # Конвертация множителей в строки для поддержки больших чисел
        factors = [str(f) for f in raw_factors]
        # Получение лога шагов из алгоритма
        steps = algo.steps_log
        # Вычисление времени выполнения в миллисекундах
        execution_time = (time.perf_counter() - start_time) * 1000
    except asyncio.TimeoutError:
        # Обработка таймаута: возвращаем исходное число и сообщение об ошибке
        execution_time = TIMEOUT_SECONDS * 1000
        factors = [str(n)]
        steps = [{"step": "Тайм-аут", "details": {
            "message": f"Превышено время ожидания ({int(TIMEOUT_SECONDS)}с). Алгоритм не справился с числом данной разрядности."
        }}]
    except Exception as e:
        # Обработка других исключений
        raise HTTPException(status_code=500, detail=str(e))

    # Формирование и возврат ответа
    return FactorizeResponse(
        factors=factors,
        time_ms=execution_time,
        steps=steps
    )

if __name__ == "__main__":
    import uvicorn
    # Запуск сервера для локальной разработки
    # host="127.0.0.1" - доступ только с локальной машины
    # port=8000 - стандартный порт для FastAPI
    # reload=True - автоматическая перезагрузка при изменении кода
    uvicorn.run("main:app", host="127.0.0.1", port=8000, reload=True)