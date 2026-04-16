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

# Импорт алгоритмов из главы 6 учебного пособия (метод квадратичного решета)
from algorithms.quadratic_sieve_basic import QuadraticSieveBasic          # 6.1. Алгоритм Диксона (QS Basic)

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
    "qs_basic":          QuadraticSieveBasic,  # Алгоритм Диксона / QS Basic (раздел 6.1)
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