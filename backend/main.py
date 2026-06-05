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


from algorithms.pollard import PollardRho #  ρ-метод Полларда
from algorithms.pollard_p1 import PollardP1 # (p-1)-метод Полларда

from algorithms.quadratic_sieve_basic import QuadraticSieveBasic # Алгоритм Диксона 

import time
import asyncio

# Инициализация FastAPI приложения
app = FastAPI(title="ВКР: API Факторизации")

# Словарь маппинга строковых идентификаторов на классы алгоритмов
# Ключи используются frontend для выбора алгоритма
# Значения - классы, реализующие интерфейс FactorizationAlgorithm
ALGO_MAP = {
    "pollard_rho":       PollardRho,           # ρ-метод Полларда 
    "pollard_p1":        PollardP1,            # (p-1)-метод Полларда 
    "qs_basic":          QuadraticSieveBasic,  # Алгоритм Диксона
}

# Таймаут выполнения факторизации в секундах
# Предотвращает зависание на слишком больших числах
# Значение по умолчанию для чисел без указания битности
DEFAULT_TIMEOUT_SECONDS = 240.0

def get_timeout_for_number(n: int) -> float:
    """
    Возвращает таймаут в секундах в зависимости от размера числа.
    
    Логика:
    - До 80 бит: 120 секунд (достаточно для всех алгоритмов)
    - 80-85 бит: 300 секунд (5 минут)
    - 85-92 бит: 600 секунд (10 минут, Диксон может занять ~10 минут)
    - 92-100 бит: 1200 секунд (20 минут, Диксон может занять ~20 минут)
    - 100-110 бит: 2400 секунд (40 минут, Диксон может занять ~40 минут)
    - 110+ бит: 3600 секунд (60 минут, Диксон может занять ~60 минут)
    """
    bits = n.bit_length()
    if bits < 80:
        return 120.0
    elif bits < 85:
        return 300.0
    elif bits < 92:
        return 600.0
    elif bits < 100:
        return 1200.0
    elif bits < 110:
        return 2400.0
    else:
        return 3600.0

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
    
    # Определяем таймаут в зависимости от размера числа
    timeout_seconds = get_timeout_for_number(n)
    
    start_time = time.perf_counter()

    try:
        # Передаём b_override только для квадратичного решета
        kwargs = {}
        if request.algorithm == "qs_basic" and request.b_override is not None:
            kwargs["b_override"] = request.b_override

        raw_factors = await asyncio.wait_for(
            asyncio.to_thread(algo.factorize, n, **kwargs),
            timeout=timeout_seconds
        )

        # Конвертация множителей в строки для поддержки больших чисел
        factors = [str(f) for f in raw_factors]
        # Получение лога шагов из алгоритма
        steps = algo.steps_log
        
        # Вычисление времени выполнения в миллисекундах
        execution_time = (time.perf_counter() - start_time) * 1000
    except asyncio.TimeoutError:
        # Обработка таймаута: возвращаем исходное число и сообщение об ошибке
        execution_time = timeout_seconds * 1000
        factors = [str(n)]
        steps = [{"step": "Тайм-аут", "details": {
            "message": f"Превышено время ожидания ({int(timeout_seconds)}с). Алгоритм не справился с числом данной разрядности."
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
    # reload=True - автоматическая перезагрузка при изменении кода
    uvicorn.run("main:app", host="127.0.0.1", port=8000, reload=True)