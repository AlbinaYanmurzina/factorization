"""
backend/schemas/models.py

Pydantic модели для валидации запросов и ответов API.

Определяет структуру данных для:
- FactorizeRequest: входные данные для факторизации
- FactorizeResponse: результат факторизации с детальным логом

Pydantic обеспечивает:
- Автоматическую валидацию типов
- Сериализацию/десериализацию JSON
- Генерацию OpenAPI схемы для документации
"""

from pydantic import BaseModel
from typing import List, Dict, Any, Optional

class FactorizeRequest(BaseModel):
    """
    Модель запроса на факторизацию числа.
    
    Attributes:
        number (str): Число для факторизации в виде строки.
        algorithm (str): Идентификатор алгоритма из ALGO_MAP.
        b_override (Optional[int]): Ручная граница B для квадратичного решета.
                                    None → автоматический расчёт по L-нотации.
    """
    number: str
    algorithm: str
    b_override: Optional[int] = None

class FactorizeResponse(BaseModel):
    """
    Модель ответа с результатами факторизации.
    
    Attributes:
        factors (List[str]): Список простых множителей числа в виде строк.
                            Множители отсортированы по возрастанию.
                            Пример: ["2", "2", "3", "7"] для числа 84
                            Строковый формат для поддержки больших чисел.
        
        time_ms (float): Время выполнения факторизации в миллисекундах.
                        Измеряется через time.perf_counter() для точности.
                        Включает время работы алгоритма, но не сериализацию.
        
        steps (List[Dict[str, Any]]): Детальный лог работы алгоритма.
                                     Каждый шаг - словарь с полями:
                                     - "step" (str): Название шага
                                     - "details" (Dict): Детали шага, могут содержать:
                                       * "message" (str): Текстовое описание
                                       * "table" (List[Dict]): Табличные данные
                                       * "matrix_data" (List[List]): Матрицы
                                       * "FB" (List): Факторная база
                                       * Другие специфичные для алгоритма данные
                                     
                                     Используется frontend для пошаговой визуализации
                                     и образовательных целей.
    """
    factors: List[str]
    time_ms: float
    steps: List[Dict[str, Any]]