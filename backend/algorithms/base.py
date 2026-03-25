# backend/algorithms/base.py
from abc import ABC, abstractmethod
from typing import List, Dict, Any

class FactorizationAlgorithm(ABC):
    def __init__(self):
        # Очищаем логи при каждом новом запуске
        self.steps_log: List[Dict[str, Any]] = []

    def log_step(self, step_name: str, details: Any):
        """
        Сохраняем информацию о шаге работы алгоритма.
        Эти данные FastAPI отдаст на фронтенд Streamlit.
        """
        self.steps_log.append({
            "step": step_name,
            "details": details
        })

    def clear_logs(self):
        self.steps_log = []

    @abstractmethod
    def factorize(self, n: int) -> List[int]:
        """
        Главный метод, который должен вернуть список простых множителей числа n.
        """
        pass