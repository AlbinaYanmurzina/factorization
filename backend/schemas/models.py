from pydantic import BaseModel
from typing import List, Dict, Any

class FactorizeRequest(BaseModel):
    # Передаем как строку, чтобы при росте чисел не упереться в лимиты JSON 
    number: str  
    algorithm: str # "pollard_rho", "pollard_p_1", "qs"

class FactorizeResponse(BaseModel):
    factors: List[str] # Делители тоже возвращаем как строки
    time_ms: float
    steps: List[Dict[str, Any]]