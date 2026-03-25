import pytest
from backend.algorithms.pollard import PollardRho

def test_pollard_rho():
    algo = PollardRho()
    
    # 1. Простое число (должно вернуть само себя)
    assert algo.factorize(97) == [97]
    
    # 2. Небольшое составное число (8051 = 83 * 97)
    assert algo.factorize(8051) == [83, 97]
    
    # 3. Число с несколькими делителями (2 * 3 * 3 * 5 = 90)
    assert algo.factorize(90) == [2, 3, 3, 5]
    
    # 4. Чуть большее число
    # 10403 = 101 * 103
    assert algo.factorize(10403) == [101, 103]

def test_pollard_logs():
    algo = PollardRho()
    algo.factorize(8051)
    
    # Проверяем, что логи вообще собрались
    assert len(algo.steps_log) > 0
    # Последний шаг должен содержать финальный ответ
    assert algo.steps_log[-1]["step"] == "Факторизация завершена"