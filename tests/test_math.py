import pytest
from backend.algorithms.math_utils import is_prime, generate_primes, legendre_symbol, tonelli_shanks

def test_is_prime():
    assert is_prime(2) == True
    assert is_prime(97) == True
    assert is_prime(100) == False
    assert is_prime(104729) == True  # 10000-е простое число
    assert is_prime(104727) == False

def test_generate_primes():
    primes = generate_primes(30)
    assert primes == [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
    assert len(generate_primes(100)) == 25

def test_legendre_symbol():
    # 2^2 = 4 = 4 mod 7 -> 4 это вычет
    assert legendre_symbol(4, 7) == 1
    # Нет такого x, что x^2 = 5 mod 7 -> 5 это невычет
    assert legendre_symbol(5, 7) == -1
    assert legendre_symbol(14, 7) == 0

def test_tonelli_shanks():
    # Решаем x^2 ≡ 10 (mod 13). 
    # Ответы: 6 (т.к. 36 = 13*2 + 10) и 13 - 6 = 7 (т.к. 49 = 13*3 + 10)
    root = tonelli_shanks(10, 13)
    assert root in (6, 7)
    
    # Решаем x^2 ≡ 56 (mod 101)
    root2 = tonelli_shanks(56, 101)
    assert (root2 * root2) % 101 == 56

    # Проверка на отсутствие решения
    assert tonelli_shanks(5, 7) is None