import pytest
from backend.algorithms.math_utils import is_prime, generate_primes, legendre_symbol

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