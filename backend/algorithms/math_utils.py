import random
import math
from typing import List, Optional

def is_prime(n: int, k: int = 5) -> bool:
    """
    Вероятностный тест простоты Миллера-Рабина.
    k - количество раундов тестирования (точность).
    """
    if n < 2: return False
    if n in (2, 3): return True
    if n % 2 == 0: return False

    # Представляем n - 1 как 2^s * d, где d - нечетное
    s, d = 0, n - 1
    while d % 2 == 0:
        s += 1
        d //= 2

    # Проводим k раундов тестов
    for _ in range(k):
        a = random.randint(2, n - 2)
        x = pow(a, d, n) # a^d mod n
        
        if x == 1 or x == n - 1:
            continue
            
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            # Если цикл не прервался через break, число составное
            return False
            
    return True

def generate_primes(limit: int) -> List[int]:
    """
    Классическое Решето Эратосфена.
    Генерирует список простых чисел до заданного предела (limit).
    Используется для создания Факторной Базы.
    """
    if limit < 2:
        return []
        
    sieve = [True] * (limit + 1)
    sieve[0] = sieve[1] = False
    
    for i in range(2, int(math.isqrt(limit)) + 1):
        if sieve[i]:
            for j in range(i * i, limit + 1, i):
                sieve[j] = False
                
    return [i for i, is_p in enumerate(sieve) if is_p]

def legendre_symbol(a: int, p: int) -> int:
    """
    Вычисляет символ Лежандра (a/p) с помощью критерия Эйлера.
    Возвращает:
     1, если 'a' является квадратичным вычетом по модулю 'p'
    -1, если 'a' является квадратичным невычетом
     0, если a делится на p
    """
    if a % p == 0:
        return 0
        
    # Критерий Эйлера: a^((p-1)/2) mod p
    ls = pow(a, (p - 1) // 2, p)
    
    # В Python pow() возвращает положительные числа. 
    # Если результат p - 1, это эквивалентно -1 mod p
    return -1 if ls == p - 1 else ls

def tonelli_shanks(n: int, p: int) -> Optional[int]:
    """
    Алгоритм Тонелли-Шенкса.
    Находит x, такой что x^2 ≡ n (mod p).
    Возвращает одно из решений (второе будет p - x), либо None.
    Критически важен для инициализации решета в методе Квадратичного решета.
    """
    n = n % p
    
    # Простые случаи
    if n == 0: return 0
    if p == 2: return n
    if legendre_symbol(n, p) != 1:
        return None # Решений нет
        
    # Если p ≡ 3 (mod 4), корень вычисляется по простой формуле
    if p % 4 == 3:
        return pow(n, (p + 1) // 4, p)
        
    # Шаг 1: Находим Q и S, такие что p - 1 = Q * 2^S (где Q - нечетное)
    q = p - 1
    s = 0
    while q % 2 == 0:
        s += 1
        q //= 2
        
    # Шаг 2: Находим квадратичный невычет z по модулю p
    z = 2
    while legendre_symbol(z, p) != -1:
        z += 1
        
    # Шаг 3: Инициализация переменных
    m = s
    c = pow(z, q, p)
    t = pow(n, q, p)
    r = pow(n, (q + 1) // 2, p)
    
    # Шаг 4: Поиск цикла
    while t != 0 and t != 1:
        t2i = t
        i = 0
        for i in range(1, m):
            t2i = pow(t2i, 2, p)
            if t2i == 1:
                break
                
        # Вычисляем b = c^(2^(m-i-1)) mod p
        b = pow(c, pow(2, m - i - 1), p)
        
        m = i
        c = pow(b, 2, p)
        t = (t * c) % p
        r = (r * b) % p
        
    return r