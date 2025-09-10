'''
улучшение алгоритма для sin(t)
'''
import math
def optimized_sin(t, delta=1e-3):
    two_pi = 2 * math.pi
    t_reduced = t % two_pi
    
    if t_reduced > math.pi:
        t_reduced = two_pi - t_reduced
        sign = -1
    else:
        sign = 1
    if t_reduced > math.pi/2:
        t_reduced = math.pi - t_reduced
    x = t_reduced
    term = x 
    result = term
    n = 1
    while abs(term) > delta:
        n += 2
        term = -term * x * x / ((n-1) * n)
        result += term 
    return sign * result
t_val = 10.5
print(f"math.sin({t_val}) = {math.sin(t_val)}")
print(f"optimized_sin({t_val}) = {optimized_sin(t_val)}")
print(f"Разница: {abs(math.sin(t_val) - optimized_sin(t_val))}")