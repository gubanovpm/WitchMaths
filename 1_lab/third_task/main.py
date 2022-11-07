import numpy as np
from math import exp

# Табличные Данные(в СИ)
z = np.array([0.0189, 0.1567, 0.8244])
M = np.array([0.028013, 0.04401, 0.03007])                      # кг/моль
Pc = np.multiply(np.array([33.944, 73.866, 48.839]), 10**(-5))  # Па
Tc = np.array([126.2, 304.7, 305.43])                           # К
OmegaA0 = np.array([0.45724, 0.45724, 0.45724])
OmegaB  = np.array([0.077796, 0.077796, 0.077796])
omega = np.array([0.04, 0.225, 0.0986])

T = np.array([373, 373, 373])
P = np.array([1, 1, 1])

Tr = np.array(list(map(lambda Ti, Tci: Ti/ Tci, T, Tc)))
Pr = np.array(list(map(lambda Pi, Pci: Pi/ Pci, P, Pc)))

# 1. Находим Ki для каждого элемента
K  = np.array(list(map(lambda omegai, tri, pci: exp(5.373 * (1 - omegai))*(1-1/tri)*pci, omega, Tr, Pc)))

# 2. Решимм уравнение -//- находим a
# Для решения уравнения воспользуемся методом Ньютона
zip(z, K)
G = lambda alpha: sum(list(map(lambda zi, Ki: (zi * (Ki - 1)) / (alpha*(Ki-1) + 1, z, K))))
derG = lambda alpha: sum(list(map(lambda zi, Ki: (zi * (Ki - 1)) / (alpha*(Ki-1) + 1, z, K))))
