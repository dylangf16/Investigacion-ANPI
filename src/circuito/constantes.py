"""Constantes físicas y numéricas del modelo de circuitos con diodos.

Las referencias de los valores típicos provienen del modelo de Shockley y de
``docs/explicacion.md``.
"""
from __future__ import annotations

import math

#: Voltaje térmico a 300 K (V).  Vt = k*T/q.
VT_300K: float = 0.025852

#: Corriente de saturación inversa típica de un diodo de señal (A).
IS_DEFECTO: float = 1e-14

#: Factor de idealidad por defecto (adimensional, normalmente entre 1 y 2).
N_DEFECTO: float = 1.0

#: Argumento de ``exp`` a partir del cual se desborda la doble precisión.
#: ``exp(z)`` con ``z > 709.78`` produce ``inf`` en ``float64``.
ARG_OVERFLOW: float = 709.78


def voltaje_overflow(n: float = N_DEFECTO, vt: float = VT_300K) -> float:
    """Voltaje de unión que provoca overflow de ``exp`` en doble precisión.

    Para ``n=1`` da ~18.35 V; un paso de Newton sin amortiguar lo cruza con
    facilidad desde un arranque frío (ver ``docs/explicacion.md``, punto 2).
    """
    return ARG_OVERFLOW * n * vt
