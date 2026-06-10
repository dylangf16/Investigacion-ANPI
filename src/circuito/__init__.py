"""Mini-SPICE: punto de operación DC de circuitos no lineales con diodos.

Resuelve ``F(x) = 0`` (análisis nodal + ecuación de Shockley) mediante
Newton-Raphson amortiguado, resolviendo en cada paso una Jacobiana dispersa
sin invertirla.  Ver ``docs/explicacion.md``.

Ejemplo
-------
>>> from circuito import Red, ModeloDiodo, Ensamblador, newton_amortiguado
>>> red = Red(tierra=0)
>>> red.fuente_dc("Vs", 10.0)
>>> red.resistencia("Vs", "N1", 1e3)
>>> red.diodo("N1", 0)
>>> ens = Ensamblador(red)
>>> r = newton_amortiguado(ens)
>>> r.exito
True
"""
from __future__ import annotations

from .constantes import IS_DEFECTO, N_DEFECTO, VT_300K, voltaje_overflow
from .dispositivos import Diodo, ModeloDiodo, Resistencia
from .ensamblador import Ensamblador
from .newton import (
    Resultado,
    newton_amortiguado,
    newton_puro,
)
from .red import Red

__all__ = [
    "Red",
    "Ensamblador",
    "Resistencia",
    "Diodo",
    "ModeloDiodo",
    "Resultado",
    "newton_puro",
    "newton_amortiguado",
    "voltaje_overflow",
    "VT_300K",
    "IS_DEFECTO",
    "N_DEFECTO",
]
