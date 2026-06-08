"""Modelos de dispositivo y su *stamping* en el sistema no lineal.

Cada dispositivo sabe inyectar ("stamp", en jerga SPICE) su aporte tanto al
vector residuo ``F(x)`` como a la Jacobiana ``J(x)``.  El ensamblador
(``ensamblador.py``) recorre los dispositivos y delega en ``evaluar``.

Convención de signos
--------------------
La ecuación asociada a un nodo libre es la Ley de Corrientes de Kirchhoff
escrita como "suma de corrientes que salen del nodo = 0".  Para una corriente
``I`` que circula del nodo ``a`` al nodo ``b`` se aporta ``+I`` a la ecuación de
``a`` y ``-I`` a la de ``b``.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np

from .constantes import IS_DEFECTO, N_DEFECTO, VT_300K


class Contexto(Protocol):
    """Interfaz que el ensamblador expone a cada dispositivo."""

    def v(self, nodo) -> float:
        """Voltaje del nodo (0 si es tierra, fijo si es fuente, incógnita si libre)."""

    def aporte_F(self, nodo, valor: float) -> None:
        """Acumula ``valor`` en la ecuación KCL del nodo (ignorado si no es libre)."""

    def aporte_J(self, fila, columna, valor: float) -> None:
        """Acumula ``valor`` en ``J[fila, columna]`` (ignorado si alguno no es libre)."""


class Dispositivo(Protocol):
    """Todo dispositivo conoce sus nodos y sabe estamparse."""

    def nodos(self) -> tuple:
        ...

    def evaluar(self, ctx: Contexto) -> None:
        ...


@dataclass(frozen=True)
class Resistencia:
    """Resistencia lineal de ``R`` ohmios entre ``a`` y ``b``.

    Aporta el Laplaciano del grafo: es la parte lineal y dispersa de la
    Jacobiana, idéntica a la estructura de las matrices de diferencias finitas.
    """

    a: object
    b: object
    R: float

    def __post_init__(self) -> None:
        if self.R <= 0:
            raise ValueError(f"La resistencia debe ser positiva, recibí R={self.R}")

    def nodos(self) -> tuple:
        return (self.a, self.b)

    def evaluar(self, ctx: Contexto) -> None:
        g = 1.0 / self.R
        va, vb = ctx.v(self.a), ctx.v(self.b)
        corriente = g * (va - vb)  # de a hacia b

        ctx.aporte_F(self.a, corriente)
        ctx.aporte_F(self.b, -corriente)

        ctx.aporte_J(self.a, self.a, g)
        ctx.aporte_J(self.a, self.b, -g)
        ctx.aporte_J(self.b, self.a, -g)
        ctx.aporte_J(self.b, self.b, g)


@dataclass(frozen=True)
class ModeloDiodo:
    """Parámetros del modelo de Shockley de un diodo."""

    Is: float = IS_DEFECTO
    n: float = N_DEFECTO
    Vt: float = VT_300K

    @property
    def n_vt(self) -> float:
        return self.n * self.Vt

    @property
    def v_critico(self) -> float:
        """Voltaje de máxima curvatura, útil para *limiting* (Nagel, 1975)."""
        return self.n_vt * np.log(self.n_vt / (np.sqrt(2.0) * self.Is))


def _exp_seguro(z: np.ndarray | float) -> np.ndarray | float:
    """``exp`` que devuelve ``inf`` sin emitir ``RuntimeWarning``.

    No saturamos el resultado a propósito: el desbordamiento a ``inf`` es
    justamente el modo de falla del Newton puro que queremos exhibir.  El
    solucionador detecta el ``inf``/``nan`` y reacciona (ver ``newton.py``).
    """
    with np.errstate(over="ignore"):
        return np.exp(z)


@dataclass(frozen=True)
class Diodo:
    """Diodo de Shockley entre ánodo y cátodo.

    ``I = Is * (exp(Vd / (n*Vt)) - 1)`` con ``Vd = V(anodo) - V(catodo)``.
    Conductancia dinámica ``g = dI/dVd = (Is / (n*Vt)) * exp(Vd / (n*Vt))``.
    """

    anodo: object
    catodo: object
    modelo: ModeloDiodo = ModeloDiodo()

    def nodos(self) -> tuple:
        return (self.anodo, self.catodo)

    def corriente_y_conductancia(self, vd: float) -> tuple[float, float]:
        """Corriente y conductancia dinámica para una caída ``vd`` dada."""
        m = self.modelo
        e = _exp_seguro(vd / m.n_vt)
        corriente = m.Is * (e - 1.0)
        conductancia = (m.Is / m.n_vt) * e
        return corriente, conductancia

    def evaluar(self, ctx: Contexto) -> None:
        vd = ctx.v(self.anodo) - ctx.v(self.catodo)
        corriente, g = self.corriente_y_conductancia(vd)

        ctx.aporte_F(self.anodo, corriente)
        ctx.aporte_F(self.catodo, -corriente)

        ctx.aporte_J(self.anodo, self.anodo, g)
        ctx.aporte_J(self.anodo, self.catodo, -g)
        ctx.aporte_J(self.catodo, self.anodo, -g)
        ctx.aporte_J(self.catodo, self.catodo, g)
