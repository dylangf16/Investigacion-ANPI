"""Métodos de Newton-Raphson para ``F(x) = 0``.

Dos variantes, en orden creciente de robustez:

#. :func:`newton_puro` — el paso completo del Avance 2.  Se desborda
   (``inf``/``nan``) en cuanto un paso empuja una unión más allá de ~18 V.
#. :func:`newton_amortiguado` — Newton con búsqueda de línea (backtracking).
   Frena el paso hasta que ``||F||`` decrece; evita el overflow.

Ambas resuelven el sistema lineal con :mod:`circuito.lineal` (sin invertir).
"""
from __future__ import annotations

import time
from dataclasses import dataclass, field

import numpy as np

from . import lineal
from .ensamblador import Ensamblador


@dataclass
class Resultado:
    """Salida estándar de cualquier solucionador."""

    metodo: str
    exito: bool
    x: np.ndarray
    iteraciones: int
    historial_residuo: list[float] = field(default_factory=list)
    tiempo: float = 0.0
    mensaje: str = ""

    @property
    def residuo_final(self) -> float:
        return self.historial_residuo[-1] if self.historial_residuo else float("nan")

    def __repr__(self) -> str:
        estado = "OK" if self.exito else "FALLO"
        return (f"Resultado({self.metodo}: {estado}, "
                f"iter={self.iteraciones}, "
                f"||F||={self.residuo_final:.3e}, "
                f"t={self.tiempo*1e3:.2f} ms)")


def _norma(v: np.ndarray) -> float:
    # Con valores enormes el producto interno puede desbordar a inf; lo
    # dejamos pasar (es el síntoma del overflow) sin emitir RuntimeWarning.
    with np.errstate(over="ignore", invalid="ignore"):
        return float(np.linalg.norm(v, ord=2))


def _finito(v: np.ndarray) -> bool:
    return bool(np.all(np.isfinite(v)))


# ---------------------------------------------------------------------------
# 1. Newton puro
# ---------------------------------------------------------------------------
def newton_puro(ens: Ensamblador, x0: np.ndarray | None = None, *,
                tol: float = 1e-9, max_iter: int = 100,
                metodo_lineal: str = "lu") -> Resultado:
    """Newton sin globalizar: ``x <- x + dx``.  Diverge con diodos fríos."""
    x = np.zeros(ens.n) if x0 is None else np.array(x0, dtype=float)
    hist: list[float] = []
    t0 = time.perf_counter()

    for k in range(max_iter):
        F, J = ens.evaluar(x)
        nrm = _norma(F)
        hist.append(nrm)
        if not (_finito(F) and np.isfinite(nrm)):
            return Resultado("Newton puro", False, x, k, hist,
                             time.perf_counter() - t0,
                             "overflow: exp desbordó (||F|| no finito)")
        if nrm < tol:
            return Resultado("Newton puro", True, x, k, hist,
                             time.perf_counter() - t0, "convergió")
        try:
            dx = lineal.resolver(J, -F, metodo=metodo_lineal)
        except lineal.ErrorLineal as exc:
            return Resultado("Newton puro", False, x, k, hist,
                             time.perf_counter() - t0, str(exc))
        x = x + dx

    return Resultado("Newton puro", False, x, max_iter, hist,
                     time.perf_counter() - t0, "no convergió en max_iter")


# ---------------------------------------------------------------------------
# 2. Newton amortiguado (backtracking line search)
# ---------------------------------------------------------------------------
def newton_amortiguado(ens: Ensamblador, x0: np.ndarray | None = None, *,
                       tol: float = 1e-9, max_iter: int = 100,
                       metodo_lineal: str = "lu",
                       alpha_min: float = 1e-4,
                       c_armijo: float = 1e-4) -> Resultado:
    """Newton con búsqueda de línea por bisección sobre ``alpha``.

    En cada iteración se acepta el mayor ``alpha in (0, 1]`` tal que
    ``||F(x + alpha*dx)||`` decrezca lo suficiente (condición de Armijo sobre
    la mitad del cuadrado de la norma).  Esto es lo que evita el overflow del
    Newton puro y es la estrategia de damping de SPICE (Nagel, 1975).
    """
    x = np.zeros(ens.n) if x0 is None else np.array(x0, dtype=float)
    hist: list[float] = []
    t0 = time.perf_counter()

    F, J = ens.evaluar(x)
    nrm = _norma(F)
    hist.append(nrm)

    for k in range(max_iter):
        if nrm < tol:
            return Resultado("Newton amortiguado", True, x, k, hist,
                             time.perf_counter() - t0, "convergió")
        try:
            dx = lineal.resolver(J, -F, metodo=metodo_lineal)
        except lineal.ErrorLineal as exc:
            return Resultado("Newton amortiguado", False, x, k, hist,
                             time.perf_counter() - t0, str(exc))

        alpha = 1.0
        objetivo = 0.5 * nrm * nrm
        aceptado = False
        while alpha >= alpha_min:
            x_nuevo = x + alpha * dx
            F_nuevo = ens.F(x_nuevo)
            if _finito(F_nuevo):
                nrm_nuevo = _norma(F_nuevo)
                # Descenso suficiente (Armijo simplificado).
                if 0.5 * nrm_nuevo * nrm_nuevo <= (1.0 - c_armijo * alpha) * objetivo:
                    aceptado = True
                    break
            alpha *= 0.5

        if not aceptado:
            return Resultado("Newton amortiguado", False, x, k, hist,
                             time.perf_counter() - t0,
                             f"line search falló (alpha<{alpha_min:g})")

        x = x_nuevo
        F, J = ens.evaluar(x)
        nrm = _norma(F)
        hist.append(nrm)

    return Resultado("Newton amortiguado", False, x, max_iter, hist,
                     time.perf_counter() - t0, "no convergió en max_iter")