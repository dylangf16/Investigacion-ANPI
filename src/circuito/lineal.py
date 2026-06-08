"""Solución del sistema lineal disperso ``J dx = b`` **sin invertir** ``J``.

El enunciado del Avance 2 P3 prohíbe calcular la inversa: en cada paso de
Newton hay que resolver el sistema lineal con un solver eficiente.  Aquí se
ofrecen dos estrategias, ambas sobre la matriz dispersa:

* ``"lu"``  — factorización LU dispersa (``scipy.sparse.linalg.splu``).  Es la
  vía directa: se refactoriza la Jacobiana en cada iteración, tal como hace
  SPICE y como narra ``docs/explicacion.md`` ("rearmar y refactorizar").
* ``"gmres"`` — método iterativo de Krylov, sin factorización, útil cuando la
  rejilla es enorme y la LU se vuelve costosa.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla


class ErrorLineal(RuntimeError):
    """La Jacobiana es singular o el solver iterativo no convergió."""


def resolver(J: sp.spmatrix, b: np.ndarray, metodo: str = "lu") -> np.ndarray:
    """Resuelve ``J x = b`` y devuelve ``x``.

    Parameters
    ----------
    J : matriz dispersa cuadrada.
    b : término independiente.
    metodo : ``"lu"`` (directo) o ``"gmres"`` (iterativo).
    """
    if metodo == "lu":
        return _resolver_lu(J, b)
    if metodo == "gmres":
        return _resolver_gmres(J, b)
    raise ValueError(f"Método lineal desconocido: {metodo!r}")


def _resolver_lu(J: sp.spmatrix, b: np.ndarray) -> np.ndarray:
    try:
        factor = spla.splu(J.tocsc())
    except RuntimeError as exc:  # matriz singular
        raise ErrorLineal(f"Jacobiana singular: {exc}") from exc
    return factor.solve(np.asarray(b, dtype=float))


def _resolver_gmres(J: sp.spmatrix, b: np.ndarray) -> np.ndarray:
    Jc = J.tocsc()
    # Precondicionador ILU: clave para que GMRES converja en la Jacobiana
    # mal condicionada (conductancias que varían ~11 órdenes de magnitud).
    try:
        ilu = spla.spilu(Jc)
        M = spla.LinearOperator(Jc.shape, ilu.solve)
    except RuntimeError:
        M = None
    x, info = spla.gmres(Jc, np.asarray(b, dtype=float), M=M, rtol=1e-10, atol=0.0)
    if info != 0:
        raise ErrorLineal(f"GMRES no convergió (info={info})")
    return x


def numero_condicion(J: sp.spmatrix) -> float:
    """Estimación del número de condición-2 (denso: solo para circuitos chicos).

    Sirve para los experimentos numéricos: muestra cómo se dispara el
    condicionamiento cuando los diodos encienden.
    """
    A = J.toarray()
    return float(np.linalg.cond(A))
