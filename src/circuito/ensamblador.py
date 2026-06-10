"""Ensamblaje del residuo ``F(x)`` y de la Jacobiana dispersa ``J(x)``.

El ensamblador recorre los dispositivos de la :class:`~circuito.red.Red` y deja
que cada uno se "estampe".  La Jacobiana se construye en formato disperso COO y
se convierte a CSC para que el solucionador lineal la factorice sin invertirla.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from .red import Red


class Ensamblador:
    """Construye ``F(x)`` y ``J(x)`` por *stamping* sobre una :class:`Red`."""

    def __init__(self, red: Red) -> None:
        self.red = red
        red.finalizar()
        self.n = red.n_incognitas

        # Estado interno reutilizado durante un stamping (evita reasignar).
        self._x: np.ndarray | None = None
        self._F: np.ndarray | None = None
        self._filas: list[int] = []
        self._cols: list[int] = []
        self._datos: list[float] = []
        self._stamp_J: bool = True

    # -- API que ven los dispositivos (protocolo Contexto) ----------------
    def v(self, nodo: object) -> float:
        if nodo == self.red.tierra:
            return 0.0
        if nodo in self.red.fuentes:
            return self.red.fuentes[nodo]
        return float(self._x[self.red.indice[nodo]])

    def aporte_F(self, nodo: object, valor: float) -> None:
        idx = self.red.indice.get(nodo)
        if idx is not None:
            self._F[idx] += valor

    def aporte_J(self, fila: object, columna: object, valor: float) -> None:
        if not self._stamp_J:
            return
        i = self.red.indice.get(fila)
        if i is None:
            return
        j = self.red.indice.get(columna)
        if j is None:
            return
        self._filas.append(i)
        self._cols.append(j)
        self._datos.append(valor)

    # -- evaluación --------------------------------------------------------
    def F(self, x: np.ndarray) -> np.ndarray:
        """Solo el residuo (más barato; no arma la Jacobiana)."""
        self._x = np.asarray(x, dtype=float)
        self._F = np.zeros(self.n)
        self._stamp_J = False
        for d in self.red.dispositivos:
            d.evaluar(self)
        self._stamp_J = True
        return self._F

    def evaluar(self, x: np.ndarray) -> tuple[np.ndarray, sp.csc_matrix]:
        """Residuo y Jacobiana a la vez (comparte las evaluaciones de ``exp``)."""
        self._x = np.asarray(x, dtype=float)
        self._F = np.zeros(self.n)
        self._filas, self._cols, self._datos = [], [], []
        self._stamp_J = True

        for d in self.red.dispositivos:
            d.evaluar(self)

        J = sp.coo_matrix(
            (self._datos, (self._filas, self._cols)),
            shape=(self.n, self.n),
        ).tocsc()
        # COO suma duplicados automáticamente al convertir; perfecto para el
        # patrón de stamping donde varios dispositivos tocan la misma entrada.
        return self._F, J

    def J(self, x: np.ndarray) -> sp.csc_matrix:
        return self.evaluar(x)[1]

    # -- utilidades --------------------------------------------------------
    def voltajes_nodales(self, x: np.ndarray) -> dict[object, float]:
        """Diccionario nodo -> voltaje, incluyendo tierra y fuentes."""
        self._x = np.asarray(x, dtype=float)
        resultado: dict[object, float] = {self.red.tierra: 0.0}
        for nodo, V in self.red.fuentes.items():
            resultado[nodo] = V
        for nodo, idx in self.red.indice.items():
            resultado[nodo] = float(self._x[idx])
        return resultado
