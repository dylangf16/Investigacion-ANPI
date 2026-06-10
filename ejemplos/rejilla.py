"""Rejilla resistiva 2D con diodos: el circuito "enorme" del documento.

Malla de ``n x n`` nodos libres conectados por resistencias (Laplaciano de
5 puntos, idéntico a las matrices de diferencias finitas de los Avances 1/3/5).
La columna izquierda se alimenta de una fuente ``Vs`` y la derecha va a tierra,
ambas a través de resistencias.  En una fracción de los nodos se cuelga un
diodo a tierra, que introduce la no linealidad exponencial.

Para ``n=15`` da 225 incógnitas (tamaño exigente).  Subiendo ``n`` se
obtienen miles de incógnitas con Jacobiana dispersa.
"""
from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from circuito import (  # noqa: E402
    Ensamblador,
    ModeloDiodo,
    Red,
    newton_amortiguado,
)


def construir_rejilla(n: int = 15, *, Vs: float = 5.0, R: float = 1e3,
                      R_borde: float = 1e3, fraccion_diodos: float = 0.5,
                      modelo: ModeloDiodo | None = None) -> Red:
    """Construye la rejilla ``n x n`` con diodos de fuga a tierra.

    Parameters
    ----------
    n : lado de la malla (``n*n`` incógnitas).
    Vs : voltaje de la fuente que alimenta la columna izquierda.
    R : resistencia de cada rama interna de la malla.
    R_borde : resistencia de acople a la fuente y a tierra.
    fraccion_diodos : proporción de nodos con diodo de fuga (patrón regular).
    """
    modelo = modelo or ModeloDiodo()
    red = Red(tierra=0)
    red.fuente_dc("Vs", Vs)

    def nodo(i: int, j: int) -> tuple[int, int]:
        return (i, j)

    paso_diodo = max(1, round(1.0 / fraccion_diodos)) if fraccion_diodos > 0 else 0
    k = 0
    for i in range(n):
        for j in range(n):
            # Ramas hacia la derecha y hacia abajo (evita duplicar).
            if j + 1 < n:
                red.resistencia(nodo(i, j), nodo(i, j + 1), R)
            if i + 1 < n:
                red.resistencia(nodo(i, j), nodo(i + 1, j), R)
            # Acoples de borde: izquierda -> fuente, derecha -> tierra.
            if j == 0:
                red.resistencia("Vs", nodo(i, j), R_borde)
            if j == n - 1:
                red.resistencia(nodo(i, j), 0, R_borde)
            # Diodos de fuga en un patrón regular.
            if paso_diodo and (k % paso_diodo == 0):
                red.diodo(nodo(i, j), 0, modelo)
            k += 1
    return red


def main() -> None:
    for n in (4, 10, 30):
        red = construir_rejilla(n)
        ens = Ensamblador(red)
        r = newton_amortiguado(ens)
        print(f"n={n:>3}  incógnitas={red.n_incognitas:>4}  {r}")


if __name__ == "__main__":
    main()
