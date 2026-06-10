"""Escalera R-diodo de N nodos: Jacobiana NxN tridiagonal (por defecto N=50).

Netlist SPICE equivalente (mostrado para N=4; aquí se generaliza a N nodos)::

    Vs in 0 5
    R0 in n1 1k
    Rk n(k) n(k+1) 1k     (cadena de resistencias)
    Dk n(k) 0 DMOD        (diodo de fuga en cada nodo)
    .model DMOD D(IS=1e-14 N=1 RS=0)
    .op

La cadena de resistencias ``in -> n1 -> n2 -> ... -> nN`` acopla únicamente nodos
adyacentes, y cada diodo de fuga aporta solo a la diagonal.  Por eso la Jacobiana
es exactamente **tridiagonal** (la misma estructura de banda de las matrices de
diferencias finitas 1D del curso), pero el sistema es no lineal por los diodos.
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
    newton_puro,
)

# Parámetros del netlist.
VS = 5.0
R = 1e3
MODELO = ModeloDiodo(Is=1e-14, n=1.0)  # DMOD: IS=1e-14, N=1


def construir(N: int = 50) -> Red:
    """Escalera de ``N`` nodos: fuente -> R -> n1 -> R -> ... -> nN, con diodo
    de fuga a tierra en cada nodo."""
    red = Red(tierra=0)
    red.fuente_dc("in", VS)
    anterior = "in"
    for k in range(1, N + 1):                  # cadena resistiva R0..R(N-1)
        red.resistencia(anterior, f"n{k}", R)
        anterior = f"n{k}"
    for k in range(1, N + 1):                  # diodos de fuga D1..DN
        red.diodo(f"n{k}", 0, MODELO)
    return red


def main() -> None:
    red = construir(50)
    ens = Ensamblador(red)
    print(red)
    print("Orden de incógnitas:", red.nodos_libres)

    # La Jacobiana es tridiagonal: sin entradas fuera de la banda |i-j|<=1.
    import numpy as np
    J = ens.J(np.zeros(ens.n)).toarray()
    fuera_banda = np.abs(J[np.abs(np.subtract.outer(range(ens.n), range(ens.n))) > 1])
    print(f"¿Jacobiana tridiagonal? {np.all(fuera_banda == 0)}")

    print("\n--- Newton puro ---")
    print(newton_puro(ens))

    print("\n--- Newton amortiguado ---")
    r = newton_amortiguado(ens)
    print(r)
    if r.exito:
        print("\nVoltajes nodales (V):")
        for nodo, V in ens.voltajes_nodales(r.x).items():
            print(f"  {str(nodo):>4}: {V: .6f}")


if __name__ == "__main__":
    main()
