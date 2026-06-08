"""Circuito pequeño de 4 nodos de ``docs/explicacion.md``.

    Vs=10V --[R1=1k]-- N1 --[D1]--> N3
                       |             |
                     [R2=1k]       [R4=2k]
                       |             |
                       N2          tierra
                       |
                     [R3=1k]
                       |
                       N3 --[D3]--> N4 --[R5=1k]-- tierra
                       N2 --[D2]--> N4

Construye exactamente el sistema F(x)=0 y la Jacobiana descritos en el documento
(diodos D1: N1->N3, D2: N2->N4, D3: N3->N4).  Útil para validar a mano.
"""
from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from circuito import Ensamblador, Red, newton_amortiguado, newton_puro  # noqa: E402


def construir(Vs: float = 10.0) -> Red:
    red = Red(tierra=0)
    red.fuente_dc("Vs", Vs)
    red.resistencia("Vs", "N1", 1e3)   # R1
    red.resistencia("N1", "N2", 1e3)   # R2
    red.resistencia("N2", "N3", 1e3)   # R3
    red.resistencia("N3", 0, 2e3)      # R4
    red.resistencia("N4", 0, 1e3)      # R5
    red.diodo("N1", "N3")              # D1
    red.diodo("N2", "N4")              # D2
    red.diodo("N3", "N4")              # D3
    return red


def main() -> None:
    red = construir()
    ens = Ensamblador(red)
    print(red)

    print("\n--- Newton puro (debería desbordar desde arranque frío) ---")
    print(newton_puro(ens))

    print("\n--- Newton amortiguado ---")
    r = newton_amortiguado(ens)
    print(r)
    if r.exito:
        print("\nVoltajes nodales (V):")
        for nodo, V in ens.voltajes_nodales(r.x).items():
            print(f"  {nodo:>4}: {V: .6f}")


if __name__ == "__main__":
    main()
