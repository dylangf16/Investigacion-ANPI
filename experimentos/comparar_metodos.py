"""Experimento numérico de la rúbrica: Newton puro vs amortiguado.

Genera la tabla que pide la sección de experimentos del paper: para cada método
reporta convergencia, número de iteraciones, residuo final y tiempo.  También
muestra cómo el número de condición de la Jacobiana se dispara cuando los
diodos encienden (el mal condicionamiento del documento).
"""
from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402

from circuito import (  # noqa: E402
    Ensamblador,
    newton_amortiguado,
    newton_puro,
)
from circuito import lineal  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "ejemplos"))
from exigente import construir as construir_exigente  # noqa: E402
from pequeno import construir as construir_pequeno  # noqa: E402
from rejilla import construir_rejilla  # noqa: E402


def _fila(r) -> str:
    estado = "converge" if r.exito else "FALLA"
    return (f"{r.metodo:<22} {estado:<9} "
            f"{r.iteraciones:>5} {r.residuo_final:>12.3e} {r.tiempo*1e3:>10.2f}")


def comparar(nombre: str, red) -> None:
    ens = Ensamblador(red)
    print(f"\n=== {nombre}  ({red.n_incognitas} incógnitas) ===")
    print(f"{'Método':<22} {'Estado':<9} {'iter':>5} {'||F||':>12} {'t (ms)':>10}")
    print("-" * 62)
    for solver in (newton_puro, newton_amortiguado):
        print(_fila(solver(ens)))


def condicionamiento(red) -> None:
    """Número de condición de la Jacobiana: frío vs en el punto de operación."""
    ens = Ensamblador(red)
    x_frio = np.zeros(ens.n)
    cond_frio = lineal.numero_condicion(ens.J(x_frio))
    r = newton_amortiguado(ens)
    cond_op = lineal.numero_condicion(ens.J(r.x))
    print(f"\nNúmero de condición de J  (red {red.n_incognitas}x{red.n_incognitas}):")
    print(f"  arranque frío (diodos apagados): {cond_frio:.3e}")
    print(f"  punto de operación (encendidos): {cond_op:.3e}")
    print(f"  factor de empeoramiento:         {cond_op / cond_frio:.1f}x")


def main() -> None:
    comparar("Circuito pequeño (4 nodos)", construir_pequeno())
    comparar("Circuito exigente (6 nodos, 5 diodos)", construir_exigente())
    comparar("Rejilla 4x4", construir_rejilla(4))
    comparar("Rejilla 20x20", construir_rejilla(20))
    condicionamiento(construir_exigente())


if __name__ == "__main__":
    main()
