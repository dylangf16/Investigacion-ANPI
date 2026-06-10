"""Circuito exigente de 6 nodos y 5 diodos (multiplicador con diodo de realimentación).

Topología (ver imagen del enunciado):

    Vs --[Rs]-- N1 --|D1>-- N2 --|D2>-- N3 --|D3>-- N4 --|D4>-- N5 (Vout)
                 |          |             |             |          |
               [R1]       [R2]          [R3]          [R4]       [R5]
                |          |             |             |          |
               gnd        N6 (D5 en inversa, de N6 a N2, paralelo a R2)
                          [R6]
                           |
                          gnd

* Cadena principal D1..D4: acopla nodos adyacentes -> Jacobiana tridiagonal.
* D5 (ánodo N6, cátodo N2) cierra una realimentación.  Como el divisor R2/R6
  mantiene ``N6 = N2 * R6/(R2+R6) = N2/3 < N2``, D5 queda en inversa y su
  conductancia dinámica es ~1e-52 S (prácticamente cero).  Sus aportes a la
  Jacobiana son las entradas (2,6) y (6,2): las dos únicas fuera de la banda
  tridiagonal.  Si una iteración sobreestima N6, D5 intenta conducir y esas
  entradas explotan: ahí se rompe el Newton puro.
* Arranque frío x=0 con Vs=15 V: el primer paso de Newton sin amortiguar cruza
  el umbral de overflow de la exponencial y produce Inf en la primera iteración.

Es el banco de pruebas ideal: si el Newton amortiguado converge aquí, el solver
está bien.
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

# Parámetros del enunciado.
VS = 15.0
RS = 50.0
R1, R2, R3, R4, R5, R6 = 10e3, 1e3, 5e3, 2e3, 1e3, 500.0
MODELO = ModeloDiodo(Is=1e-14, n=1.0, Vt=0.02585)


def construir() -> Red:
    red = Red(tierra=0)
    red.fuente_dc("Vs", VS)
    # Cadena principal primero, para fijar el orden de incógnitas N1..N5.
    red.resistencia("Vs", "N1", RS)        # Rs
    red.diodo("N1", "N2", MODELO)          # D1
    red.diodo("N2", "N3", MODELO)          # D2
    red.diodo("N3", "N4", MODELO)          # D3
    red.diodo("N4", "N5", MODELO)          # D4
    # Cargas a tierra.
    red.resistencia("N1", 0, R1)
    red.resistencia("N3", 0, R3)
    red.resistencia("N4", 0, R4)
    red.resistencia("N5", 0, R5)
    # Rama de realimentación (introduce N6 como 6.ª incógnita).
    red.resistencia("N2", "N6", R2)        # R2
    red.diodo("N6", "N2", MODELO)          # D5 (retro, en inversa)
    red.resistencia("N6", 0, R6)           # R6
    return red


def main() -> None:
    red = construir()
    ens = Ensamblador(red)
    print(red)
    print("Orden de incógnitas:", red.nodos_libres)

    print("\n--- Newton puro (debe desbordar en la 1.ª iteración) ---")
    print(newton_puro(ens))

    print("\n--- Newton amortiguado ---")
    r = newton_amortiguado(ens)
    print(r)

    if r.exito:
        print("\nVoltajes nodales (V):")
        v = ens.voltajes_nodales(r.x)
        for nodo in ("Vs", "N1", "N2", "N3", "N4", "N5", "N6"):
            print(f"  {nodo:>4}: {v[nodo]: .6f}")
        print(f"\n  Vout = V(N5) = {v['N5']:.6f} V")
        print(f"  Comprobación divisor: N6/N2 = {v['N6']/v['N2']:.4f}  (esperado ~0.3333)")


if __name__ == "__main__":
    main()
