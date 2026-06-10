"""Pruebas del circuito exigente de 6 nodos y 5 diodos (ejemplos/exigente.py)."""
import numpy as np
import pytest

from circuito import Ensamblador, newton_amortiguado, newton_puro
from circuito.simbolico import verificar_jacobiana


def test_newton_puro_desborda_en_primera_iteracion():
    from exigente import construir
    ens = Ensamblador(construir())
    r = newton_puro(ens)
    assert not r.exito
    # Arranque frío con Vs=15 V: la 1.ª iteración ya produce Inf.
    assert r.iteraciones == 1
    assert not np.isfinite(r.residuo_final)


def test_newton_amortiguado_converge():
    from exigente import construir
    ens = Ensamblador(construir())
    r = newton_amortiguado(ens)
    assert r.exito
    assert np.linalg.norm(ens.F(r.x)) < 1e-8


def test_divisor_mantiene_d5_en_inversa():
    from exigente import construir
    ens = Ensamblador(construir())
    r = newton_amortiguado(ens)
    v = ens.voltajes_nodales(r.x)
    # N6 = N2 * R6/(R2+R6) = N2/3  =>  D5 (ánodo N6, cátodo N2) en inversa.
    assert v["N6"] / v["N2"] == pytest.approx(1.0 / 3.0, rel=1e-3)
    assert v["N6"] < v["N2"]


def test_jacobiana_casi_tridiagonal():
    from exigente import construir
    ens = Ensamblador(construir())
    r = newton_amortiguado(ens)
    J = ens.J(r.x).toarray()
    # Las únicas entradas fuera de la banda tridiagonal son (2,6) y (6,2)
    # (0-based: (1,5) y (5,1)), aportadas por R2 y D5 entre N2 y N6.
    n = J.shape[0]
    extra = set()
    for i in range(n):
        for j in range(n):
            if abs(i - j) > 1 and abs(J[i, j]) > 0:
                extra.add((i, j))
    assert extra == {(1, 5), (5, 1)}


def test_jacobiana_simbolica_vs_estampada():
    from exigente import construir
    assert verificar_jacobiana(construir(), puntos=4)
