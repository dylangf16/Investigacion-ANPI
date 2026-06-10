"""Pruebas de la escalera R-diodo (ejemplos/escalera.py)."""
import numpy as np
import pytest

from circuito import Ensamblador, newton_amortiguado, newton_puro


def test_jacobiana_tridiagonal():
    from escalera import construir
    ens = Ensamblador(construir(4))
    J = ens.J(np.zeros(ens.n)).toarray()
    idx = np.abs(np.subtract.outer(range(ens.n), range(ens.n)))
    assert np.all(J[idx > 1] == 0)        # nada fuera de la banda tridiagonal
    assert np.all(np.diag(J) != 0)        # diagonal no nula


def test_amortiguado_converge_y_es_fisico():
    from escalera import construir
    ens = Ensamblador(construir(4))
    r = newton_amortiguado(ens)
    assert r.exito
    assert np.linalg.norm(ens.F(r.x)) < 1e-9
    v = ens.voltajes_nodales(r.x)
    # Diodos de fuga en conducción (~0.5 a 0.75 V) y voltaje decreciente.
    assert 0.5 < v["n1"] < 0.75
    assert v["n1"] > v["n2"] > v["n3"] > v["n4"] > 0


@pytest.mark.parametrize("N", [4, 8, 16])
def test_escala_a_N_nodos(N):
    from escalera import construir
    red = construir(N)
    assert len(red.nodos_libres) == N
    r = newton_amortiguado(Ensamblador(red))
    assert r.exito


def test_newton_puro_falla():
    from escalera import construir
    r = newton_puro(Ensamblador(construir(4)))
    assert not r.exito
