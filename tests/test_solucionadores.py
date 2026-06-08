"""Pruebas de los solucionadores de Newton y la verificación simbólica."""
import numpy as np
import pytest

from circuito import (
    Ensamblador,
    Red,
    newton_amortiguado,
    newton_continuacion,
    newton_puro,
)
from circuito.simbolico import punto_operacion_un_diodo, verificar_jacobiana


def _un_diodo(Vs=5.0, R=1e3):
    red = Red(0)
    red.fuente_dc("Vs", Vs)
    red.resistencia("Vs", "N", R)
    red.diodo("N", 0)
    return red


def test_un_diodo_coincide_con_lambert_w():
    V_exacto, _ = punto_operacion_un_diodo(5.0, 1e3)
    ens = Ensamblador(_un_diodo())
    r = newton_amortiguado(ens)
    assert r.exito
    assert ens.voltajes_nodales(r.x)["N"] == pytest.approx(V_exacto, rel=1e-7)


def test_newton_puro_falla_circuito_pequeno():
    from pequeno import construir
    ens = Ensamblador(construir())
    r = newton_puro(ens)
    assert not r.exito  # se desborda / Jacobiana singular


def test_newton_amortiguado_resuelve_circuito_pequeno():
    from pequeno import construir
    ens = Ensamblador(construir())
    r = newton_amortiguado(ens)
    assert r.exito
    # Residuo de KCL realmente cerca de cero.
    assert np.linalg.norm(ens.F(r.x)) < 1e-8
    # Caídas de diodo en rango físico (~0.6 a 0.75 V).
    v = ens.voltajes_nodales(r.x)
    assert 0.5 < v["N1"] - v["N3"] < 0.8


def test_continuacion_resuelve_donde_puro_falla():
    from pequeno import construir
    ens = Ensamblador(construir())
    r = newton_continuacion(ens)
    assert r.exito
    assert np.linalg.norm(ens.F(r.x)) < 1e-7


def test_amortiguado_y_continuacion_dan_misma_solucion():
    from pequeno import construir
    ens = Ensamblador(construir())
    ra = newton_amortiguado(ens)
    rc = newton_continuacion(ens)
    assert np.allclose(ra.x, rc.x, atol=1e-6)


def test_jacobiana_simbolica_vs_estampada():
    from pequeno import construir
    assert verificar_jacobiana(construir(), puntos=5)


@pytest.mark.parametrize("n", [4, 10, 20])
def test_rejilla_converge(n):
    from rejilla import construir_rejilla
    ens = Ensamblador(construir_rejilla(n))
    r = newton_amortiguado(ens)
    assert r.exito
    assert np.linalg.norm(ens.F(r.x)) < 1e-7


def test_solver_gmres_coincide_con_lu():
    from pequeno import construir
    ens = Ensamblador(construir())
    r_lu = newton_amortiguado(ens, metodo_lineal="lu")
    r_gmres = newton_amortiguado(ens, metodo_lineal="gmres")
    assert r_gmres.exito
    assert np.allclose(r_lu.x, r_gmres.x, atol=1e-6)
