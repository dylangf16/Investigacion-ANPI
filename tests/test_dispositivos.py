"""Pruebas unitarias de los modelos de dispositivo y el ensamblaje."""
import numpy as np
import pytest

from circuito import Diodo, Ensamblador, ModeloDiodo, Red, Resistencia


def test_resistencia_rechaza_valor_no_positivo():
    with pytest.raises(ValueError):
        Resistencia("a", "b", 0.0)


def test_divisor_resistivo_lineal():
    # Vs=10 - R1=1k - N - R2=1k - tierra  =>  N = 5 V exacto.
    red = Red(0)
    red.fuente_dc("Vs", 10.0)
    red.resistencia("Vs", "N", 1e3)
    red.resistencia("N", 0, 1e3)
    ens = Ensamblador(red)
    F = ens.F(np.array([5.0]))
    assert abs(F[0]) < 1e-12


def test_diodo_corriente_y_conductancia():
    d = Diodo("a", "k", ModeloDiodo(Is=1e-14, n=1.0, Vt=0.025852))
    # En Vd=0: corriente nula, conductancia = Is/(n Vt).
    i0, g0 = d.corriente_y_conductancia(0.0)
    assert i0 == pytest.approx(0.0, abs=1e-20)
    assert g0 == pytest.approx(1e-14 / 0.025852, rel=1e-9)
    # La conductancia crece monótonamente con Vd.
    _, g_on = d.corriente_y_conductancia(0.7)
    assert g_on > g0 * 1e9


def test_exp_no_lanza_overflow_sino_inf():
    d = Diodo("a", "k")
    i, g = d.corriente_y_conductancia(50.0)  # exp(~1934) -> inf, sin excepción
    assert np.isinf(i) and np.isinf(g)


def test_estructura_dispersa_de_la_rejilla():
    # La Jacobiana de una rejilla debe ser dispersa (5 puntos), no densa.
    from rejilla import construir_rejilla
    red = construir_rejilla(10)
    ens = Ensamblador(red)
    J = ens.J(np.zeros(ens.n))
    densidad = J.nnz / (ens.n * ens.n)
    assert densidad < 0.05  # muy dispersa
