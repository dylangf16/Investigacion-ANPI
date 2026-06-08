r"""Equivalencia entre ``circuito.newton_puro`` y el método aprobado por el profe.

El archivo ``docs/p3_av2.m`` (Avance 2, Pregunta 3) implementa Newton-Raphson
así, en cada iteración:

    er = norm(F(x_k), 2)          % error = norma 2 del residuo
    if er < tol: break            % criterio de parada
    y  = J(x_k) \ F(x_k)          % sistema lineal SIN invertir (mldivide)
    x_k = x_k - y                 % actualización

Que es idéntico a ``dx = resolver(J, -F); x = x + dx`` con parada ``||F||2 < tol``.
Aquí se transcribe fielmente la referencia del profe en NumPy y se comprueba que
mi ``newton_puro`` produce la MISMA trayectoria (iterado a iterado) sobre los
tres sistemas de prueba del ``.m``.
"""
import numpy as np
import pytest
import scipy.sparse as sp

from circuito.newton import newton_puro


# --- Referencia fiel al p3_av2.m -------------------------------------------
def newton_profe(x0, F_func, J_func, tol, iter_max):
    """Transcripción 1:1 del ``newton_raphson`` de docs/p3_av2.m."""
    x = np.array(x0, dtype=float)
    er = []
    i = 0
    for i in range(1, iter_max + 1):          # for i = 1:iterMax
        f = F_func(x)
        er.append(np.linalg.norm(f, 2))       # er_k(end+1) = norm(f_eval, 2)
        if er[-1] < tol:                      # if er_k(end) < tol: break
            break
        y = np.linalg.solve(J_func(x), f)     # y = mldivide(J, f)
        x = x - y                             # x_k = x_k - y
    return x, i, er


# --- Adaptador para que newton_puro consuma una F/J genérica ---------------
class SistemaGenerico:
    """Expone la interfaz que espera ``newton_puro`` (``n``, ``evaluar``, ``F``)."""

    def __init__(self, F_func, J_func, n):
        self.n = n
        self._F = F_func
        self._J = J_func

    def evaluar(self, x):
        x = np.asarray(x, dtype=float)
        return self._F(x), sp.csc_matrix(self._J(x))

    def F(self, x):
        return self._F(np.asarray(x, dtype=float))


# --- Los tres sistemas de docs/p3_av2.m ------------------------------------
def _sistema_1():
    def F(x):
        x1, x2 = x
        return np.array([x1**2 - 2*x1 - x2 + 0.5, x1**2 + 4*x2**2 - 4])

    def J(x):
        x1, x2 = x
        return np.array([[2*x1 - 2, -1.0], [2*x1, 8*x2]])
    return [3.0, 2.0], F, J


def _sistema_2():
    def F(x):
        x1, x2 = x
        return np.array([np.sin(x1) + x2*np.cos(x1), x1 - x2])

    def J(x):
        x1, x2 = x
        return np.array([[np.cos(x1) - x2*np.sin(x1), np.cos(x1)], [1.0, -1.0]])
    return [1.2, -1.5], F, J


def _sistema_3():
    def F(x):
        x1, x2, x3, x4 = x
        return np.array([
            x2*x3 + x4*(x2 + x3),
            x1*x3 + x4*(x1 + x3),
            x1*x2 + x4*(x1 + x2),
            x1*x2 + x1*x3 + x2*x3 - 1,
        ])

    def J(x):
        x1, x2, x3, x4 = x
        return np.array([
            [0.0,      x3 + x4,  x2 + x4,  x2 + x3],
            [x3 + x4,  0.0,      x1 + x4,  x1 + x3],
            [x2 + x4,  x1 + x4,  0.0,      x1 + x2],
            [x2 + x3,  x1 + x3,  x1 + x2,  0.0],
        ])
    return [-1.0, -1.0, -1.0, -1.0], F, J


SISTEMAS = {"sistema_1": _sistema_1, "sistema_2": _sistema_2, "sistema_3": _sistema_3}


@pytest.mark.parametrize("nombre", list(SISTEMAS))
def test_newton_puro_igual_que_profe(nombre):
    tol, iter_max = 1e-5, 1000
    x0, F, J = SISTEMAS[nombre]()

    x_ref, i_ref, er_ref = newton_profe(x0, F, J, tol, iter_max)

    sistema = SistemaGenerico(F, J, n=len(x0))
    r = newton_puro(sistema, x0, tol=tol, max_iter=iter_max, metodo_lineal="lu")

    # Mismo punto final (ambos resuelven los mismos sistemas lineales exactos).
    assert np.allclose(r.x, x_ref, atol=1e-10), f"{r.x} != {x_ref}"
    # Misma trayectoria de error (norma 2 del residuo en cada iteración).
    assert len(r.historial_residuo) == len(er_ref)
    assert np.allclose(r.historial_residuo, er_ref, rtol=1e-10, atol=1e-12)
    # Convención de conteo: el profe cuenta desde 1, newton_puro desde 0.
    assert r.iteraciones == i_ref - 1
    # Ambos llegan por debajo de la tolerancia.
    assert r.exito and er_ref[-1] < tol
