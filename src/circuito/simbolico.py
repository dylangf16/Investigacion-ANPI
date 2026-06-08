"""Apoyo simbólico con SymPy.

Dos usos, ambos para *validar* la maquinaria numérica:

#. :func:`construir_simbolico` / :func:`verificar_jacobiana` — arma ``F(x)`` y
   su Jacobiana de forma exacta con SymPy y la compara contra la Jacobiana
   estampada a mano por el :class:`~circuito.ensamblador.Ensamblador`.  Así se
   garantiza que el stamping (propenso a errores de signo) es correcto.
#. :func:`punto_operacion_un_diodo` — solución en **forma cerrada** del caso de
   un único diodo con resistencia serie vía la función ``W`` de Lambert
   (Banwell & Jaski, 2000).  Es el "caso trivial" que sí tiene fórmula y que
   sirve de patrón de oro para el solver de Newton.
"""
from __future__ import annotations

import numpy as np
import sympy as sp

from .dispositivos import Diodo, ModeloDiodo, Resistencia
from .ensamblador import Ensamblador
from .red import Red


def construir_simbolico(red: Red):
    """Devuelve ``(variables, F, J)`` simbólicos de la red.

    ``variables`` es la lista de símbolos de los nodos libres en orden de
    incógnita; ``F`` es un ``sympy.Matrix`` columna; ``J`` su Jacobiana exacta.
    """
    red.finalizar()
    simbolos = {nodo: sp.Symbol(f"V_{i}", real=True)
                for nodo, i in red.indice.items()}

    def v(nodo):
        if nodo == red.tierra:
            return sp.Integer(0)
        if nodo in red.fuentes:
            return sp.Float(red.fuentes[nodo])
        return simbolos[nodo]

    F = [sp.Integer(0)] * red.n_incognitas

    def aporte(nodo, expr):
        idx = red.indice.get(nodo)
        if idx is not None:
            F[idx] += expr

    for d in red.dispositivos:
        if isinstance(d, Resistencia):
            g = sp.Rational(1) / sp.Float(d.R)
            corriente = g * (v(d.a) - v(d.b))
            aporte(d.a, corriente)
            aporte(d.b, -corriente)
        elif isinstance(d, Diodo):
            m: ModeloDiodo = d.modelo
            vd = v(d.anodo) - v(d.catodo)
            corriente = sp.Float(m.Is) * (sp.exp(vd / sp.Float(m.n_vt)) - 1)
            aporte(d.anodo, corriente)
            aporte(d.catodo, -corriente)
        else:
            raise TypeError(f"Dispositivo sin modelo simbólico: {type(d).__name__}")

    variables = [simbolos[nodo] for nodo in red.indice]
    Fvec = sp.Matrix(F)
    J = Fvec.jacobian(variables)
    return variables, Fvec, J


def verificar_jacobiana(red: Red, puntos: int = 5, tol: float = 1e-6,
                        semilla: int = 0) -> bool:
    """Compara la Jacobiana estampada con la simbólica en puntos aleatorios.

    Devuelve ``True`` si coinciden dentro de ``tol`` (relativo) en todos los
    puntos; de lo contrario lanza ``AssertionError`` con el detalle.
    """
    variables, _, J = construir_simbolico(red)
    J_func = sp.lambdify(variables, J, modules="numpy")

    ens = Ensamblador(red)
    rng = np.random.default_rng(semilla)

    for _ in range(puntos):
        # Voltajes modestos: evita el overflow simbólico-numérico de exp.
        x = rng.uniform(-0.5, 0.7, size=red.n_incognitas)
        J_num = np.atleast_2d(np.asarray(J_func(*x), dtype=float))
        J_stamp = ens.J(x).toarray()
        if not np.allclose(J_num, J_stamp, rtol=tol, atol=tol):
            diff = np.max(np.abs(J_num - J_stamp))
            raise AssertionError(
                f"Jacobiana estampada difiere de la simbólica (máx |Δ|={diff:.3e})\n"
                f"x={x}\nsimbólica=\n{J_num}\nstamp=\n{J_stamp}")
    return True


def punto_operacion_un_diodo(Vs: float, R: float,
                             modelo: ModeloDiodo | None = None) -> tuple[float, float]:
    """Forma cerrada (Lambert W) de ``Vs -[R]- nodo -[diodo]- tierra``.

    Resuelve ``(V - Vs)/R + Is*(exp(V/(n*Vt)) - 1) = 0`` con la identidad

        ``I = (n*Vt/R) * W( (Is*R/(n*Vt)) * exp((Vs + Is*R)/(n*Vt)) ) - Is``

    Devuelve ``(V_nodo, I_diodo)``.  Sirve como referencia exacta para validar
    el Newton numérico en el caso de un solo diodo.
    """
    m = modelo or ModeloDiodo()
    a = m.n_vt
    arg = (m.Is * R / a) * sp.exp(sp.Float((Vs + m.Is * R) / a))
    I = float((a / R) * sp.LambertW(arg) - m.Is)
    V = Vs - I * R
    return V, I
