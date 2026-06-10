Esta implementación usa la variante con Armijo, que se considera más robusta, pero sigue siendo Newton Amortiguado


ENTRADA:
    F        : función residuo, x ↦ F(x) ∈ Rⁿ
    J        : Jacobiana, x ↦ J(x) ∈ Rⁿˣⁿ
    x₀       : aproximación inicial
    tol      : tolerancia (criterio de parada)
    max_iter : máximo de iteraciones
    α_min    : factor de amortiguamiento mínimo (p.ej. 1e-4)
    c        : constante de Armijo (p.ej. 1e-4)

SALIDA:
    x (solución), éxito (sí/no), nº de iteraciones

INICIO
    x ← x₀
    F_x ← F(x)
    norm ← ‖F_x‖₂

    PARA k = 0, 1, ..., max_iter-1 HACER

        # --- criterio de parada ---
        SI norm < tol ENTONCES
            DEVOLVER (x, éxito, k)          # convergió

        # --- paso de Newton: resolver J·Δx = −F  SIN invertir J ---
        Δx ← resolver_sistema_lineal( J(x), −F_x )

        # --- búsqueda de línea (backtracking + Armijo) ---
        α ← 1
        objetivo ← ½ · norm²
        aceptado ← falso
        MIENTRAS α ≥ α_min HACER
            x_nuevo ← x + α·Δx
            F_nuevo ← F(x_nuevo)
            SI F_nuevo es finito  Y  ½·‖F_nuevo‖₂² ≤ (1 − c·α)·objetivo ENTONCES
                aceptado ← verdadero
                SALIR DEL MIENTRAS
            SINO
                α ← α / 2               # frena el paso (amortiguamiento)

        SI NO aceptado ENTONCES
            DEVOLVER (x, fracaso, k)        # la línea de búsqueda no redujo ‖F‖

        # --- actualización ---
        x   ← x_nuevo
        F_x ← F_nuevo
        norm ← ‖F_x‖₂

    DEVOLVER (x, fracaso, max_iter)          # no convergió en max_iter
FIN
