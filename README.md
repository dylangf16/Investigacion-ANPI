# Investigacion-ANPI — Mini-SPICE de diodos

Punto de operación DC de circuitos no lineales con diodos mediante
**Newton-Raphson amortiguado**, resolviendo en cada paso una Jacobiana dispersa
**sin invertirla**. Es, en pequeño, lo que hace SPICE por dentro.

El planteamiento matemático, la motivación y las referencias están en
[docs/explicacion.md](docs/explicacion.md). En una frase: el análisis nodal con
la ecuación de Shockley produce un sistema no lineal `F(x) = 0`; la no
linealidad exponencial lo vuelve mal condicionado y hace que el Newton "puro"
del Avance 2 se desborde (overflow). El aporte es estabilizarlo (amortiguamiento
/ continuación) y resolver la Jacobiana dispersa en cada iteración.

## Estructura

```text
src/circuito/          Paquete principal
  constantes.py        Vt, Is, n, umbral de overflow de exp
  dispositivos.py      Resistencia y Diodo (Shockley) + stamping en F y J
  red.py               Netlist: nodos, fuentes, mapeo de incógnitas
  ensamblador.py       Construye F(x) y la Jacobiana dispersa (COO -> CSC)
  lineal.py            Solver lineal disperso sin invertir (LU dispersa / GMRES+ILU)
  newton.py            Newton puro, amortiguado (line search) y continuación
  simbolico.py         SymPy: verifica la Jacobiana y forma cerrada (Lambert W)
ejemplos/
  pequeno.py           Circuito de 4 nodos del documento
  exigente.py          Multiplicador de 6 nodos / 5 diodos con diodo de realimentación
  rejilla.py           Rejilla resistiva 2D con diodos (16 a miles de incógnitas)
experimentos/
  comparar_metodos.py  Tabla puro vs amortiguado vs continuación + condicionamiento
tests/                 Pruebas con pytest (incluye validación simbólica)
```

## Instalación

Crea y activa un entorno virtual, luego instala las dependencias.

**Windows (PowerShell):**

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
```

Si PowerShell bloquea la activación, habilítala una vez con
`Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`.

**Linux / macOS (bash):**

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Para salir del entorno: `deactivate`.

Dependencias: NumPy (vectores y residuo), SciPy (matrices dispersas y LU/GMRES),
SymPy (verificación de la Jacobiana y la solución cerrada de Lambert W),
Matplotlib (gráficas) y pytest.

## Uso rápido

```python
import sys; sys.path.insert(0, "src")
from circuito import Red, Ensamblador, newton_amortiguado

red = Red(tierra=0)
red.fuente_dc("Vs", 10.0)          # nodo de voltaje fijo
red.resistencia("Vs", "N1", 1e3)   # R en ohmios
red.diodo("N1", 0)                 # diodo ánodo -> cátodo

ens = Ensamblador(red)
r = newton_amortiguado(ens)
print(r)                                  # iteraciones, ||F||, tiempo
print(ens.voltajes_nodales(r.x))          # voltaje de cada nodo
```

Ejecutar los ejemplos y el experimento:

```bash
python ejemplos/pequeno.py
python ejemplos/rejilla.py
python experimentos/comparar_metodos.py
pytest -q
```

## Los tres métodos

| Función | Técnica numérica | Comportamiento |
| --- | --- | --- |
| `newton_puro` | Newton-Raphson del Avance 2 (paso completo) | Se desborda desde un arranque frío: la Jacobiana se vuelve singular (`inf`) |
| `newton_amortiguado` | Newton con búsqueda de línea (backtracking, Armijo) | Converge; evita el overflow frenando el paso |
| `newton_continuacion` | Source stepping / homotopía | Red de seguridad: rampa las fuentes de 0 a su valor real |

El sistema lineal de cada iteración se resuelve con `lineal.resolver` (LU
dispersa por defecto, GMRES+ILU para mallas muy grandes); nunca se calcula la
inversa de la Jacobiana.

## Resultado del experimento (resumen)

Newton puro **falla** (overflow → Jacobiana singular), mientras amortiguado
converge en pocas iteraciones y la continuación también; el número de condición
de la Jacobiana empeora ~90× entre los diodos apagados y el punto de operación.
La rejilla escala a miles de incógnitas (p. ej. 6400 nodos en <1 s) gracias a la
Jacobiana dispersa de 5 puntos.

## Notas

- El valor de `R2` (entre N1 y N2) no aparece en el documento; se fijó en 1 kΩ.
  Las ecuaciones `F` y la Jacobiana implementadas coinciden con las del documento.
- La verificación simbólica (`circuito.simbolico.verificar_jacobiana`) compara la
  Jacobiana estampada a mano contra la derivada exacta por SymPy en puntos
  aleatorios; es la garantía de que el stamping no tiene errores de signo.
