## Resumen

- El argumento de venta para un profe de análisis numérico no es el circuito, es esto: el sistema F(x)=0 tiene no linealidad exponencial, lo que lo vuelve mal condicionado y hace que el Newton-Raphson "puro" del Avance 2 falle por overflow. Resolverlo bien exige Newton amortiguado y resolver un sistema lineal disperso por iteración. Esa es la dificultad numérica real.
- Tres ganchos concretos y cuantificables: (1) las conductancias que entran a la Jacobiana varían ~11 ordenes de magnitud durante la iteracion, (2) un paso de Newton que empuje una unión por encima de ~18.4 V produce Inf en doble precisión, (3) la Jacobiana tiene la misma estructura dispersa que sus matrices de diferencias finitas, pero el sistema es no lineal.
- Conexión directa con el curso: es el método del Avance 2 P3, pero llevado a un régimen donde el método base se rompe, lo que da material genuino de "experimentos numéricos" (comparar Newton puro vs amortiguado vs continuación).

## El problema matemático (lo que el profe debe ver primero)

El análisis nodal con la ecuación de Shockley produce un sistema no lineal:

```
F(x) = 0,   x = (V1, V2, ..., VN) en R^N
```

donde cada ecuación es la ley de corrientes de Kirchhoff en un nodo, y la corriente de cada diodo es:

```
I = Is * ( exp( V / (n*Vt) ) - 1 )
```

| Parámetro | Valor típico | Unidad |
|---|---|---|
| Vt (voltaje térmico, 300 K) | 0.025852 | V |
| Is (corriente de saturación) | 1e-14 | A |
| n (factor de idealidad) | 1 a 2 | adimensional |

Por qué no lo resuelve matemática general: un solo diodo con una resistencia tiene forma cerrada vía la función W de Lambert (función especial, no elemental). Con dos o más diodos acoplados desaparece toda forma cerrada: hay varios exponenciales sumándose en distintas ecuaciones y no se pueden despejar. Queda obligatoriamente un sistema no lineal NxN, que es exactamente el objeto del Avance 2 P3.

## Drivers de complejidad numérica (cada uno atado al curso)

### 1. Mal condicionamiento por la exponencial

La derivada del diodo (la conductancia que va a la Jacobiana) es:

```
dI/dV = (Is / (n*Vt)) * exp( V / (n*Vt) )
```

Esto varía enormemente según el estado del diodo:

| Estado del diodo | Conductancia dI/dV | Derivación |
|---|---|---|
| Apagado (V cercano a 0) | ~3.9e-13 S | Is/(n*Vt) |
| Conduciendo a 5 mA | ~0.193 S | I/(n*Vt) |

Inferencia (derivada del modelo, no medida): los elementos diagonales de la Jacobiana abarcan cerca de 11 ordenes de magnitud durante la iteración. Eso eleva el número de condición y es lo que obliga a tratar el sistema lineal interno con cuidado. Conecta directo con la discusión de condicionamiento de los Avances 1 y 3.

### 2. El Newton puro del Avance 2 se rompe (overflow)

En doble precisión, `exp(z)` desborda a Inf cuando `z > 709.78`. Traducido a voltaje:

```
V_overflow = 709.78 * n * Vt = 709.78 * 0.025852 ≈ 18.35 V   (con n=1)
```

Un diodo real opera en ~0.7 V. Como la función es casi plana antes de la unión y explosivamente empinada después, un paso de Newton completo desde un arranque frío sobrepasa fácilmente ese umbral por un factor de decenas y devuelve Inf o NaN. Resultado: el método tal cual lo implementaron en el Avance 2 no converge. Esto no es un detalle, es el corazón del aporte.

### 3. La solución correcta requiere globalizar Newton

Lo que se necesita encima del Avance 2:

| Técnica | Qué es en términos numéricos | Rol |
|---|---|---|
| Amortiguamiento / limiting | Búsqueda de línea sobre el paso de Newton | Evita el overflow del punto 2 |
| Source stepping / gmin stepping | Continuación u homotopía | Lleva la solución desde un problema fácil al real |

Ambas son técnicas numéricas legítimas y explicables. Son justo lo que usa SPICE internamente. Aquí está el material de "experimentos numéricos" de la rúbrica: comparar Newton puro (falla), Newton amortiguado (converge) y, si quieren subir nivel, continuación.

### 4. Resolver el sistema lineal por iteración, sin invertir

El enunciado del Avance 2 P3 ya lo exige: no calcular la inversa de la Jacobiana, sino resolver el sistema lineal en cada paso con un solver eficiente. Aquí eso se vuelve obligatorio porque la Jacobiana es dispersa.

### 5. La Jacobiana tiene estructura (gancho fuerte para este profe)

Si arman el circuito como una malla de resistencias 2D (rejilla n x n) con diodos en algunas ramas, la parte lineal de la Jacobiana es el Laplaciano del grafo, que es la misma estructura dispersa de 5 puntos que las matrices de diferencias finitas de los Avances 1, 3 y 5. El mensaje para el profe:

> "Profesor, la Jacobiana de cada iteración tiene la misma estructura pentadiagonal por bloques que las matrices de diferencias finitas que vimos, pero ahora el sistema es no lineal y hay que rearmar y refactorizar esa matriz en cada paso de Newton."

## Circuito concreto para fijar el tamaño

Para que sea claramente mayor a 3x3 y no descalificable como trivial:

| Topología | Nodos / incógnitas | Comentario |
|---|---|---|
| Rejilla resistiva 2D, malla 4x4 interior, con diodos en ramas | 16 | Jacobiana dispersa tipo 5 puntos |
| Multiplicador de voltaje (escalera de diodos), 8 a 12 etapas | 8 a 12 | Aplicación CE real: fuentes, RFID, energy harvesting |
| Puente rectificador + carga + filtro resistivo multi-nodo | 6 a 10 | Clásico y vistoso para el video |

Recomendación: la rejilla 2D con diodos da el sistema más grande, la estructura dispersa más reconocible y la mejor narrativa numérica. Tamaño sugerido: malla de al menos 4x4 (16 incógnitas) para no quedar corto.

## Cómo presentarlo al profe en una frase

"El punto de operación de una red con varios diodos es un sistema no lineal F(x)=0 sin solución cerrada. Aplicamos Newton-Raphson para sistemas (Avance 2 P3), pero la no linealidad exponencial vuelve el problema mal condicionado y hace que el Newton estándar desborde numéricamente. El aporte es estabilizarlo con amortiguamiento y resolver en cada iteración un sistema lineal disperso con la estructura de las matrices de diferencias finitas del curso."

## Referencias para la parte de investigación

| Elemento | Referencia | Qué citas |
|---|---|---|
| Problema general (circuito) | A. Sedra, K. Smith, Microelectronic Circuits | Modelo del diodo (Shockley) y análisis nodal |
| Por qué numérico / sin forma cerrada | L. O. Chua, P. M. Lin, Computer-Aided Analysis of Electronic Circuits | Formulación de redes no lineales y solución por Newton |
| Globalización de Newton | L. W. Nagel, "SPICE2" (1975, UC Berkeley memo ERL-M520) | Damping, gmin y source stepping en simuladores reales |
| Método numérico de clase | Cualquier texto de análisis numérico que use en CE-1111 | Newton-Raphson para sistemas |

Aquí va, en cristiano:

## Qué haríamos

Calcular el punto de operación (los voltajes en cada nodo) de un circuito con varios diodos. Es lo que hace un simulador tipo SPICE por dentro. No se puede a mano ni con mate de cole: hay que resolverlo con métodos numéricos.

## El circuito

Una rejilla de resistencias (malla tipo 4x4) con diodos metidos en varias ramas. Eso nos da 16 nodos, o sea 16 incógnitas (16 voltajes que no conocemos).

Cada diodo se comporta con la ecuación de Shockley:

```
I = Is * (exp(V / (n*Vt)) - 1)
```

Esa exponencial es la que hace todo el problema difícil.

## Qué usamos del curso

- Newton-Raphson para sistemas (Avance 2, Pregunta 3). El mismo que ya programamos.
- En cada iteración hay que resolver un sistema lineal grande sin invertir la matriz (lo que pide ese mismo avance).
- La matriz que sale (la Jacobiana) es dispersa, con la misma estructura de las matrices de diferencias finitas que vimos en los Avances 1, 3 y 5.

O sea: reusamos código que ya tenemos, pero lo llevamos a un terreno donde se rompe y hay que arreglarlo. Ahí está la nota.

## Por qué NO es sencillo (esto es lo que vendemos)

1. No tiene solución cerrada. Con un solo diodo existe una fórmula rara (función de Lambert). Con varios diodos acoplados, nada: queda obligatoriamente un sistema no lineal de 16x16.

2. El Newton normal se revienta. La exponencial crece tan rápido que un solo paso de Newton manda el voltaje a un valor enorme, y `exp(...)` se desborda a infinito en el computador (pasa arriba de ~18 V). El método tal cual lo hicimos en el Avance 2 NO converge.

3. Hay que estabilizarlo. Le agregamos amortiguamiento (frenar el paso de Newton para que no se dispare). Eso es justo lo que hace SPICE de verdad.

4. La matriz está mal condicionada. Las conductancias que entran a la Jacobiana cambian como 11 órdenes de magnitud entre un diodo apagado y uno encendido. Eso vuelve el sistema lineal numéricamente delicado.

## El experimento para el paper

Comparar tres corridas: Newton puro (falla), Newton amortiguado (converge), y opcional continuación. Tabla con iteraciones, error y tiempo. Eso llena la sección de experimentos numéricos que pide la rúbrica.

## La frase para el profe

"Es Newton-Raphson para sistemas del Avance 2, pero la no linealidad exponencial lo vuelve mal condicionado y hace que el método estándar desborde. El aporte es estabilizarlo y resolver una Jacobiana dispersa en cada paso."

Resumen de una línea para los amigos: armamos un mini-SPICE, usamos el Newton que ya programamos, pero el circuito lo hace explotar y nosotros lo arreglamos. Eso es lo difícil y lo que da nota.

Buena pregunta, y te voy a ser directo: hay referencias que te puedo dar con alta confianza y algunas que debes verificar antes de poner en el paper. Te digo cuáles son cuáles.

## Referencias con alta confianza (puedes citar directamente)

### Modelo del diodo (Shockley)
**W. Shockley**, "The Theory of p-n Junctions in Semiconductors and p-n Junction Transistors," *Bell System Technical Journal*, vol. 28, no. 3, pp. 435-489, July 1949.

- Premio Nobel de Física 1956. La ecuación que usamos viene de aquí. No hay referencia más autoritativa.

### SPICE (el simulador real que hace exactamente nuestro problema)
**L. W. Nagel y D. O. Pederson**, "SPICE (Simulation Program with Integrated Circuit Emphasis)," Technical Report UCB/ERL M382, EECS Department, University of California, Berkeley, April 1973.

**L. W. Nagel**, "SPICE2: A Computer Program to Simulate Semiconductor Circuits," Technical Report UCB/ERL M520, EECS Department, University of California, Berkeley, May 1975.

- Estos son reportes técnicos de Berkeley, accesibles online. El segundo describe explícitamente el source stepping y el gmin stepping, que es la globalización de Newton que nosotros implementamos.

### Formulación nodal del sistema (por qué queda un sistema NxN)
**C. W. Ho, A. E. Ruehli, y P. A. Brennan**, "The Modified Nodal Approach to Network Analysis," *IEEE Transactions on Circuits and Systems*, vol. 22, no. 6, pp. 504-509, June 1975.

- Paper clásico de IEEE que formaliza exactamente cómo KCL + modelos de dispositivos producen el sistema F(x)=0 que nosotros resolvemos. Muy citado en EDA (Electronic Design Automation).

### Libro de circuitos (para el modelo del diodo en contexto CE)
**A. S. Sedra y K. C. Smith**, *Microelectronic Circuits*, Oxford University Press, 7ma edición, 2014.

- El libro de electrónica que usa TEC en cursos de CE. Le da contexto al profe de qué es el circuito.

### Libro de análisis numérico de circuitos
**L. O. Chua y P. M. Lin**, *Computer-Aided Analysis of Electronic Circuits: Algorithms and Computational Techniques*, Prentice-Hall, 1975.

- El libro canónico que conecta análisis de circuitos con métodos numéricos. Tiene la formulación del Newton-Raphson para redes no lineales.

## Referencia para la función W de Lambert (por qué un solo diodo sí tiene solución cerrada y varios no)

**T. C. Banwell y A. Jaski**, "Exact analytical solution for current flow through diode with series resistance," *Electronics Letters*, vol. 36, no. 4, pp. 291-292, Feb. 2000.

- Muestra que un diodo con una resistencia se resuelve con la función W de Lambert. Esto es oro para el paper: lo citas para argumentar que la solución cerrada existe en el caso trivial y desaparece cuando hay varios diodos acoplados.
- **Verifica el volumen y páginas** en IEEE Xplore antes de citar. El título y autores son correctos.

## Referencia para el método numérico del curso

**W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery**, *Numerical Recipes: The Art of Scientific Computing*, 3ra edición, Cambridge University Press, 2007.

- Cubre Newton-Raphson para sistemas, condicionamiento y backtracking (el amortiguamiento). Disponible online en numericipes.com. Te sirve para la sección de métodos numéricos del paper.

## Cómo usar estas referencias en el paper IEEE

| Sección del paper | Qué citar |
|---|---|
| Introducción | Nagel & Pederson 1973 (motivación: SPICE resuelve este problema) |
| Problema general | Sedra & Smith (modelo del circuito y diodo) |
| Problema matemático | Shockley 1949 (ecuación del diodo), Ho et al. 1975 (por qué queda un sistema NxN) |
| Por qué numérico | Banwell & Jaski 2000 (solución cerrada existe solo para 1 diodo, aquí no aplica) |
| Método numérico | Chua & Lin 1975 (Newton para redes), Nagel 1975 (damping/source stepping) |
| Implementación | Numerical Recipes (Newton amortiguado) |


## Lo que le dices al profe, en orden

### 1. El problema de ingeniería

"El problema es calcular el punto de operación de un circuito con múltiples diodos: encontrar los voltajes en cada nodo cuando el circuito está en estado estable. Es exactamente lo que hace un simulador SPICE internamente."

### 2. El modelo matemático

"Cada diodo sigue la ecuación de Shockley:

```
I = Is * (exp(V / (n*Vt)) - 1)
```

Al aplicar KCL en cada nodo del circuito, esa ecuación aparece dentro de cada balance de corrientes. El resultado es un sistema no lineal:

```
F(x) = 0,    x en R^N
```

donde cada incógnita es el voltaje en un nodo y N es la cantidad de nodos del circuito."

### 3. Por qué no lo resuelve matemática general

"Con un solo diodo y una resistencia existe una solución cerrada mediante la función W de Lambert (Banwell & Jaski, 2000). En cuanto hay dos o más diodos acoplados en distintos nodos, los exponenciales aparecen en ecuaciones distintas y no se pueden despejar. No existe ninguna forma analítica. La única salida es un método numérico."

### 4. Qué método del curso usamos y por qué es difícil

"Aplicamos Newton-Raphson para sistemas no lineales, que es exactamente el Avance 2 Pregunta 3. Pero aquí el método base se rompe por dos razones numéricas concretas:

La primera es overflow. La derivada del diodo que entra a la Jacobiana es:

```
dI/dV = (Is / (n*Vt)) * exp(V / (n*Vt))
```

Un paso de Newton completo desde condiciones iniciales frías puede empujar el voltaje más allá de ~18 V, y `exp(709)` desborda a infinito en doble precisión. El método estándar no converge.

La segunda es mal condicionamiento. Esa misma derivada varía aproximadamente 11 órdenes de magnitud entre un diodo apagado y uno encendido. Los elementos de la Jacobiana cambian drásticamente entre iteraciones, lo que hace el sistema lineal interno numéricamente delicado."

### 5. El aporte del proyecto

"La solución es agregar amortiguamiento al paso de Newton: en lugar de tomar el paso completo, lo escalamos hasta que la norma de F disminuya. Eso estabiliza el método. El experimento numérico del paper compara tres variantes: Newton puro (falla), Newton amortiguado (converge) y, opcionalmente, continuación por parámetro. La tabla de resultados reporta iteraciones, error y tiempo por método."

### 6. Conexión con el curso (esto es lo que más le importa)

"La Jacobiana que resolvemos en cada iteración tiene la misma estructura dispersa de las matrices de diferencias finitas del Avance 1 y el Avance 3. En cada iteración resolvemos ese sistema lineal sin invertir la matriz, exactamente como lo exige el Avance 2 Pregunta 3. El proyecto toma el método del curso y lo lleva a un régimen donde el método base falla, lo que obliga a estabilizarlo. Esa es la dificultad numérica real."

### 7. Las referencias

"Las referencias clave son:

- Shockley (1949) en Bell System Technical Journal para la ecuación del diodo.
- Ho, Ruehli y Brennan (1975) en IEEE Transactions on Circuits and Systems para la formulación nodal que produce el sistema F(x)=0.
- Nagel (1975), reporte técnico UCB/ERL M520 de Berkeley, que describe exactamente el amortiguamiento y la continuación que nosotros implementamos.
- Banwell y Jaski (2000) en Electronics Letters para argumentar que la solución cerrada existe solo en el caso trivial de un diodo."

## El circuito

```
         R1          D1
Vs ----[1kΩ]---- N1 --|>|-- N3
                  |          |
                 [D2]       [R4=2kΩ]
                  |          |
                  N2        GND
                  |    [D3]
                [R3=1kΩ] |>|
                  |          |
                  N3 --------N4
                             |
                           [R5=1kΩ]
                             |
                            GND
```

Mas claro:

```
Vs=10V ---[R1=1k]--- N1 ---[D1]--- N3
                     |              |
                    [D2]          [R4=2k]
                     |              |
                     N2            GND
                     |
                   [R3=1k]
                     |
                     N3 ---[D3]--- N4
                                   |
                                 [R5=1k]
                                   |
                                  GND
```

4 nodos desconocidos: N1, N2, N3, N4.

## El sistema F(x) = 0

Aplicando KCL en cada nodo, con la ecuación de Shockley en cada diodo (Is = 1e-14 A, Vt = 0.02585 V):

**KCL en N1:**
```
F1 = (N1-Vs)/R1 + (N1-N2)/R2 + Is*(exp((N1-N3)/Vt) - 1) = 0
```

**KCL en N2:**
```
F2 = (N2-N1)/R2 + (N2-N3)/R3 + Is*(exp((N2-N4)/Vt) - 1) = 0
```

**KCL en N3:**
```
F3 = (N3-N2)/R3 + N3/R4
     - Is*(exp((N1-N3)/Vt) - 1)
     + Is*(exp((N3-N4)/Vt) - 1) = 0
```

**KCL en N4:**
```
F4 = N4/R5
     - Is*(exp((N2-N4)/Vt) - 1)
     - Is*(exp((N3-N4)/Vt) - 1) = 0
```

## La Jacobiana (lo que se resuelve en cada iteración de Newton)

Definiendo las conductancias dinámicas de cada diodo:

```
gD1 = (Is/Vt) * exp((N1-N3)/Vt)
gD2 = (Is/Vt) * exp((N2-N4)/Vt)
gD3 = (Is/Vt) * exp((N3-N4)/Vt)
```

La Jacobiana queda:

```
     | N1                | N2        | N3                    | N4        |
-----|-------------------|-----------|-----------------------|-----------|
F1   | 1/R1+1/R2+gD1     | -1/R2     | -gD1                  | 0         |
F2   | -1/R2             | 1/R2+1/R3+gD2 | -1/R3            | -gD2      |
F3   | -gD1              | -1/R3     | 1/R3+gD1+gD3+1/R4    | -gD3      |
F4   | 0                 | -gD2      | -gD3                  | 1/R5+gD2+gD3 |
```

## Por qué es complejo y no trivial

**1. Los gD cambian radicalmente entre iteraciones:**

| Diodo | Apagado (V~0) | Encendido (V~0.7V) |
|---|---|---|
| gD | ~3.9e-13 S | ~0.19 S |

Eso es 11 órdenes de magnitud de diferencia. La Jacobiana cambia completamente de una iteración a la siguiente.

**2. El Newton puro explota.** Si en alguna iteración N1-N3 sube a 1V, el término `exp(1/0.02585)` = `exp(38.7)` = 6.2e16. Si sube a 0.5 más, desborda a Inf.

**3. El sistema no es lineal en ningún sentido.** No se puede simplificar, no se puede separar, no tiene forma cerrada.

## Por eso Newton-Raphson modificado

En cada iteración, en vez de:
```
x = x + dx
```

Se hace:
```
alpha = 1.0
while ||F(x + alpha*dx)|| > ||F(x)||:
    alpha = alpha / 2
x = x + alpha * dx
```

## Tema puntual del proyecto

**Título sugerido:**

> "Análisis del Punto de Operación de Circuitos No Lineales con Diodos mediante Newton-Raphson Amortiguado"

En inglés para el paper IEEE (que lo pide en ese formato):

> "DC Operating Point Analysis of Nonlinear Diode Circuits Using Damped Newton-Raphson"

Ese es el término exacto que usa la literatura: **DC Operating Point** o **DC Operating Point Analysis**. Si el profe busca eso en IEEE Xplore, encuentra cientos de papers.

## Método usado (nombre formal)

**Newton-Raphson con búsqueda de línea (Damped Newton-Raphson / Newton-Raphson with backtracking line search)**

No es solo "Newton-Raphson modificado", el nombre preciso importa para las referencias.

## Referencias exactas

### Referencia 1 - De donde sale el sistema F(x) = 0

Esta es la más importante. Es donde se formaliza que aplicar KCL con modelos de dispositivos produce un sistema no lineal NxN:

> C. W. Ho, A. E. Ruehli, y P. A. Brennan, "The Modified Nodal Approach to Network Analysis," **IEEE Transactions on Circuits and Systems**, vol. 22, no. 6, pp. 504-509, Jun. 1975.

- Búscalo en **IEEE Xplore** con ese título exacto, lo encontrás de inmediato.
- En ese paper está explícitamente la formulación nodal que produce el sistema que nosotros resolvemos.
- **Alta confianza** en volumen, número y páginas.

### Referencia 2 - De donde sale la ecuación de Shockley

> W. Shockley, "The Theory of p-n Junctions in Semiconductors and p-n Junction Transistors," **Bell System Technical Journal**, vol. 28, no. 3, pp. 435-489, Jul. 1949.

- Disponible online en los archivos de Bell Labs / Nokia Bell Labs.
- Es la fuente original de la ecuación del diodo. No hay nada más autoritativo.
- **Alta confianza.**

### Referencia 3 - De donde sale el Newton amortiguado para circuitos

> L. W. Nagel, "SPICE2: A Computer Program to Simulate Semiconductor Circuits," Technical Report **UCB/ERL M520**, EECS Department, University of California, Berkeley, May 1975.

- Disponible gratis en: **https://www2.eecs.berkeley.edu/Pubs/TechRpts/1975/9602.html**
- En ese reporte está descrito explícitamente el damping y el source stepping que nosotros implementamos.
- **Alta confianza.** Es un reporte técnico, no una revista, pero es la referencia canónica del campo.

### Referencia 4 - Por qué no hay solución cerrada con varios diodos

> T. C. Banwell y A. Jaski, "Exact analytical solution for current flow through diode with series resistance," **Electronics Letters**, vol. 36, no. 4, pp. 291-292, Feb. 2000.

- Búscalo en **IEEE Xplore** con el título exacto.
- Muestra que un solo diodo con una resistencia tiene forma cerrada via función W de Lambert. Con varios diodos acoplados eso desaparece.
- **Verifica páginas en IEEEXplore antes de citar**, el título y autores son correctos.

### Referencia 5 - Para el método numérico del curso

> W. H. Press, S. A. Teukolsky, W. T. Vetterling, B. P. Flannery, **Numerical Recipes: The Art of Scientific Computing**, 3ra ed., Cambridge University Press, 2007.

- Capítulo 9: "Root Finding and Nonlinear Sets of Equations"
- Sección 9.7: "Newton-Raphson Method for Nonlinear Systems of Equations"
- Sección 9.7 también cubre el backtracking. Disponible en **numerical.recipes**
- **Alta confianza.**

## Cómo citarlas en el paper

| Sección del paper | Referencia | Qué dice exactamente |
|---|---|---|
| Problema matemático | Ho et al. (1975) | "La formulación de análisis nodal modificado [1] produce el sistema no lineal F(x)=0" |
| Ecuación del diodo | Shockley (1949) | "La corriente del diodo sigue el modelo de Shockley [2]" |
| Por qué no hay solución cerrada | Banwell & Jaski (2000) | "Para un solo diodo existe solución via función W de Lambert [3]; con múltiples diodos acoplados no existe forma analítica" |
| Newton amortiguado | Nagel (1975) | "El amortiguamiento del paso de Newton es la estrategia implementada en SPICE [4]" |
| Método numérico | Numerical Recipes | "La búsqueda de línea se implementó conforme a [5]" |
