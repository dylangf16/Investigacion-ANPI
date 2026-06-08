"""Descripción del circuito (netlist) y mapeo de incógnitas.

Un nodo puede ser de tres tipos:

* **Tierra** (``ground``): referencia a 0 V, no es incógnita.
* **Fuente DC**: voltaje fijo conocido, no es incógnita.  Alimenta el resto a
  través de las resistencias que se le conecten (igual que ``Vs`` en el
  ejemplo de ``docs/explicacion.md``).
* **Libre**: su voltaje es una incógnita del sistema ``F(x) = 0``.
"""
from __future__ import annotations

from .dispositivos import Diodo, Dispositivo, ModeloDiodo, Resistencia


class Red:
    """Contenedor del circuito: nodos, fuentes y dispositivos.

    Los nombres de nodo son cualquier objeto *hashable* (enteros, cadenas o
    tuplas ``(i, j)`` para una rejilla).  El orden de las incógnitas es el de
    primera aparición, para reproducibilidad.
    """

    def __init__(self, tierra: object = 0) -> None:
        self.tierra = tierra
        self.fuentes: dict[object, float] = {}
        self.dispositivos: list[Dispositivo] = []
        self._orden: list[object] = []
        self._vistos: set[object] = set()
        self.indice: dict[object, int] = {}
        self.n_incognitas: int = 0
        self._finalizada = False

    # -- registro de nodos -------------------------------------------------
    def _registrar(self, *nodos: object) -> None:
        for nodo in nodos:
            if nodo not in self._vistos:
                self._vistos.add(nodo)
                self._orden.append(nodo)
        self._finalizada = False

    # -- construcción del circuito ----------------------------------------
    def fuente_dc(self, nodo: object, voltaje: float) -> "Red":
        """Fija el voltaje de un nodo (fuente independiente a tierra)."""
        self._registrar(nodo)
        self.fuentes[nodo] = voltaje
        return self

    def resistencia(self, a: object, b: object, R: float) -> "Red":
        self._registrar(a, b)
        self.dispositivos.append(Resistencia(a, b, R))
        return self

    def diodo(self, anodo: object, catodo: object,
              modelo: ModeloDiodo | None = None) -> "Red":
        self._registrar(anodo, catodo)
        self.dispositivos.append(Diodo(anodo, catodo, modelo or ModeloDiodo()))
        return self

    def agregar(self, dispositivo: Dispositivo) -> "Red":
        """Agrega un dispositivo arbitrario que cumpla el protocolo."""
        self._registrar(*dispositivo.nodos())
        self.dispositivos.append(dispositivo)
        return self

    # -- finalización ------------------------------------------------------
    def finalizar(self) -> "Red":
        """Calcula el mapeo nodo libre -> índice de incógnita."""
        libres = [nodo for nodo in self._orden
                  if nodo != self.tierra and nodo not in self.fuentes]
        self.indice = {nodo: i for i, nodo in enumerate(libres)}
        self.n_incognitas = len(libres)
        self._finalizada = True
        return self

    @property
    def nodos_libres(self) -> list[object]:
        if not self._finalizada:
            self.finalizar()
        return list(self.indice.keys())

    def es_libre(self, nodo: object) -> bool:
        return nodo != self.tierra and nodo not in self.fuentes

    def __repr__(self) -> str:
        if not self._finalizada:
            self.finalizar()
        return (f"Red(incognitas={self.n_incognitas}, "
                f"fuentes={len(self.fuentes)}, "
                f"dispositivos={len(self.dispositivos)})")
