"""Configuración de pytest: pone ``src`` y ``ejemplos`` en el path."""
import pathlib
import sys

RAIZ = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RAIZ / "src"))
sys.path.insert(0, str(RAIZ / "ejemplos"))
