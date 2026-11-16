import os
import sys

sys.dont_write_bytecode = True
PYTHON_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MATRIX_DIR = os.path.abspath(os.path.join(PYTHON_DIR, ".."))
if PYTHON_DIR not in sys.path:
    sys.path.insert(0, PYTHON_DIR)
