import os
import re
from typing import List, Dict
from pathlib import Path
import sys

# ─────────────────────── detección de carpetas ───────────────────────
# Si se ejecuta como .exe (PyInstaller) → EXE_DIR = carpeta del ejecutable
# Si se ejecuta como script .py          → EXE_DIR = carpeta del fuente
if getattr(sys, "frozen", False):          # ← ejecutable .exe
    EXE_DIR = Path(sys.executable).resolve().parent
else:                                       # ← script .py
    EXE_DIR = Path(__file__).resolve().parent

# Punto de partida para resolver rutas relativas
BASE_DIR = EXE_DIR


def find_root(start: Path, target_subfolder: str = "datos", max_levels: int = 3) -> Path:
    """Sube directorios hasta encontrar la carpeta que contiene <target_subfolder>."""
    current = start
    for _ in range(max_levels + 1):
        if (current / target_subfolder).is_dir():
            return current
        current = current.parent
    raise FileNotFoundError(f"No se encontró la carpeta '{target_subfolder}' desde {start}")


# ───────── rutas relativas al proyecto (…/AutomatizacionWord) ─────────
ROOT_DIR       = find_root(BASE_DIR)
DATA_DIR       = ROOT_DIR / "datos"
PLANTILLAS_DIR = ROOT_DIR / "plantillas"
GENERADOS_DIR  = ROOT_DIR / "generados"
PLANOS_DIR     = ROOT_DIR / "planos"

DATOS_XLSX     = DATA_DIR / "plantilla_automatización.xlsx"
EXCEL_MT       = DATA_DIR / "CALCULOS ELECTRICOS MT-v00.xlsx"

# Ruta al PDF ESyS_Planos (relativa al proyecto)
ESYS_PDF = PLANTILLAS_DIR / "Doc6 - Estudios de Entidad Propia" / "AAAA-SOL-AL-PE-ESS-0001-A01_00 ESyS_Planos.pdf"

"""juntar_PDF – versión 8 (planos ordenados por prefijo 00‑01‑02…)
================================================================
Une en un único PDF los documentos generados previamente.
Sólo se han tocado rutas para que el script sea 100 % portátil.
"""

# ---------------------------------------------------------------------------
# PikePDF (stub para v<7)
# ---------------------------------------------------------------------------
try:
    from pikepdf import Pdf, PdfMerger
except ImportError:  # compatibilidad mínima
    from pikepdf import Pdf

    class PdfMerger:  # noqa: N801
        def __init__(self):
            self._out = Pdf.new()

        def append(self, pdf: Pdf):  # type: ignore[override]
            self._out.pages.extend(pdf.pages)

        def write(self, path: str):  # type: ignore[override]
            self._out.save(path)

        def close(self):  # type: ignore[override]
            self._out.close()

# ------------------------------------------------------
# Carpeta generada más reciente
# ------------------------------------------------------
carpetas = [d for d in os.listdir(GENERADOS_DIR) if (GENERADOS_DIR / d).is_dir()]
if not carpetas:
    raise RuntimeError("No se encontró ninguna carpeta en 'generados'. Ejecuta los pasos previos.")
ultima_carpeta = GENERADOS_DIR / max(carpetas, key=lambda d: (GENERADOS_DIR / d).stat().st_mtime)

# -----------------------------------------------------------
# Deducir sigla de planta desde Doc0 – Portada
# -----------------------------------------------------------
DOC0_DIR = ultima_carpeta / "Doc0 - Portada + Índice"
code_planta = None
for f in DOC0_DIR.iterdir():
    if f.suffix.lower() == ".pdf":
        m = re.match(r"^[0-9]{0,2}[_-]?([A-Z0-9]+)-SOL-", f.name, re.IGNORECASE)
        if m:
            code_planta = m.group(1).upper()
            break
if not code_planta:
    raise RuntimeError("No se pudo deducir la sigla de planta (GAMM, …)")
print(f"Carpeta: {ultima_carpeta}\nSigla: {code_planta}\n")

# ------------------------------------------------------------------
# Sufijos en orden corporativo (sin sigla) – documentos no‑planos
# ------------------------------------------------------------------
principal: List[str] = [
    "SOL-AL-PE-PRY-0001",
    "SOL-AL-PE-MEM-0001",
    "SOL-AL-PE-MEM-0001-A01",
    "SOL-AL-PE-MEM-0001-A02",
    "SOL-AL-PE-MEM-0001-A03",
    "SOL-AL-PE-MEM-0001-A04",
    "SOL-AL-PE-CAL-0001",
    "SOL-AL-PE-CAL-0001-A01",
    "SOL-AL-PE-CAL-0001-A02",
    "SOL-AL-PE-CAL-0001-A03",
    "SOL-AL-PE-CAL-0001-A04",
    "SOL-AL-PE-CAL-0001-A05",
    "SOL-AL-PE-CAL-0001-A05_00_MDT",
    "SOL-AL-PE-CAL-0001-A05_00_MDT_GAMM",
    "SOL-AL-PE-DWG-0001",
]

despues: List[str] = [
    "SOL-AL-PE-PCT-0001",
    "SOL-AL-PE-PRS-0001",
    "SOL-AL-PE-EST-0001",
    "SOL-AL-PE-ESS-0001-A01",  # después insertamos ESyS_Planos
    "SOL-AL-PE-EGR-0001-A02",
    "SOL-AL-PE-PDM-0001",
]

# -----------------------------------------------------------
# Regex helpers y deduplicación
# -----------------------------------------------------------
appended: set[str] = set()

SUF_PAT = r"(?:_[0-9]{2})?(?:[._-].*)?\.pdf$"


def _regex(sufijo: str) -> re.Pattern:
    return re.compile(fr"^[0-9]{{0,2}}[_-]?{code_planta}-{re.escape(sufijo)}{SUF_PAT}", re.IGNORECASE)


def _revision(path: str) -> int:
    m = re.search(r"_([0-9]{2})(?=[._-].*\.pdf$)", os.path.basename(path))
    return int(m.group(1)) if m else 0


def _latest(candidates: List[str]) -> str:
    return max(candidates, key=_revision)


def _buscar(suf: str, base: Path, recursive: bool = True) -> str | None:
    pat = _regex(suf)
    walker = os.walk(base) if recursive else [(base, [], os.listdir(base))]
    cands = [os.path.join(r, f) for r, _, fs in walker for f in fs if f.lower().endswith(".pdf") and pat.match(f)]
    return _latest(cands) if cands else None


def _append(path: str | None, merger):
    if not path:
        return False
    if path in appended:
        print("Omitido duplicado:", os.path.basename(path))
        return False
    merger.append(Pdf.open(path))
    appended.add(path)
    return True

# -----------------------------------------------------------
# Listar y ordenar planos por prefijo 00,01,…
# -----------------------------------------------------------
planos_paths: Dict[int, List[str]] = {}
planos_regex = re.compile(r"^([0-9]{2})[-_][A-Z0-9]+-SOL-.*-(DWG|DRW)-[0-9]{4}.*\.pdf$", re.IGNORECASE)
for root, _, files in os.walk(PLANOS_DIR):
    for f in files:
        if f.lower().endswith(".pdf") and code_planta in f:
            m = planos_regex.match(f)
            if m:
                idx = int(m.group(1))
                planos_paths.setdefault(idx, []).append(os.path.join(root, f))

planos_ordered: List[str] = []
for idx in sorted(planos_paths.keys()):
    planos_ordered.append(_latest(planos_paths[idx]))

# --------------------------
# Unir PDFs
# --------------------------
merger = PdfMerger()

for suf in principal:
    _append(_buscar(suf, ultima_carpeta), merger)

for p in planos_ordered:
    _append(p, merger)

for suf in despues:
    _append(_buscar(suf, ultima_carpeta), merger)
    if suf == "SOL-AL-PE-ESS-0001-A01":
        if not _append(str(ESYS_PDF) if ESYS_PDF.exists() else None, merger):
            print("No encontrado ESyS_Planos:", ESYS_PDF)

# ---------------------------
# Guardar PDF final
# ---------------------------
SALIDA = GENERADOS_DIR / f"{code_planta}-SOL-AL-PE-PRY-0001.pdf"
if SALIDA.exists():
    try:
        SALIDA.unlink()
    except PermissionError as e:
        raise RuntimeError(f"No se pudo sobrescribir {SALIDA}. Ciérralo en el visor y reintenta") from e

merger.write(str(SALIDA))
merger.close()
print("\nPDF final generado en:", SALIDA)
