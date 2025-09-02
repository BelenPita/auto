#!/usr/bin/env python
"""generar_pdfs.py ─ Paso intermedio del flujo documental
Convierte los *.docx* a PDF en la MISMA carpeta que el Word.
"""

from __future__ import annotations
import subprocess, sys
from pathlib import Path
from typing import List

# ───────── localizar raíz proyecto (soporta .exe) ─────────
if getattr(sys, "frozen", False):          # ← ejecutable .exe
    EXE_DIR  = Path(sys.executable).resolve().parent
else:                                      # ← script .py
    EXE_DIR  = Path(__file__).resolve().parent

# BASE_DIR será la carpeta donde está el .exe o el .py.
BASE_DIR = EXE_DIR


def find_root(start: Path, target_subfolder: str = "datos", max_levels: int = 3) -> Path:
    """Sube hasta max_levels para hallar la carpeta que contiene <target_subfolder>."""
    current = start
    for _ in range(max_levels + 1):
        if (current / target_subfolder).is_dir():
            return current
        current = current.parent
    raise FileNotFoundError(f"No se encontró la carpeta '{target_subfolder}' desde {start}")

# ───────── rutas relativas al proyecto ─────────
ROOT_DIR       = find_root(BASE_DIR)      # …/AutomatizacionWord
DATA_DIR       = ROOT_DIR / "datos"
PLANTILLAS_DIR = ROOT_DIR / "plantillas"
GENERADOS_DIR  = ROOT_DIR / "generados"
PLANOS_DIR     = ROOT_DIR / "planos"

DATOS_XLSX     = DATA_DIR / "plantilla_automatización.xlsx"
EXCEL_MT       = DATA_DIR / "CALCULOS ELECTRICOS MT-v00.xlsx"

# ───────── back-ends de conversión ─────────
try:
    from docx2pdf import convert as docx2pdf_convert  # type: ignore
except ImportError:
    docx2pdf_convert = None


def convert_with_docx2pdf(docx: Path) -> None:
    docx2pdf_convert(str(docx))                 # crea foo.pdf junto a foo.docx


def convert_with_libreoffice(docx: Path) -> None:
    subprocess.run(
        [
            "soffice", "--headless", "--convert-to", "pdf",
            "--outdir", str(docx.parent), str(docx)
        ],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    )


# Selección de backend --------------------------------------------------------
if sys.platform.startswith(("win", "darwin")) and docx2pdf_convert:
    convert = convert_with_docx2pdf
else:
    convert = convert_with_libreoffice


# ───────── utilidades ─────────

def collect_docx(root: Path) -> List[Path]:
    return [p for p in root.rglob("*.docx") if not p.name.startswith("~$")]


# ───────── main ─────────

def main() -> None:
    if not GENERADOS_DIR.exists():
        print("[ERROR] La carpeta 'generados' no existe. Ejecuta primero generar_docs.py")
        sys.exit(1)

    # ── busca en TODOS los proyectos de 'generados' ──
    project_dirs = [d for d in GENERADOS_DIR.iterdir() if d.is_dir()]
    docx_files: List[Path] = []
    for d in project_dirs:
        docx_files.extend(collect_docx(d))

    if not docx_files:
        print("No se encontraron archivos .docx para convertir.")
        return

    print(f"Convirtiendo DOCX→PDF en {len(project_dirs)} proyectos…\n")

    total = len(docx_files)
    for idx, docx in enumerate(docx_files, 1):
        pdf = docx.with_suffix(".pdf")
        if pdf.exists():
            print(f"[{idx}/{total}] ✔︎ Ya existe {pdf.relative_to(docx.parents[2])} (omitido)")
            continue
        try:
            convert(docx)
            status = "✓"
        except Exception as e:
            status = "✗"
            print(f"[{idx}/{total}] {status} Error al convertir {docx.name}: {e}")
            continue
        print(f"[{idx}/{total}] {status} {pdf.relative_to(docx.parents[2])}")

    print("\nConversión terminada ✔️")


if __name__ == "__main__":
    main()