#!/usr/bin/env python3
"""
dwg_PDF.py – Publica en PDF todos los *layouts* de cada DWG situado en la
carpeta ``dwg/`` usando AutoCAD Core Console (accoreconsole.exe).

v2025‑07‑30 c – Corrección «Archivo no válido»
------------------------------------------------
* El script .SCR ya **no envía la ruta entre comillas** (AutoCAD ES la trata
  como texto literal → «Archivo no válido»).
* El .DSD se guarda con finales de línea *CRLF* (`\r\n`) para seguir la
  convención de Windows.
* Mantiene comandos para silenciar diálogos (`FILEDIA 0`, etc.).
* Constantes DEFAULT_* permiten ejecutar sin argumentos.

Uso rápido
~~~~~~~~~~
```
python dwg_PDF.py
```
Si quieres sobreescribir rutas:
```
python dwg_PDF.py --dwg-dir X:\planos --output-dir Z:\pdf --accore C:\Ruta\accoreconsole.exe
```
"""
from __future__ import annotations

import argparse
import subprocess
import tempfile
import uuid
from pathlib import Path
from typing import List

try:
    import win32com.client  # type: ignore
except ImportError:
    win32com = None  # type: ignore

# -----------------------------------------------------------------------------
# 1)  RUTAS POR DEFECTO --------------------------------------------------------
# -----------------------------------------------------------------------------
ROOT_DIR = Path(__file__).resolve().parent
DEFAULT_ACCORE = Path(r"C:\Program Files\Autodesk\AutoCAD 2026\accoreconsole.exe")
DEFAULT_DWG_DIR = ROOT_DIR / "dwg"
DEFAULT_OUTPUT_DIR = ROOT_DIR / "planos"
DEFAULT_PLOTTER = "DWG To PDF.pc3"
DEFAULT_PAGE_SETUP = ""

# -----------------------------------------------------------------------------
# 2)  HELPERS ------------------------------------------------------------------
# -----------------------------------------------------------------------------

def list_layouts(dwg_path: Path) -> List[str]:
    """Devuelve los layouts (excepto Model) usando COM; vacío si falla."""
    if win32com is None:
        return []
    try:
        acad = win32com.client.Dispatch("AutoCAD.Application")  # type: ignore
        acad.Visible = False
        acad.DisplayAlerts = False
        doc = acad.Documents.Open(str(dwg_path))
        layouts = [lay.Name for lay in doc.Layouts if not lay.ModelType]
        doc.Close(False)
        return layouts
    except Exception:
        return []


def build_dsd(dwg_paths: List[Path], out_dir: Path, plotter: str, page_setup: str) -> Path:
    """Genera un .DSD válido (CRLF)."""

    def add_sheet(buf: List[str], idx: int, dwg: Path, layout: str) -> None:
        buf.extend([
            f"[Sheet{idx}]",
            f"DWG={dwg}",
            f"Layout={layout}",
            ("Setup=" if not page_setup else f"Setup={page_setup}"),
            f"OriginalSheetPath={dwg.parent}",
            "HasPlotSupport=TRUE",
            "Include=TRUE",
        ])

    lines: List[str] = ["[DSD]", "Version=1", "SheetType=0"]
    idx = 0
    for dwg in dwg_paths:
        layouts = list_layouts(dwg)
        if layouts:
            for lay in layouts:
                idx += 1
                add_sheet(lines, idx, dwg, lay)
        else:
            idx += 1
            add_sheet(lines, idx, dwg, "")

    lines.extend([
        "[Options]",
        "Type=6",
        "OutputType=1",
        "PromptForFileName=FALSE",
        f"Device={plotter}",
        f"OutputLocation={out_dir}",
    ])

    dsd_path = Path(tempfile.gettempdir()) / f"publish_{uuid.uuid4().hex}.dsd"
    # unir con CRLF para máxima compatibilidad
    dsd_path.write_bytes("\r\n".join(line.encode() for line in lines) + b"\r\n")
    return dsd_path


def build_scr(dsd_path: Path) -> Path:
    """Scr minimalista sin comillas y sin líneas vacías."""
    scr_lines = [
        "FILEDIA 0",
        "CMDDIA 0",
        "PROXYNOTICE 0",
        "PROXYSHOW 0",
        "SECURELOAD 0",
        "_-PUBLISH",
        str(dsd_path),  # ← sin comillas
        "y",
        "_QUIT",
    ]
    scr_path = Path(tempfile.gettempdir()) / f"publish_{uuid.uuid4().hex}.scr"
    scr_path.write_text("\n".join(scr_lines), encoding="utf-8")
    return scr_path


def run_publish(accore: Path, scr: Path) -> None:
    subprocess.run([str(accore), "/s", str(scr)], check=True)

# -----------------------------------------------------------------------------
# 3)  CLI ----------------------------------------------------------------------
# -----------------------------------------------------------------------------

def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Batch DWG → PDF (todos los layouts)")
    p.add_argument("--accore", type=Path, default=DEFAULT_ACCORE)
    p.add_argument("--dwg-dir", type=Path, default=DEFAULT_DWG_DIR)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    p.add_argument("--plotter", default=DEFAULT_PLOTTER)
    p.add_argument("--page-setup", default=DEFAULT_PAGE_SETUP)
    return p.parse_args()

# -----------------------------------------------------------------------------
# 4)  MAIN ---------------------------------------------------------------------
# -----------------------------------------------------------------------------

def main() -> None:
    args = parse_cli()

    if not args.accore.is_file():
        raise SystemExit(f"❌ No existe accoreconsole en {args.accore}")

    dwg_files = sorted(args.dwg_dir.glob("*.dwg"))
    if not dwg_files:
        raise SystemExit(f"❌ No se encontraron DWG en {args.dwg_dir}")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    dsd = build_dsd(dwg_files, args.output_dir, args.plotter, args.page_setup)
    scr = build_scr(dsd)

    print(f"Publicando {len(dwg_files)} DWG (todos los layouts)…")
    try:
        run_publish(args.accore, scr)
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"❌ Core Console terminó con error {e.returncode}")

    print("✓ PDF generados en:", args.output_dir)


if __name__ == "__main__":
    main()
