#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Generador masivo de .docx + PDF optimizado (solo imagen de tabla MT)
────────────────────────────────────────────────────────────────────────
• Rellena todas las plantillas .docx (subcarpetas incl.) con datos/plantilla automatización.xlsx
• Inserta imagen del rango MT (Excel) donde se ponga {{tabla_mt}} en la plantilla.
• Replica jerarquía en generados\
"""

from datetime import datetime
from decimal import Decimal, ROUND_HALF_UP
import unicodedata, re, os, pathlib, sys, shutil, tempfile
import pandas as pd
from docxtpl import DocxTemplate, InlineImage
from jinja2 import Environment, Undefined, TemplateSyntaxError
from pathlib import Path, PurePath
from docx.shared import Cm

# ─── localizar la raíz real del proyecto ─────────────────────────────
BASE_DIR = Path(sys.executable if getattr(sys, "frozen", False) else __file__).resolve().parent

def find_root(start: Path, target_subfolder: str = "datos", max_levels: int = 3) -> Path:
    current = start
    for _ in range(max_levels + 1):
        if (current / target_subfolder).is_dir():
            return current
        current = current.parent
    raise FileNotFoundError(f"No se encontró la carpeta '{target_subfolder}' desde {start}")

# ─── rutas ──────────────────────────────────────────────────────────
ROOT_DIR       = find_root(BASE_DIR)
DATA_DIR       = ROOT_DIR / "datos"
PLANTILLAS_DIR = ROOT_DIR / "plantillas"
GENERADOS_DIR  = ROOT_DIR / "generados"
PLANOS_DIR     = ROOT_DIR / "planos"

DATOS_XLSX     = DATA_DIR / "plantilla_automatización.xlsx"
EXCEL_MT       = DATA_DIR / "CALCULOS ELECTRICOS MT-v00.xlsx"

# ─────────── Excel → PNG ────────────

def exportar_tabla_mt_png(excel: str, hoja: str, rng: str) -> str:
    """Exporta la tabla a un PNG temporal y retorna la ruta."""
    import win32com.client as win32
    import tempfile
    xlPicture = -4147
    xl = win32.DispatchEx("Excel.Application")
    xl.Visible = False
    tmp_png = tempfile.NamedTemporaryFile(delete=False, suffix=".png").name
    try:
        wb = xl.Workbooks.Open(excel, ReadOnly=True)
        ws = wb.Sheets(hoja)
        ws.Activate()
        ws.Range(rng).CopyPicture(Format=xlPicture)
        width = ws.Range(rng).Width
        height = ws.Range(rng).Height
        chart = ws.ChartObjects().Add(10, 10, width, height)
        chart.Activate()
        chart.Chart.Paste()
        chart.Chart.Export(tmp_png, "PNG")
        chart.Delete()
    finally:
        wb.Close(False)
        xl.Quit()
    return tmp_png

# ─────────── utilidades ────────────

def slug(txt: str) -> str:
    txt = unicodedata.normalize("NFD", txt).encode("ascii", "ignore").decode()
    txt = re.sub(r"[^A-Za-z0-9_-]+", "_", txt)
    return txt.strip("_")

def _fmt_con_dec(v, *, simbolo=False):
    try:
        d = Decimal(str(v)).quantize(Decimal("0.01"), ROUND_HALF_UP)
        miles, cents = f"{d:,.2f}".split(".")
        miles = miles.replace(",", ".")
        txt = miles if cents == "00" else (f"{miles},{cents[0]}" if cents[1]=="0" else f"{miles},{cents}")
        return f"{txt} €" if simbolo else txt
    except Exception:
        return str(v)

def _fmt_sin_dec(v, *, simbolo=False):
    try:
        d = Decimal(str(v)).quantize(Decimal("1"), ROUND_HALF_UP)
        txt = f"{d:,.0f}".replace(",", ".")
        return f"{txt} €" if simbolo else txt
    except Exception:
        return str(v)

def formato1(v): return _fmt_con_dec(v, simbolo=False)
def formato(v):  return _fmt_sin_dec(v, simbolo=False)
def euros1(v):   return _fmt_con_dec(v, simbolo=True)
def euros(v):    return _fmt_sin_dec(v, simbolo=True)
def fecha(val):
    if isinstance(val, datetime):
        return val.strftime("%d/%m/%Y")
    dt = pd.to_datetime(val, errors="coerce", dayfirst=True)
    return dt.strftime("%d/%m/%Y") if pd.notnull(dt) else str(val)

# ─────────── procesar plantilla ────────────

def procesar_plantilla(path_docx: str, ctx: dict, destino_base: Path, img_mt_path: str = None):
    tpl = DocxTemplate(path_docx)
    # Solo insertar la imagen si la plantilla es de MT y hay imagen generada
    if img_mt_path and 'MT' in Path(path_docx).stem.upper():
        ctx = ctx.copy()  # por si se reusa el dict
        ctx["tabla_mt"] = InlineImage(tpl, img_mt_path, width=Cm(16))
    env = Environment(undefined=Undefined)
    env.filters.update({
        "formato": formato, "formato1": formato1,
        "euros": euros, "euros1": euros1,
        "fecha": fecha,
    })
    try:
        tpl.render(ctx, jinja_env=env)
    except TemplateSyntaxError as err:
        print(f"⚠️  Marcador mal escrito en {Path(path_docx).name} → {err}")
        return
    rel = Path(path_docx).relative_to(PLANTILLAS_DIR)
    destino = destino_base / rel.parent
    destino.mkdir(parents=True, exist_ok=True)
    base = Path(path_docx).stem.replace("AAAA", ctx["siglas_planta"]).replace("_MDT", " MDT")
    tpl.save(destino / f"{slug(base)}.docx")

# ─────────── main ────────────

def main():
    if not DATOS_XLSX.exists():
        raise FileNotFoundError(f"No existe {DATOS_XLSX}")

    df = pd.read_excel(DATOS_XLSX, engine="openpyxl")
    if df.empty:
        raise SystemExit("El Excel está vacío.")

    plantillas = [str(p) for p in pathlib.Path(PLANTILLAS_DIR).rglob("*.docx")
                  if not p.name.startswith("~$")]
    if not plantillas:
        raise SystemExit("No hay .docx en la carpeta plantillas.")

    for _, fila in df.iterrows():
        if pd.isna(fila.iloc[1]):
            continue
        ctx = {str(k).strip().lower(): v for k, v in fila.dropna().items()}
        for k, v in list(ctx.items()):
            if isinstance(v, (pd.Timestamp, datetime)):
                ctx[k] = v.strftime("%d/%m/%Y")
        nom_planta = ctx.get("nombre_planta") or "SIN_NOMBRE"
        siglas = ctx.get("codigo_planta") or ctx.get("siglas_planta")
        if not siglas:
            siglas = slug(nom_planta)[:4]
        siglas = slug(str(siglas)).upper()
        ctx["siglas_planta"] = siglas
        carpeta = os.path.join(
            GENERADOS_DIR,
            f"{slug(nom_planta)}_{datetime.now():%Y-%m-%d_%H%M%S}"
        )
        os.makedirs(carpeta, exist_ok=True)
        print(f"\n— Procesando: {nom_planta}  |  Plantillas: {len(plantillas)}")
        # Si alguna plantilla contiene '{{tabla_mt}}', genera la imagen solo una vez y pon su ruta en el contexto
        if any("{{tabla_mt}}" in open(p, encoding="utf-8", errors="ignore").read() for p in plantillas):
            calc_dir = os.path.join(carpeta, "Doc2 - Cálculos")
            os.makedirs(calc_dir, exist_ok=True)
            img_mt_path = os.path.join(calc_dir, "tabla_mt.png")
            try:
                excel_range_to_png(
                    excel=str(EXCEL_MT),
                    hoja="CALCULOS ELECTRICOS MT-v00",
                    rng="E26:AB76",
                    destino_png=img_mt_path
                )
                print("Imagen tabla_mt generada:", os.path.exists(img_mt_path), img_mt_path)
                ctx["tabla_mt_path"] = img_mt_path
            except Exception as e:
                print("⚠️  No se pudo exportar tabla MT como imagen:", e)
        for tpl in plantillas:
            print("  →", os.path.basename(tpl))
            procesar_plantilla(tpl, ctx, carpeta)
        print("  ✓ Documentos y Word listos en:", carpeta)
    print("\n  Proceso completado.")

# ──────────────────────────────────────────────
if __name__ == "__main__":
    main()
