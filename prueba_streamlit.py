# -*- coding: utf-8 -*-
import math
import numpy as np
from math import pi
from math import sqrt as raiz
from math import log
import cmath
import io
import sys
import ast
import os
import unicodedata
import streamlit as st
import pandas as pd




print = st.write

buffer = io.StringIO()
sys.stdout = buffer

# Coeficientes
mu0 = 4e-7 * pi # H/m
epsilon0 = 8.854e-12  # F/m



#-----------------------------------------------------------------------------------------------------------------------------------------------
# LECTURA DE VALORES EXCEL
#-----------------------------------------------------------------------------------------------------------------------------------------------

def _normalize_key(s):
    if s is None:
        return ""
    s = str(s)
    s = unicodedata.normalize('NFKD', s)
    s = ''.join(c for c in s if not unicodedata.combining(c))
    s = ''.join(ch for ch in s if ch.isalnum())
    return s.lower()

def load_parameters_from_excel(file_name="Cálculos eléctricos.xlsx"):
    # Lee un archivo de Excel con dos columnas (nombre, valor)
    base_dir = os.path.dirname(__file__)
    candidates = [file_name, file_name + ".xlsx", file_name + ".xls"]
    path = None
    for c in candidates:
        p = os.path.join(base_dir, c)
        if os.path.exists(p):
            path = p
            break
    if path is None:
        return {}

    mapping = {
        "tensionnominal": "tension_nominal",
        "tension_nominal": "tension_nominal",
        "tensionmaselevada": "tensionmaselevada",
        "tension_mas_elevada": "tensionmaselevada",
        "ndecircuitos": "ndecircuitos",
        "n_de_circuitos": "ndecircuitos",
        "ndeconductoresporfase": "ndeconductoresporfase",
        "n_de_conductores_por_fase": "ndeconductoresporfase",
        "coeficientetemperatura": "coef_temp",
        "coeficiente_temperatura": "coef_temp",
        "coeficienttemperatura": "coef_temp",
        "frecuencia": "frecuencia",
        "longitud": "longitud",
        "resistividad": "resistividad",
        "cosphi": "cos_phi",
        "cos_phi": "cos_phi",
        "factordepotencia": "cos_phi",
        "factor_de_potencia": "cos_phi",
        "potenciaatransportar": "potencia_transportada",
        "potencia_a_transportar": "potencia_transportada",
        "tiempoaccionamientoproteccion": "tiempo_accionamiento_proteccion",
        "altitudmedia": "altitud_media",
        "tipodeconductor": "conductor",
        "tipo_de_conductor": "conductor",
        "material": "material",
        "temperaturainvierno": "temperatura_invierno",
        "temperaturaverano": "temperatura_verano",
        "temperaturaconductor": "Tc",
        "emisividadconductor": "emisividad_conductor",
        "emisividad_conductor": "emisividad_conductor",
        "coeficienteabsorcion": "coeficiente_absorcion",
        "coeficiente_absorcion": "coeficiente_absorcion",
        "radiacioninvierno": "radiacion_invierno",
        "radiacionverano": "radiacion_verano",
        "velocidadviento": "velocidad_viento",
        "mc": "mc",
        "mtinvierno": "mt_invierno",
        "mtverano": "mt_verano",
        "resistenciatierra": "resistencia_tierra",
        "secciontierra": "seccion_tierra",
        "diametrotierra": "diametro_tierra",
        "coeficientetemperaturatierra": "coef_temp_tierra",
        "coeficiente_temperatura_tierra": "coef_temp_tierra",
        "coeftemptierra": "coef_temp_tierra",
        "uno": "uno",
        "dos": "dos",
        "tres": "tres",
        "cuatro": "cuatro",
        "cinco": "cinco",
        "seis": "seis",
        "tierra": "tierra",
        "distanciaconductores": "distancia_conductores"
    }

    def _parse_value(v):
        if v is None:
            return None
        if isinstance(v, (int, float, complex)):
            return v
        if isinstance(v, str):
            s = v.strip()
            try:
                val = ast.literal_eval(s)
                return val
            except Exception:
                try:
                    if "," in s and "." not in s:
                        s = s.replace(',', '.')
                    return float(s)
                except Exception:
                    return s
        return v

    result = {}
    try:
        import pandas as pd
        df = pd.read_excel(path, header=None, usecols=[0, 1])
        for idx, row in df.iterrows():
            key = row.iloc[0]
            val = row.iloc[1]
            nk = _normalize_key(key)
            if nk in mapping:
                varname = mapping[nk]
                parsed = _parse_value(val)
                result[varname] = parsed
    except Exception:
        try:
            from openpyxl import load_workbook
            wb = load_workbook(path, data_only=True)
            ws = wb.active
            for row in ws.iter_rows(min_row=1, max_col=2, values_only=True):
                key, val = row[0], row[1]
                nk = _normalize_key(key)
                if nk in mapping:
                    varname = mapping[nk]
                    parsed = _parse_value(val)
                    result[varname] = parsed
        except Exception:
            return {}
    return result


# Sobreescribir valores desde el Excel "Cálculos eléctricos" si existe
excel_params = load_parameters_from_excel("Cálculos eléctricos.xlsx")
if not excel_params:
    try:
        excel_params = load_parameters_from_excel("Cálculos eléctricos")
    except Exception:
        excel_params = {}

# Asignar valores leídos del Excel, sobrescribiendo los por defecto
for varname, value in excel_params.items():
    globals()[varname] = value
    # También asignar al scope local para que el resto del código use estos valores
    if varname == 'coef_temp':
        coef_temp = value
    elif varname == 'frecuencia':
        frecuencia = value
    elif varname == 'longitud':
        longitud = value
    elif varname == 'resistividad':
        resistividad = value
    elif varname == 'tension_nominal':
        tension_nominal = value
    elif varname == 'tensionmaselevada':
        tensionmaselevada = value
    elif varname == 'ndecircuitos':
        ndecircuitos = value
    elif varname == 'ndeconductoresporfase':
        ndeconductoresporfase = value
    elif varname == 'conductor':
        conductor = value
    elif varname == 'cos_phi':
        cos_phi = value
    elif varname == 'potencia_transportada':
        potencia_transportada = value
    elif varname == 'tiempo_accionamiento_proteccion':
        tiempo_accionamiento_proteccion = value
    elif varname == 'altitud_media':
        altitud_media = value
    elif varname == 'conductor':
        conductor = value
    elif varname == 'material':
        material = value
    elif varname == 'temperatura_invierno':
        temperatura_invierno = value
    elif varname == 'temperatura_verano':
        temperatura_verano = value
    elif varname == 'Tc':
        Tc = value
    elif varname == 'emisividad_conductor':
        emisividad_conductor = value
    elif varname == 'coeficiente_absorcion':
        coeficiente_absorcion = value
    elif varname == 'radiacion_invierno':
        radiacion_invierno = value
    elif varname == 'radiacion_verano':
        radiacion_verano = value
    elif varname == 'velocidad_viento':
        velocidad_viento = value
    elif varname == 'mc':
        mc = value
    elif varname == 'mt_invierno':
        mt_invierno = value
    elif varname == 'mt_verano':
        mt_verano = value
    elif varname == 'resistencia_tierra':
        resistencia_tierra = value
    elif varname == 'seccion_tierra':
        seccion_tierra = value
    elif varname == 'diametro_tierra':
        diametro_tierra = value
    elif varname == 'coef_temp_tierra':
        coef_temp_tierra = value
    elif varname == 'uno':
        uno = value
    elif varname == 'dos':
        dos = value
    elif varname == 'tres':
        tres = value
    elif varname == 'cuatro':
        cuatro = value
    elif varname == 'cinco':
        cinco = value
    elif varname == 'seis':
        seis = value
    elif varname == 'tierra':
        tierra = value
    elif varname == 'distancia_conductores':
        distancia_conductores = value

# Asegurar que las variables de puntos sean (x, y)
def _ensure_point(p):
    try:
        if isinstance(p, (list, tuple)) and len(p) >= 2:
            return (float(p[0]), float(p[1]))
        if isinstance(p, str):
            s = p.strip().replace(';', ',')
            try:
                v = ast.literal_eval(s)
                if isinstance(v, (list, tuple)) and len(v) >= 2:
                    return (float(v[0]), float(v[1]))
            except Exception:
                parts = s.split(',')
                if len(parts) >= 2:
                    return (float(parts[0]), float(parts[1]))
    except Exception:
        pass
    return p

uno = _ensure_point(uno)
dos = _ensure_point(dos)
tres = _ensure_point(tres)
cuatro = _ensure_point(cuatro)
cinco = _ensure_point(cinco)
seis = _ensure_point(seis)
tierra = _ensure_point(tierra)
puntos = [uno, dos, tres, cuatro, cinco, seis, tierra]
def _xy(pt):
    try:
        if isinstance(pt, (list, tuple)) and len(pt) >= 2:
            return float(pt[0]), float(pt[1])
        if isinstance(pt, str):
            s = pt.strip().replace(';', ',')
            try:
                v = ast.literal_eval(s)
                if isinstance(v, (list, tuple)) and len(v) >= 2:
                    return float(v[0]), float(v[1])
            except Exception:
                parts = s.split(',')
                if len(parts) >= 2:
                    return float(parts[0]), float(parts[1])
        if isinstance(pt, (int, float)):
            return float(pt), 0.0
    except Exception:
        pass
    raise ValueError(f"Formato de punto inválido: {pt}")

uno_x, uno_y = _xy(uno)
dos_x, dos_y = _xy(dos)
tres_x, tres_y = _xy(tres)
cuatro_x, cuatro_y = _xy(cuatro)
cinco_x, cinco_y = _xy(cinco)
seis_x, seis_y = _xy(seis)
tierra_x, tierra_y = _xy(tierra)


print("\nDATOS")
print("\nDatos línea")
print(f"Tensión nominal: {tension_nominal} kV")
print(f"Frecuencia: {frecuencia} Hz")
print(f"Tensión más elevada: {tensionmaselevada} kV")
print(f"Número de circuitos: {ndecircuitos}")
print(f"Número de conductores por fase: {ndeconductoresporfase}")
if ndeconductoresporfase==2:
    print(f"Distancia entre conductores: {distancia_conductores}")
print(f"Tipo de conductor: {conductor}")
print(f"Longitud: {longitud} km")
print(f"Altitud media: {altitud_media}")
print(f"Factor de potencia: {cos_phi}")
print(f"Potencia a transportar: {potencia_transportada} MVA")
print(f"Coeficiente de temperatura: {coef_temp} 1/ºC")

print(f"Tiempo accionamiento protección: {tiempo_accionamiento_proteccion} s")
print(f"Resistividad terreno: {resistividad} Ω·m")
print(f"Temperatura conductor: {Tc} ºC")

print(f"\nEmisividad conductor: {emisividad_conductor}")
print(f"Coeficiente absorción: {coeficiente_absorcion}")
print(f"Velocidad viento: {velocidad_viento} m/s")

print("\nDatos conductor tierra")
print(f"Sección tierra: {seccion_tierra} mm²")
print(f"Diámetro tierra: {diametro_tierra} mm")
print(f"Resistencia tierra: {resistencia_tierra}")
print(f"Coeficiente teperatura tierra: {coef_temp_tierra:2e}")



print("\nUbicación conductores")
if ndecircuitos==1:
    print(f"Punto 1: {uno}")
    print(f"Punto 2: {dos}")
    print(f"Punto 3: {tres}")
    print(f"Punto tierra: {cuatro}")
elif ndecircuitos==2:
    print(f"Punto 1: {uno}")
    print(f"Punto 2: {dos}")
    print(f"Punto 3: {tres}")
    print(f"Punto 4: {cuatro}")    
    print(f"Punto 5: {cinco}")  
    print(f"Punto 6: {seis}")  
    print(f"Punto tierra: {tierra}")  


print("\nDatos localización")
print(f"Temperatura invierno: {temperatura_invierno} ºC")
print(f"Temperatura verano: {temperatura_verano} ºC")
print(f"Radiación invierno: {radiacion_invierno} W/m²")
print(f"Radiación verano: {radiacion_verano} W/m²")


# Comparación los valores de diametro y resistencia de la parte superior con los valores de la normativa a partir del dato de conductor
diametro_normativa = {"LA-30": 7.14,"LA-56": 9.45,"LA-78": 11.34,"LA-110": 14.0,"LA-145": 15.75,"LA-180": 17.5,"LA-280": 21.8,"LA-380": 25.38,"LA-455": 27.72,"LA-545": 30.42,"LA-635": 32.85}
resistencia_normativa = {"LA-30": 1.0749, "LA-56": 0.6136, "LA-78": 0.4261,"LA-110": 0.3066,"LA-145": 0.2422,"LA-180": 0.1962,"LA-280": 0.1194,"LA-380": 0.0857,"LA-455": 0.0718,"LA-545": 0.0596,"LA-635": 0.0511}
diametro_alambre_ext_normativa = {"LA-30": 2.38,"LA-56": 3.15,"LA-78": 3.78,"LA-110": 2,"LA-145": 2.25,"LA-180": 2.5,"LA-280": 2.68,"LA-380": 2.82,"LA-455": 3.08,"LA-545": 3.38,"LA-635": 3.65}
seccion_normativa = {"LA-30": 31.1,"LA-56": 54.6,"LA-78": 78.6,"LA-110": 116.2,"LA-145": 147.1,"LA-180": 181.6,"LA-280": 281.1,"LA-380": 381,"LA-455": 454.5,"LA-545": 547.3,"LA-635": 636.6}
composicion_normativa = {"LA-30": "6+1","LA-56": "6+1","LA-78": "6+1","LA-110": "30+7","LA-145": "30+7","LA-180": "30+7","LA-280": "26+7","LA-380": "54+7","LA-455": "54+7","LA-545": "54+7","LA-635": "54+19"}
print(f"\nVALORES SEGÚN EL CONDUCTOR {conductor}")
diametro_calculado = diametro_normativa.get(conductor, None)
resistencia_calculada = resistencia_normativa.get(conductor, None)
diametro_alambre_ext_calculado = diametro_alambre_ext_normativa.get(conductor, None)
seccion_calculada = seccion_normativa.get(conductor, None)
composicion_calculada = composicion_normativa.get(conductor, None)
print(f"Diámetro: {diametro_calculado} mm")
print(f"Resistencia: {resistencia_calculada} Ω/km")
print(f"Diámetro alambre: {diametro_alambre_ext_calculado} mm")
print(f"Sección: {seccion_calculada} mm²")
print(f"Composición: {composicion_calculada}")

diametro = diametro_calculado
seccion = seccion_calculada
diametro_alambre_ext = diametro_alambre_ext_calculado
composicion = composicion_calculada
resistencia = resistencia_calculada
n_i = ndeconductoresporfase











#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# SIMPLE CIRCUITO SIMPLEX
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
if ndecircuitos==1 and ndeconductoresporfase==1:      
        
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CÁLCULO RESISTENCIA
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Definir función para cálculo de resistencia a temperatura
    def resistencia_a_temp(resistencia_20C, temp):
        return resistencia_20C * (1 + coef_temp * (temp - 20))
    # Cálculo de la resistencia a 85ºC y en ca
    st.header("\n1. RESISTENCIA ELÉCTRICA DE LA LÍNEA")
    st.markdown("---")
    resistencia_85C = resistencia_a_temp(resistencia, Tc)
    print(f"Resistencia a {Tc:.2f}ºC: {resistencia_85C:.6f} Ω/km")
    reactancia_pelicular_85C = 10**-7 * 8 * pi * frecuencia * (1 / (resistencia_a_temp(resistencia, Tc) / 1000))
    print(f"Reactancia por efecto pelicular a {Tc:.2f}ºC: {reactancia_pelicular_85C:.6f} Ω/km")
    ys = reactancia_pelicular_85C ** 2 / (192 + 0.8*reactancia_pelicular_85C**2)
    print(f"Factor de efecto pelicular ys: {ys:.6f}")
    resistencia_ca = resistencia_85C * (1 + ys)
    print(f"Resistencia en CA a 85ºC: {resistencia_ca:.6f} Ω/km")
    r_ca_longitud = resistencia_ca * longitud
    print(f"Resistencia en CA para {longitud} km: {r_ca_longitud:.6f} Ω")

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CÁLCULO MATRIZ IMPEDANCIAS
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # Cálculo matriz de impedancias
    st.header("\n2. MATRIZ DE IMPEDANCIAS")
    st.markdown("---")
    st.write("""
    El cálculo de la matriz de impedancias se realiza mediante la teoría de Carson. 
    De este modo se calcula la impedancia propia y mutua de todos los conductores 
    que forman la línea teniendo en cuenta el terreno.
    """)
    st.write("Impedancia propia de un conductor:")
    st.latex(r"""
    Z_{i,i} = R_i
    + \frac{j \mu_0}{2 \pi} \, \omega 
    \left(
    \frac{\mu_r}{4 n_i'} + 
    \ln \left( \frac{D'_{i,i}}{r_{eq,i}} \right)
    \right)
    + \frac{\mu_0}{\pi} \, \omega \left( P_{i,i} + j Q_{i,i} \right)
    """)
    st.write("Impedancia mutua entre dos conductores:")
    st.latex(r"""
    Z_{i,j} = 
    \frac{j \mu_0}{2 \pi} \, \omega 
    \ln \left( \frac{D'_{i,j}}{D_{i,j}} \right)
    + \frac{\mu_0}{\pi} \, \omega \left( P_{i,j} + j Q_{i,j} \right)
    """)
    st.write("""
    Donde:
    - $R_i$ es la resistencia eléctrica del conductor $i$  
    - $D_{i,j}$ es la distancia entre el conductor $i$ y el $j$  
    - $D'_{i,j}$ es la distancia entre el conductor $i$ y la imagen del conductor $j$ respecto al suelo  
    - $r_{eq,i}$ es el radio equivalente del conductor $i$  
    - $\omega$ es la pulsación $2 \pi f$  
    - $\mu_0$ es la permeabilidad del vacío  
    - $\mu_r$ es la permeabilidad relativa del conductor  

    Los valores de $P$ y $Q$ son una serie infinita de términos que describen 
    el comportamiento del terreno. Para cálculos a frecuencias industriales, los 
    términos de la serie a partir del segundo pueden ser despreciados. Por tanto:
    """)
    st.latex(r"""
    P_{i,j} = \frac{\pi}{8} - \frac{k_{i,j} \cdot \cos(\theta_{i,j})}{3 \sqrt{2}}
    """)
    st.latex(r"""
    Q_{i,j} = \frac{1}{2} \cdot \ln \left( \frac{1.85138}{k_{i,j}} \right)
    + \frac{k_{i,j} \cdot \cos(\theta_{i,j})}{3 \sqrt{2}}
    """)
    penetracion_terreno = raiz(resistividad/(pi*frecuencia*mu0))
    print(f"Penetración terreno: {penetracion_terreno:.4f} m")
    # Matriz de distancias
    puntos = [uno, dos, tres, cuatro]
    n_puntos = len(puntos)
    matriz_distancias = [[0]*n_puntos for _ in range(n_puntos)]
    def _xy(pt):
        try:
            if isinstance(pt, (list, tuple)) and len(pt) >= 2:
                return float(pt[0]), float(pt[1])
            if isinstance(pt, str):
                s = pt.strip().replace(';', ',')
                try:
                    v = ast.literal_eval(s)
                    if isinstance(v, (list, tuple)) and len(v) >= 2:
                        return float(v[0]), float(v[1])
                except Exception:
                    parts = s.split(',')
                    if len(parts) >= 2:
                        return float(parts[0]), float(parts[1])
            if isinstance(pt, (int, float)):
                return float(pt), 0.0
        except Exception:
            pass
        raise ValueError(f"Formato de punto inválido: {pt}")

    for i in range(n_puntos):
        for j in range(n_puntos):
            x1, y1 = _xy(puntos[i])
            x2, y2 = _xy(puntos[j])
            distancia = raiz((x2-x1)**2 + (y2-y1)**2)
            matriz_distancias[i][j] = distancia
    print("\nMatriz de distancias (km):")
    df=pd.DataFrame(matriz_distancias).round(4)
    latex_matrix = r"D = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)

    # Matriz D_prima: distancia entre puntos y sus espejos respecto al suelo
    matriz_D_prima = [[0]*n_puntos for _ in range(n_puntos)]
    for i in range(n_puntos):
        for j in range(n_puntos):
            x1, y1 = _xy(puntos[i])
            # Espejo del punto j respecto al suelo
            x2, y2 = _xy(puntos[j])
            x2_espejo, y2_espejo = x2, -y2
            distancia = raiz((x2_espejo-x1)**2 + (y2_espejo-y1)**2)
            matriz_D_prima[i][j] = distancia
    print("\nMatriz D' (distancias entre puntos y espejos) (km):")
    df=pd.DataFrame(matriz_D_prima).round(4)
    latex_matrix = r"D' = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # kij = raiz(2)*D_prima_ij / penetracion_terreno
    matriz_kij = [[0]*n_puntos for _ in range(n_puntos)]
    for i in range(n_puntos):
        for j in range(n_puntos):
            matriz_kij[i][j] = raiz(2) * matriz_D_prima[i][j] / penetracion_terreno  
    st.write("\nMatriz k_ij:")
    df=pd.DataFrame(matriz_kij).round(4)
    latex_matrix = r"k_{ij} = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # P_ij = pi/8 - k_ij*cos(tetha_ij)/(3*raiz(2)) siendo cos(tetha_ij) = (h_i + h_j )/ D_prima_ij
    matriz_Pij = [[0]*n_puntos for _ in range(n_puntos)]
    for i in range(n_puntos):
        for j in range(n_puntos):
            h_i = _xy(puntos[i])[1]
            h_j = _xy(puntos[j])[1]
            D_prima_ij = matriz_D_prima[i][j]
            if D_prima_ij != 0:
                cos_tetha_ij = (h_i + h_j) / D_prima_ij
            else:
                cos_tetha_ij = 0
            k_ij = matriz_kij[i][j]
            P_ij = (pi / 8) - (k_ij * cos_tetha_ij) / (3 * raiz(2))
            matriz_Pij[i][j] = P_ij
    print("\nMatriz P_ij:")
    df=pd.DataFrame(matriz_Pij).round(4)
    latex_matrix = r"P_{ij} = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # Q_ij = 0.5*log neperiano(1.85138/k_ij) + k_ij*cos(tetha_ij)/(3*raiz(2))
    matriz_Qij = [[0]*n_puntos for _ in range(n_puntos)]
    for i in range(n_puntos):
        for j in range(n_puntos):
            k_ij = matriz_kij[i][j]
            h_i = _xy(puntos[i])[1]
            h_j = _xy(puntos[j])[1]
            D_prima_ij = matriz_D_prima[i][j]
            if D_prima_ij != 0:
                cos_tetha_ij = (h_i + h_j) / D_prima_ij
            else:
                cos_tetha_ij = 0
            if k_ij != 0:
                Q_ij = 0.5 * log(1.85138 / k_ij) + (k_ij * cos_tetha_ij) / (3 * raiz(2))           
            else:
                Q_ij = 0
            matriz_Qij[i][j] = Q_ij
    print("\nMatriz Q_ij:")
    df=pd.DataFrame(matriz_Qij).round(4)
    latex_matrix = r"Q_{ij} = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # Resistencia de cada uno de los conductores (ultimo a tierra con dato de resistencia a tierra directamente)
    resistencias_conductores = [resistencia_ca, resistencia_ca, resistencia_ca, resistencia_tierra]
    print("\nResistencias de los conductores:")
    df=pd.DataFrame(resistencias_conductores).round(4)
    latex_matrix = r"R = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    latex_matrix += r"\text{Ω/km}"
    st.latex(latex_matrix)
    # Radio de cada conductor en mm
    radio_conductores_mm = [
        diametro / 2,
        diametro / 2,
        diametro / 2,
        diametro_tierra / 2
    ]
    print("\nRadio de los conductores (mm):")
    df=pd.DataFrame(radio_conductores_mm).round(4)
    latex_matrix = r"Radio = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    latex_matrix += r"\text{mm}"
    st.latex(latex_matrix)
    # Matriz de impedancias con:
    # Z_ii = R_i + j*mu0*frecuencia*((1/(4*n_i)+log(D_prima_ii/r_i)))+mu0*2*frecuencia*(P_ii + j*Q_ii) 
    # Z_ij = j*mu0*frecuencia*(log(D_prima_ij/D_ij)) + mu0*2*frecuencia*(P_ij + j*Q_ij) para i != j
    matriz_impedancias = [[0]*n_puntos for _ in range(n_puntos)]
    for i in range(n_puntos):
        for j in range(n_puntos):
            R_i = resistencias_conductores[i]
            D_ij = matriz_distancias[i][j]
            D_prima_ij = matriz_D_prima[i][j]
            P_ij = matriz_Pij[i][j]
            Q_ij = matriz_Qij[i][j]
            r_i = radio_conductores_mm[i] / 1000  # Convertir mm a m
            if i == j:
                n_i = 1  # Número de conductores en fase
                Z_ii = R_i + (1j * mu0 * frecuencia * ((1 / (4 * n_i)) + log(D_prima_ij / r_i)) + mu0 * 2 * frecuencia * (P_ij + 1j * Q_ij))*1000
                matriz_impedancias[i][j] = Z_ii
            else:
                Z_ij = (1j * mu0 * frecuencia * (log(D_prima_ij / D_ij)) + mu0 * 2 * frecuencia * (P_ij + 1j * Q_ij))*1000
                matriz_impedancias[i][j] = Z_ij
    print("\nMatriz de impedancias (Ω/km):")
    df=pd.DataFrame(matriz_impedancias).round(4)
    latex_matrix = r"Z = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # Como la linea posee cable de tierra, es necesario realizar un análisis matricial para eliminarlos y obtener una matriz 3*3 que representa las impedancias por fase
    # Z=[Zf, Zft; Ztf, Zt]
    # Zfas = Zf - Zft * inv(Zt) * Ztf

    Z = np.array(matriz_impedancias)
    Zf = Z[0:3, 0:3]
    Zt = Z[3:4, 3:4]
    Zft = Z[0:3, 3:4]
    Ztf = Z[3:4, 0:3]
    Zt_inv = np.linalg.inv(Zt)
    Zfas = Zf - Zft @ Zt_inv @ Ztf
    print("\nMatriz de impedancias por fase (Ω/km):")
    df=pd.DataFrame(Zfas).round(4)
    latex_matrix = r"Z_{fas} = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # Finalmente, multiplicar por la longitud de la línea para obtener las impedancias totales
    Zfas_total = Zfas * longitud
    print("\nMatriz de impedancias por fase total para la longitud dada (Ω):")
    df=pd.DataFrame(Zfas_total).round(4)
    latex_matrix = r"Z_{fas} = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # IMPEDANCIAS DE SECUENCIA
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    st.header("\n3. IMPEDANCIAS DE SECUENCIA")
    st.markdown("---")
    A = np.array([[1, 1, 1],
                [complex(-0.5, -raiz(3)/2), complex(-0.5, raiz(3)/2), 1],
                [complex(-0.5, raiz(3)/2), complex(-0.5, -raiz(3)/2), 1]])
    A_inv = np.linalg.inv(A)
    Z_seq = A_inv @ Zfas @ A
    print("\nMatriz de impedancias de secuencia (Ω):")
    df=pd.DataFrame(Z_seq).round(4)
    latex_matrix = r"Z_{012} = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # Matriz anterior de impedancias de secuencia pero con argumento y angulo
    print("\nMatriz de impedancias de secuencia con magnitud y ángulo (Ω):")
    latex_matrix = r"Z_{012} = \begin{bmatrix}"
    for fila in Z_seq:
        elementos = []
        for Z in fila:
            magnitud = abs(Z)
            angulo = math.degrees(math.atan2(Z.imag, Z.real))
            elementos.append(
                rf"{magnitud:.4f} \angle {angulo:.2f}^\circ"
            )
        
        latex_matrix += " & ".join(elementos)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # Impedancia homopolar de la linea (Z0)
    Z0 = Z_seq[2, 2]
    print(f"\nImpedancia homopolar de la línea (Z0): {Z0:.4f} Ω")
    # Impedancia directa e inversa de la linea (Z1)
    Z1 = Z_seq[1, 1]
    print(f"Impedancia directa e inversa de la línea (Z1): {Z1:.4f} Ω")
    #Teniendo en cuenta la longitud de la línea
    Z_seq_total = Z_seq * longitud
    print("\nMatriz de impedancias de secuencia total para la longitud dada (Ω):")
    df=pd.DataFrame(Z_seq_total).round(4)
    latex_matrix = r"Z_{012} = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # Matriz anterior de impedancias de secuencia pero con argumento y angulo
    print("\nMatriz de impedancias de secuencia total con magnitud y ángulo para la longitud dada (Ω):")
    latex_matrix = r"Z_{012} = \begin{bmatrix}"
    for fila in Z_seq_total:
        elementos = []
        for Z in fila:
            magnitud = abs(Z)
            angulo = math.degrees(math.atan2(Z.imag, Z.real))
            elementos.append(
                rf"{magnitud:.4f} \angle {angulo:.2f}^\circ"
            )
        latex_matrix += " & ".join(elementos)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    Z0_total = Z_seq_total[2, 2]
    print(f"\nImpedancia homopolar de la línea total (Z0): {Z0_total:.4f} Ω")
    Z1_total = Z_seq_total[1, 1]
    print(f"Impedancia directa e inversa de la línea total (Z1): {Z1_total:.4f} Ω")
    # Resistencias R0/R1
    R0 = Z0.real
    R1 = Z1.real
    print(f"\nR0/R1: {R0/R1:.3f}")
    # Reactancias X0/X1
    X0 = Z0.imag
    X1 = Z1.imag
    print(f"X0/X1: {X0/X1:.3f}")

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # MATRIZ DE CAPACIDADES
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    st.header("\n4. MATRIZ DE CAPACIDADES")
    st.markdown("---")
    st.write("""
    El cálculo de la matriz de capacidades de la línea se realiza mediante la matriz de coeficientes de potencial que liga el potencial eléctrico con la carga. La ecuación general en forma matricial es:
    """)
    st.latex(r"V = P \, q")
    st.markdown("""
    Donde:
    - $V$ es el vector de potenciales de la línea  
    - $P$ es la matriz de coeficientes de potencial  
    - $q$ es el vector de cargas  
    """)
    st.latex(r"q = C \, V = P^{-1} \, V")
    st.markdown("""
    Siendo $C$ la matriz de capacidades de la línea. Por tanto, la susceptancia capacitiva de la línea será:
    """)
    st.latex(r"B = \omega \, P^{-1}")
    st.markdown("Los elementos de la matriz $P$ de coeficientes de capacidades se pueden calcular:")
    st.write("Para elementos de la diagonal principal:")
    st.latex(r"P_{i,i} = \frac{1}{2 \pi \epsilon_0} \ln \left( \frac{D_{i,i}'}{r_\mathrm{eq}} \right)")
    st.write("Para elementos fuera de la diagonal:")
    st.latex(r"P_{i,j} = \frac{1}{2 \pi \epsilon_0} \ln \left( \frac{D_{i,i}'}{D_{i,j}} \right), \quad i \neq j")
    st.markdown("""
    Donde:  
    - $D_{i,j}$ es la distancia entre el conductor i y el j  
    - $D'_{i,j}$ es la distancia entre el conductor i y la imagen del conductor j respecto al suelo  
    - $r_\mathrm{eq}$ es el radio equivalente del conductor  
    - $\\epsilon_0$ es la permitividad eléctrica del vacío  

    Las matrices de distancias entre los conductores y entre los conductores y el espejo de los mismos respecto al suelo se han evaluado en un apartado anterior.
    """)

    st.write("La matriz de coeficientes de potencial de fase de la línea incluyendo cables de tierra será:")
    # Coeficientes de potencial
    # C_ii = (1/(2*pi*epsilon0))*log(D_prima_ii/r_i)
    # C_ij = (1/(2*pi*epsilon0))*log(D_prima_ij/D_ij) para i != j
    matriz_capacidades = [[0]*n_puntos for _ in range(n_puntos)]
    for i in range(n_puntos):
        for j in range(n_puntos):
            D_ij = matriz_distancias[i][j]/1000   # Convertir m a km
            D_prima_ij = matriz_D_prima[i][j]/1000   # Convertir m a km
            r_i = radio_conductores_mm[i] / 1000000  # Convertir mm a km
            if i == j:
                C_ii = ((1 / (2 * pi * epsilon0)) * log(D_prima_ij / r_i))/1000000000
                matriz_capacidades[i][j] = C_ii
            else:
                C_ij = (1 / (2 * pi * epsilon0)) * log(D_prima_ij / D_ij)/1000000000
                matriz_capacidades[i][j] = C_ij

    df=pd.DataFrame(matriz_capacidades).round(4)
    latex_matrix = r"P = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    latex_matrix += r"\;\text{km}/\mu\text{F}"
    st.latex(latex_matrix)
    # Matriz de capacitancias por fase (eliminando tierra)
    st.write("""
    Como la línea posee cables de tierra es necesario realizar un análisis matricial 
    para eliminarlos y obtener una matriz 3x3 que representa los coeficientes 
    de potencial de cada fase.
    """)
    st.latex(r"""
    P =
    \begin{pmatrix}
    P_f & P_{ft} \\
    P_{tf} & P_{tt}
    \end{pmatrix}
    """)
    st.write("""
    Aplicando la eliminación matricial (reducción de Kron), la matriz equivalente 
    de coeficientes de potencial de fase viene dada por:
    """)
    # Reducción de Kron
    st.latex(r"""
    P_{fas} =
    P_f
    -
    P_{ft}
    P_{tt}^{-1}
    P_{tf}
    """)
    st.write("""
    El resultado del cálculo es la matriz reducida de coeficientes de potencial de fase es:
    """)
    C = np.array(matriz_capacidades)
    Cf = C[0:3, 0:3]
    Ct = C[3:4, 3:4]
    Cft = C[0:3, 3:4]
    Ctf = C[3:4, 0:3]
    Ct_inv = np.linalg.inv(Ct)
    Cfas = Cf - Cft @ Ct_inv @ Ctf
    df=pd.DataFrame(Cfas).round(4)
    latex_matrix = r"P_{fas} = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    latex_matrix += r"\;\text{km}/\mu\text{F}"
    st.latex(latex_matrix)

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # MATRIZ DE SUSCEPTANCIAS
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    st.header("\n5. MATRIZ DE SUSCEPTANCIAS")
    st.markdown("---")
    st.write("""La matriz de susceptancias de la línea será el producto de la pulsación por la inversa de la matriz de coeficientes de potencial:""")
    st.latex(r"B = jwP^{-1}")
    # Matriz de susceptancias por km (B=j*2*pi*frecuencia*C^-1) solo parte imaginaria, no mostrando la parte real
    Cfas_inv = np.linalg.inv(Cfas)
    Bfas = 1j * 2 * pi * frecuencia * Cfas_inv
    df=pd.DataFrame(Bfas).round(4)
    latex_matrix = r"B = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    latex_matrix += r"\;\mu\text{S}/\text{km}"
    st.latex(latex_matrix)
    # Matriz de susceptancias total para la longitud dada
    Bfas_total = Bfas * longitud
    st.markdown(f"""Con la longitud de la línea de {longitud} km, se obtiene:""")
    df=pd.DataFrame(Bfas_total).round(4)
    latex_matrix = r"B = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    latex_matrix += r"\;\mu\text{S}"
    st.latex(latex_matrix)
    # Susceptancias de secuencia
    B_seq = A_inv @ Bfas @ A
    print("\nMatriz de susceptancias de secuencia:")
    df=pd.DataFrame(B_seq).round(4)
    latex_matrix = r"B_{012} = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    latex_matrix += r"\;\mu\text{S}"

    st.latex(latex_matrix)
    B0 = B_seq[2, 2]
    print(f"\nSusceptancia homopolar de la línea (B0): {B0:.4f} μS/km")
    B1 = B_seq[1, 1]
    print(f"Susceptancia directa e inversa de la línea (B1): {B1:.4f} μS/km")
    B0_total = B_seq[2, 2] * longitud
    print(f"\nSusceptancia homopolar de la línea total (B0): {B0_total:.4f} μS")
    B1_total = B_seq[1, 1] * longitud
    print(f"Susceptancia directa e inversa de la línea total (B1): {B1_total:.4f} μS")

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CÁLCULO IMPEDANCIA CARACTERÍSTICA Y CONSTANTE DE PROPAGACIÓN
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    st.header("\n6. IMPEDANCIA CARACTERÍSTICA Y CONSTANTE DE PROPAGACIÓN")
    st.markdown("---")
    # Impedancia característica (Zc=raiz((R+jX)/(jB)))
    Zc = cmath.sqrt(Z1 / (1 * B1*1e-6))
    magnitud_Zc = abs(Zc)
    angulo_Zc = math.degrees(math.atan2(Zc.imag, Zc.real))
    print(f"\nImpedancia característica de secuencia directa (Zc): {Zc:.4f} = {magnitud_Zc:.4f} ∠ {angulo_Zc:.2f}° Ω")
    # Constante de propagación (gamma=raiz((R+jX)*(jB)))
    gamma = cmath.sqrt(Z1 * (1 * B1*1e-6))
    magnitud_gamma = abs(gamma)
    angulo_gamma = math.degrees(math.atan2(gamma.imag, gamma.real))
    print(f"\nConstante de propagación (γ): {gamma:.6f} = {magnitud_gamma:.6f} ∠ {angulo_gamma:.2f}° 1/km")
    # Constante de propagación total para la longitud dada
    gamma_total = gamma * longitud
    magnitud_gamma_total = abs(gamma_total)
    angulo_gamma_total = math.degrees(math.atan2(gamma_total.imag, gamma_total.real))
    print(f"\nConstante de propagación para la longitud dada: {gamma_total:.6f} = {magnitud_gamma_total:.6f} ∠ {angulo_gamma_total:.2f}°")

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # POTENCIA CARACTERÍSTICA
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    st.header("\n7. POTENCIA CARACTERÍSTICA")
    st.markdown("---")
    print("\nLa potencia característica de la línea es función de la tensión y de la impedancia característica a través de la siguiente expresión:")
    st.latex(r"P_c=\frac{U^2}{Z_c}")
    # Pc = tension_nominal**2 / Zc
    Pc = (tension_nominal) ** 2 / Zc
    magnitud_Pc = abs(Pc)
    angulo_Pc = math.degrees(math.atan2(Pc.imag, Pc.real))
    print(f"\nPotencia característica Pc: {Pc:.4f} = {magnitud_Pc:.4f} ∠ {angulo_Pc:.2f}° MVA")
    st.write("\nPara los valores de la línea:")
    st.latex(
        rf"P_c=\frac{{U^2}}{{Z_c}}"
        rf"=\frac{{{tension_nominal}^2}}{{{Zc:.4f}}}"
        rf"={Pc:.4f}"
        rf"={magnitud_Pc:.4f}\angle {angulo_Pc:.2f}^\circ \,\text{{MVA}}"
    )
    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CAIDA DE TENSIÓN
    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    st.header("\n8. CAIDA DE TENSIÓN")
    st.markdown("---")
    print("\nLa caída de tensión por resistencia y reactancia de una línea (despreciendo la influencia de la capacidad) viene dada por las fórmulas:")
    st.latex(r"\Delta U = \frac{P \cdot L \cdot (R \cdot \cos(\varphi) + X \cdot \sin(\varphi))}{U^2 \cdot 10 \cdot \cos(\varphi)}")
    # AU=potencia_transportada(MW)*longitud*(R*cos(phi)+X*sin(phi))/((tension_nominal**2)*10*cos(phi))
    potencia_transportada_MW = potencia_transportada * cos_phi  # Convertir MVA a MW usando el factor de potencia
    ΔU = potencia_transportada_MW * 1000 * longitud * (Z1.real * cos_phi + Z1.imag * math.sin(math.acos(cos_phi))) / ((tension_nominal ** 2) * 10 * cos_phi)
    st.latex(
        rf"\Delta U = \frac{{{potencia_transportada_MW:.2f} \cdot {longitud} \cdot ({Z1.real:.4f} \cdot {cos_phi} + {Z1.imag:.4f} \cdot {math.sin(math.acos(cos_phi)):.3f})}}{{{tension_nominal}^2 \cdot 10 \cdot {cos_phi}}}"
        rf"={ΔU:.4f} %"
    )
    print(f"\nCaída de tensión ΔU: {ΔU:.4f} %")
    # Poner si la caida de tension es inferior al 5%: La caída de tensión es inferior al 5%; y si no lo es: La caída de tensión es superior al 5%
    if ΔU < 5:
        print("La caída de tensión es inferior al 5%")
    else:
        print("La caída de tensión es superior al 5%")

    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # DENSIDAD MÁXIMA DE CORRIENTE E INTENSIDAD MÁXIMA POR CABLE
    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    st.header("\n9. DENSIDAD MÁXIMA DE CORRIENTE E INTENSIDAD MÁXIMA POR CABLE")
    st.markdown("---")
    # Densidades de corriente
    secciones = [15, 25, 35, 50, 70, 95, 125, 160, 200, 250, 300, 400, 500, 600]
    densidades = [6, 5, 4.55, 4.0, 3.55, 3.2, 2.9, 2.7, 2.5, 2.3, 2.15, 1.95, 1.8, 1.65]  # A/mm2
    def densidad_corriente(seccion):
        if seccion < secciones[0] or seccion > secciones[-1]:
            return None  # Fuera de rango
        for i in range(len(secciones) - 1):
            if secciones[i] <= seccion <= secciones[i + 1]:
                # Regla de tres
                densidad = densidades[i] + (densidades[i + 1] - densidades[i]) * (seccion - secciones[i]) / (secciones[i + 1] - secciones[i])
                return densidad
        return None
    densidad_maxima = densidad_corriente(seccion)
    print(f"\nDensidad máxima de corriente para la sección {seccion} mm²: {densidad_maxima:.4f} A/mm²")
    # Coeficientes reductores
    composiciones = ["30+7", "6+1", "26+7", "54+7", "45+7"]
    coef_reductores = [0.916, 0.937, 0.937, 0.95, 0.97]
    def coeficiente_reductor(composicion):
        if composicion in composiciones:
            index = composiciones.index(composicion)
            return coef_reductores[index]
        return None
    coef_reductor = coeficiente_reductor(composicion)
    print(f"Coeficiente reductor para la composición {composicion}: {coef_reductor:.4f}")
    densidad_max_con_reduccion = densidad_maxima * coef_reductor
    print(f"Densidad máxima de corriente con coeficiente reductor: {densidad_max_con_reduccion:.4f} A/mm²")
    # Intensidad máxima
    intensidad_maxima_conductor = densidad_max_con_reduccion * seccion  # A
    print(f"Corriente máxima admisible del conductor: {intensidad_maxima_conductor:.2f} A")

    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # POTENCIA MÁXIMA ADMISIBLE POR INTENSIDAD
    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    st.header("\n10. POTENCIA MÁXIMA ADMISIBLE POR INTENSIDAD")
    st.markdown("---")
    potencia_maxima_admisible = (intensidad_maxima_conductor * tension_nominal * raiz(3)) * cos_phi 
    print(f"\nPotencia máxima admisible por intensidad: {potencia_maxima_admisible/1000:.2f} MW")
    # Comparar potencia máxima admisible con potencia transportada
    if potencia_maxima_admisible/1000 >= potencia_transportada_MW:
        print("La potencia es mayor que la potencia a transportar.")
    else:
        print("La potencia es menor que la potencia a transportar.")

    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # POTENCIA DE TRANSPORTE EN FUNCIÓN DE CONDICIONES METEOROLÓGICAS
    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    st.header("\n11. CÁLCULO POTENCIA DE TRANSPORTE EN FUNCIÓN DE CONDICIONES METEOROLÓGICAS")
    st.markdown("---")
    # Calor aportado por radiación solar para invierno y verano Qs = coeficiente_absorcion * radiacion * diametro/1000
    st.subheader("\n   11.1. CALOR APORTADO POR RADIACIÓN SOLAR")
    Qs_invierno = coeficiente_absorcion * radiacion_invierno * diametro/1000 # W/m
    Qs_verano = coeficiente_absorcion * radiacion_verano * diametro/1000 # W/m
    print(f"\nCalor aportado por radiación solar en invierno: {Qs_invierno:.3f} W/m")
    print(f"Calor aportado por radiación solar en verano: {Qs_verano:.3f} W/m")
    # Calor cedido por radiación Qr = pi * diametro/1000 * emisividad_conductor * cte_boltzmann * ((Tc+273.15)^4 - (Ta+273.15)^4)

    st.subheader("\n   11.2. CALOR CEDIDO POR RADIACIÓN SOLAR")
    cte_boltzmann = 5.6704e-8 # W/m2K4
    Qr_invierno = pi * (diametro/1000) * emisividad_conductor * cte_boltzmann * ((Tc + 273.15)**4 - (temperatura_invierno + 273.15)**4) # W/m
    Qr_verano = pi * (diametro/1000) * emisividad_conductor * cte_boltzmann * ((Tc + 273.15)**4 - (temperatura_verano + 273.15)**4) # W/m
    print(f"\nCalor cedido por radiación en invierno: {Qr_invierno:.3f} W/m")
    print(f"Calor cedido por radiación en verano: {Qr_verano:.3f} W/m")

    # Calor cedido por convección Qc = pi * conductividad_termica * (Tc - Ta) * Nu
    st.subheader("\n   11.3. CALOR CEDIDO POR CONVECCIÓN")
    # 1.  Convección natural: Nu = A * (Gr*Pr)^m
    print("\n       11.3.1. CONVECCIÓN NATURAL")
    # Valores de A y m dependiendo de Gr*Pr
    def coeficientes_conveccion_natural(Gr_Pr):
        if 0.1 <= Gr_Pr < 100:
            return 1.02, 0.148
        elif 100 <= Gr_Pr < 10000:
            return 0.85, 0.188
        elif 10000 <= Gr_Pr < 1e7:
            return 0.48, 0.25
        elif 1e7 <= Gr_Pr < 1e12:
            return 0.125, 0.333
        else:
            return None, None
    # Gr = (diametro/1000)^3 * (Tc - Ta) * g / ((Tav + 27.15) * (viscosidad_cinematica)^2)
    g = 9.81  # m/s2
    # viscosidad_cinematica = 1.32e-5 + 9.5e-8 * Tav
    def viscosidad_cinematica(Tav):
        return 1.32e-5 + 9.5e-8 * Tav
    viscosidad_cinematica_invierno = viscosidad_cinematica((Tc + temperatura_invierno) / 2)
    viscosidad_cinematica_verano = viscosidad_cinematica((Tc + temperatura_verano) / 2)
    print(f"\nViscosidad cinemática en invierno: {viscosidad_cinematica_invierno:.4e} m²/s")
    print(f"Viscosidad cinemática en verano: {viscosidad_cinematica_verano:.4e} m²/s")
    # densidad_relativa_aire = exp(-1.16e-4* h)
    def Gr(Ta, viscosidad_cinematica):
        Tav = (Tc + Ta) / 2
        return ((diametro/1000)**3 * (Tc - Ta) * g) / ((Tav + 273.15) * (viscosidad_cinematica**2))
    Gr_invierno = Gr(temperatura_invierno, viscosidad_cinematica_invierno)
    Gr_verano = Gr(temperatura_verano, viscosidad_cinematica_verano)
    print(f"\nNúmero de Grashof en invierno: {Gr_invierno:.2f}")
    print(f"Número de Grashof en verano: {Gr_verano:.2f}")
    # Pr = calor_esp_aire * viscosidad_dinamica / conductividad_termica_aire
    calor_esp_aire = 1005  # J/kgK
    # viscosidad_dinamica = (17.239 + 4.625e-2 * Tav - 2.03e-5 * Tav **2)e-6 kg/m s
    def viscosidad_dinamica(Tav):
        return (17.239 + 4.625e-2 * Tav - 2.03e-5 * Tav **2) * 1e-6  # kg/m s
    viscosidad_dinamica_invierno = viscosidad_dinamica((Tc + temperatura_invierno) / 2)
    viscosidad_dinamica_verano = viscosidad_dinamica((Tc + temperatura_verano) / 2)
    print(f"\nViscosidad dinámica en invierno: {viscosidad_dinamica_invierno:.6e} kg/m s")
    print(f"Viscosidad dinámica en verano: {viscosidad_dinamica_verano:.6e} kg/m s")
    # conductividad_termica_aire = 2.368e-2 + 7.23e-5 * Tav - 2.763e-8 * Tav**2 W/K·m
    def conductividad_termica_aire(Tav):
        return 2.368e-2 + 7.23e-5 * Tav - 2.763e-8 * Tav**2  # W/K·m
    conductividad_termica_aire_invierno = conductividad_termica_aire((Tc + temperatura_invierno) / 2)
    conductividad_termica_aire_verano = conductividad_termica_aire((Tc + temperatura_verano) / 2)
    print(f"\nConductividad térmica del aire en invierno: {conductividad_termica_aire_invierno:.6e} W/K·m")
    print(f"Conductividad térmica del aire en verano: {conductividad_termica_aire_verano:.6e} W/K·m")
    Pr_invierno = calor_esp_aire * viscosidad_dinamica_invierno / conductividad_termica_aire_invierno
    Pr_verano = calor_esp_aire * viscosidad_dinamica_verano / conductividad_termica_aire_verano
    print(f"\nNúmero de Prandtl en invierno: {Pr_invierno:.4f}")
    print(f"Número de Prandtl en verano: {Pr_verano:.4f}")
    Gr_Pr_invierno = Gr_invierno * Pr_invierno
    Gr_Pr_verano = Gr_verano * Pr_verano
    print(f"\nNúmero Gr*Pr en invierno: {Gr_Pr_invierno:.2f}")
    print(f"\nNúmero Gr*Pr en verano: {Gr_Pr_verano:.2f}")
    Nu_invierno = 0
    Nu_verano = 0
    A_invierno, m_invierno = coeficientes_conveccion_natural(Gr_Pr_invierno)
    if A_invierno is not None:
        Nu_invierno = A_invierno * (Gr_Pr_invierno ** m_invierno)
    A_verano, m_verano = coeficientes_conveccion_natural(Gr_Pr_verano)
    if A_verano is not None:
        Nu_verano = A_verano * (Gr_Pr_verano ** m_verano)
    print(f"\nNúmero de Nusselt en invierno: {Nu_invierno:.4f}")
    print(f"Número de Nusselt en verano: {Nu_verano:.4f}")

    # 2. Convección forzada: Nu = B1 * Re^n
    print("\n       11.3.2. CONVECCIÓN FORZADA")
    # Re = (diametro/1000) * velocidad_viento / viscosidad_cinematica
    # Rugosidad: Rf = diametro_alambre_ext / (2*(diametro-diametro_alambe_ext))
    Rf = diametro_alambre_ext / (2 * (diametro - diametro_alambre_ext))
    print(f"\nRugosidad Rf: {Rf:.3f}")
    def Re(viscosidad_cinematica):
        return (diametro/1000) * velocidad_viento / viscosidad_cinematica
    Re_invierno = Re(viscosidad_cinematica_invierno)
    Re_verano = Re(viscosidad_cinematica_verano)
    print(f"\nNúmero de Reynolds en invierno: {Re_invierno:.2f}")
    print(f"Número de Reynolds en verano: {Re_verano:.2f}")
    #B1 y n según Re
    def coeficientes_conveccion_forzada(Re, Rf):
        if 100 < Re < 2.65e3:
            return 0.641, 0.471
        elif Rf <= 0.05 and 2.65e3 <= Re < 5e4:
            return 0.178, 0.633
        elif Rf > 0.05 and 2.65e3 <= Re < 5e4:
            return 0.048, 0.8
        else:
            return None, None
    Nu_invierno_forzada = 0
    Nu_verano_forzada = 0
    B1_invierno, n_invierno = coeficientes_conveccion_forzada(Re_invierno, Rf)
    if B1_invierno is not None:
        Nu_invierno_forzada = B1_invierno * (Re_invierno ** n_invierno)
    B1_verano, n_verano = coeficientes_conveccion_forzada(Re_verano, Rf)
    if B1_verano is not None:
        Nu_verano_forzada = B1_verano * (Re_verano ** n_verano)
    print(f"\nNúmero de Nusselt por convección forzada en invierno: {Nu_invierno_forzada:.4f}")
    print(f"Número de Nusselt por convección forzada en verano: {Nu_verano_forzada:.4f}")
    Nu_45_invierno = Nu_invierno_forzada * (0.42 + 0.58*math.sin(45*pi/180)**0.9)
    Nu_45_verano = Nu_verano_forzada * (0.42 + 0.58*math.sin(45*pi/180)**0.9)
    print(f"\nNúmero de Nusselt corregido para ángulo de 45° en invierno: {Nu_45_invierno:.4f}")
    print(f"Número de Nusselt corregido para ángulo de 45° en verano: {Nu_45_verano:.4f}")
    # Qc = pi* conductividad_termica * (Tc-Ta) * Nu con Nu=max(Nu_natural, Nu_45_forzada)
    Qc_invierno = pi * conductividad_termica_aire_invierno * (Tc - temperatura_invierno) * max (Nu_invierno, Nu_45_invierno)
    Qc_verano = pi * conductividad_termica_aire_verano * (Tc - temperatura_verano) * max (Nu_verano, Nu_45_verano)
    print(f"\nCalor cedido por convección en invierno: {Qc_invierno:.3f} W/m")
    print(f"Calor cedido por convección en verano: {Qc_verano:.3f} W/m")

    # Resultados corriente máxima: I = raiz((Qr+Qc-Qs)/resistencia_ca/1000)
    st.subheader("\n   11.4. RESULTADOS CORRIENTE MÁXIMA")
    print(f"\nResistencia del conductor a la temperatura de cálculo: {resistencia_ca:.6f} Ohmios/km")
    I_max_invierno = raiz((Qr_invierno + Qc_invierno - Qs_invierno)/ (resistencia_ca * 1e-3))
    I_max_verano = raiz((Qr_verano + Qc_verano - Qs_verano) / (resistencia_ca * 1e-3))
    print(f"\nCorriente máxima admisible en invierno según condiciones meteorológicas: {I_max_invierno:.2f} A")
    print(f"Corriente máxima admisible en verano según condiciones meteorológicas: {I_max_verano:.2f} A")

    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # POTENCIA MÁXIMA DE TRANSPORTE SEGÚN CONDICIONES METEOROLÓGICAS
    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    st.subheader("\n   11.5. POTENCIA MÁXIMA DE TRANSPORTE")
    potencia_maxima_invierno = (I_max_invierno * tension_nominal * raiz(3)) * cos_phi
    potencia_maxima_verano = (I_max_verano * tension_nominal * raiz(3)) * cos_phi
    print(f"\nPotencia máxima de transporte en invierno según condiciones meteorológicas: {potencia_maxima_invierno/1000:.2f} MW")
    print(f"Potencia máxima de transporte en verano según condiciones meteorológicas: {potencia_maxima_verano/1000:.2f} MW")

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # PÉRDIDAS DE POTENCIA
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    st.header("\n12. PÉRDIDAS DE POTENCIA")
    st.markdown("---")
    print(f"Pot transportada mw= {potencia_transportada_MW}")
    perdidas_potencia = ((potencia_transportada_MW * resistencia_ca * longitud) / (tension_nominal**2 * cos_phi**2))*100
    print(f"\nPérdidas de potencia en la línea: {perdidas_potencia:.5f} %")
    # En valor absoluto
    perdidas_potencia_valor = (potencia_transportada_MW * resistencia_ca * longitud) / (tension_nominal**2 * cos_phi**2) * potencia_transportada_MW
    print(f"Pérdidas de potencia en la línea en valor absoluto: {perdidas_potencia_valor:.5f} MW")

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CORTOCIRCUITO MÁXIMO
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    st.header("\n13. CORTOCIRCUITO MÁXIMO")
    st.markdown("---")
    # Constantes según material
    materiales = ["Cobre", "Aluminio-Acero", "Acero"]
    c_conductores = [390, 910, 480]  # J/kgK
    densidades_conductores = [8900, 2700, 7850]  # kg/m3
    ctes_conductores = [56e6, 34.8e6, 7.25e6]  # 1/ohmio m
    alphas_conductores = [0.0039, 0.004, 0.0045]  # 1/ºC
    def propiedades_material(material):
        if material in materiales:
            index = materiales.index(material)
            return (c_conductores[index], densidades_conductores[index], ctes_conductores[index], alphas_conductores[index])
        return (None, None, None, None)
    c_conductor, densidad_conductor, cte_conductor, alpha_conductor = propiedades_material(material)
    print(f"c_conductor, densidad_conductor, cte_conductor, alpha_conductor", c_conductor, densidad_conductor, cte_conductor, alpha_conductor )
    # Temperatura máxima recomendada para el material
    if material == "Acero":
        temperatura_max_recomendada = 300
    else:
        temperatura_max_recomendada = 200
    # Factor K
    multiplicacion = log ((1+alpha_conductor*(temperatura_max_recomendada-20))/(1+alpha_conductor*(Tc-20)))
    print(f"multiplicación:", multiplicacion)
    K = (raiz (cte_conductor*c_conductor*densidad_conductor*multiplicacion/alpha_conductor))/1000000
    print(f"Factor K: {K} A*(s)^(1/2)/mm²")
    Icc_max = K * seccion / raiz(tiempo_accionamiento_proteccion)
    print(f"Icc max: {Icc_max} A")

    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # EFECTO CORONA
    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    st.header("\n14. EFECTO CORONA ")
    st.markdown("---")
    beta = 1 # Coeficiente reductor para conductores múltiples
    presion_barometrica = 1/(log(log(76)-altitud_media/18330))*100 
    print(f"Presión barométrica {presion_barometrica} cmHg")
    # Factor corrección densidad aire: delta
    def delta (temperatura):
        return 3.92*presion_barometrica/(273+temperatura)
    delta_invierno = delta (temperatura_invierno)
    delta_verano = delta (temperatura_verano)
    print(f"Factor corrección densidad aire invierno {delta_invierno}")
    print(f"Factor corrección densidad aire verano {delta_verano}")
    # Calcula de la distancia media geométrica DMG
    DMG = 100 * ((matriz_distancias[0][1]*matriz_distancias[1][2]*matriz_distancias[2][0]*matriz_distancias[0][3]*matriz_distancias[1][3]*matriz_distancias[2][3])**(1/6))
    print(f"DMG: {DMG:.3f} cm") # con tierra
    DMG = 100 * ((matriz_distancias[0][1]*matriz_distancias[1][2]*matriz_distancias[2][0])**(1/3)) # cm
    print(f"DMG: {DMG:.1f} cm") # entre fases
    DMG = 807
    print(f"DMG: {DMG}")
    # Tensión crítica dieléctrica
    def Uc (mt,delta):
        return raiz(3)*mc*mt*(30/raiz(2))*delta*diametro/20*log(DMG/(diametro/20))
    Uc_invierno = Uc(mt_invierno, delta_invierno)
    Uc_verano = Uc(mt_verano,delta_verano)
    print(f"Tensión critica invierno {Uc_invierno:.2f} kV")
    print(f"Tensión crítica verano {Uc_verano:.2f} kV")
    if Uc_invierno < tension_nominal:
        print(f"Hay pérdidas por efecto corona en invierno")
        Perdidas_efecto_corona_invierno = 241/delta_invierno * (frecuencia + 25) * raiz (diametro/20/DMG) * ((tension_nominal - Uc_invierno)/raiz(3))**2 * 1e-5
        print(f"Pérdidas por efecto corona en invierno: {Perdidas_efecto_corona_invierno:.2f} kW·km/fase")
    else: print(f"No hay pérdidas por efecto corona en invierno")
    if Uc_verano < tension_nominal:
        print(f"Hay pérdidas por efecto corona en verano")
        Perdidas_efecto_corona_verano = 241/delta_verano * (frecuencia + 25) * raiz (diametro/20/DMG) * (tension_nominal - Uc_verano)/raiz(3) * 1e-5
        print(f"Pérdidas por efecto corona en verano: {Perdidas_efecto_corona_verano:.2f} kW·km/fase")
    else: print(f"No hay pérdidas por efecto corona en verano")

    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CAMPO ELÉCTRICO
    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    st.header("\n15. CAMPO ELÉCTRICO ")
    st.markdown("---")
    # Matriz coeficientes de potencial
    print("\nMatriz de capacidades (km/uF):")
    df=pd.DataFrame(matriz_capacidades).round(4)
    latex_matrix = r"P = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # Matriz coeficientes de potencial inversa
    matriz_capacidades_inversa = np.linalg.inv(matriz_capacidades) * 1000
    print("\nMatriz de capacidades inversa(nF/km):")
    df=pd.DataFrame(matriz_capacidades_inversa).round(4)
    latex_matrix = r"P^{-1} = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # Vector de potenciales (tensiones nominales/raiz(3) (angulos 0 -120 y 120, tension tierra)) en numeros complejos y kV
    V_phase = tension_nominal / raiz(3)
    V_vector = np.array([V_phase * cmath.rect(1, math.radians(0)),
                        V_phase * cmath.rect(1, math.radians(-120)),
                        V_phase * cmath.rect(1, math.radians(120)),
                        0])
    print("\nVector de potenciales (kV):")
    df=pd.DataFrame(V_vector).round(4)
    latex_matrix = r"U = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # Vector de cargas (matriz_capacidades_inversa @ V_vector)
    Q_vector = matriz_capacidades_inversa @ V_vector /1000
    print("\nVector de cargas (kV·nF/km):")
    df=pd.DataFrame(Q_vector).round(4)
    latex_matrix = r"q = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)

    # Valor complejo de la corriente que circula por cada conductor
    I_phase = potencia_transportada * 1e6 / (raiz(3) * tension_nominal*1e3)
    I_vector = np.array([I_phase * cmath.rect(1, math.radians(0)),
                        I_phase * cmath.rect(1, math.radians(-120)),
                        I_phase * cmath.rect(1, math.radians(120))])
    print("\nVector de corrientes (A):")
    df=pd.DataFrame(I_vector).round(4)
    latex_matrix = r"I = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
    latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)


    resistencia_ca_final = resistencia * (1+7.5*frecuencia**2*(diametro/20)**4*1e-7)
    print(f"\nResistencia:{resistencia_ca_final:.6f}")
    r=resistencia_ca_final*(1+alpha_conductor*(Tc-20))
    print(f"\nResist: {r:.6f}")


    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # EJECUCIÓN PROGRAMA
    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    sys.stdout = sys.__stdout__
    # Ahora redirigir al buffer para capturar el resto de los cálculos
    sys.stdout = buffer
    sys.stdout = sys.__stdout__
    texto_resultados = buffer.getvalue()
    html = """
    <html>
    <head>
    <meta charset="utf-8">
    <style>
    body { font-family: Arial; }
    .mayus { font-weight: bold; }
    </style>
    </head>
    <body>
    <pre>
    """
    def es_mayusculas(linea):
        letras = [c for c in linea if c.isalpha()]
        return len(letras)>0 and all(c.isupper() for c in letras)
    for linea in texto_resultados.split("\n"):
        # Detectar líneas completamente en mayúsculas
        if es_mayusculas(linea):
            html += f'<span class="mayus">{linea}</span>\n'
        else:
            html += linea + "\n"
    html += """
    </pre>
    </body>
    </html>
    """
    with open("resultados.html", "w", encoding="utf-8") as f:
        f.write(html)


    from docx import Document
    doc = Document()
    doc.add_heading("Resultados", 1)
    for linea in texto_resultados.split("\n"):
        doc.add_paragraph(linea)
    doc.save("resultados.docx")












#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DOBLE CIRCUITO SIMPLEX
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

elif ndecircuitos==2 and ndeconductoresporfase==1:
    # Cálculo de la resistencia a 85ºC y en ca
    print("\n1. RESISTENCIA ELÉCTRICA DE LA LÍNEA")
    def resistencia_a_temp(resistencia_20C, temp):
        return resistencia_20C * (1 + coef_temp * (temp - 20))
    resistencia_85C = resistencia_a_temp(resistencia, Tc)
    print(f"Resistencia a {Tc}ºC: {resistencia_85C:.6f} Ω")
    reactancia_pelicular_85C = 10**-7 * 8 * pi * frecuencia * (1 / (resistencia_a_temp(resistencia, Tc) / 1000))
    print(f"Reactancia por efecto pelicular a {Tc}ºC: {reactancia_pelicular_85C:.6f} Ω")
    ys = reactancia_pelicular_85C ** 2 / (192 + 0.8*reactancia_pelicular_85C**2)
    resistencia_ca = resistencia_85C * (1 + ys)
    print(f"Resistencia en CA a {Tc}ºC: {resistencia_ca:.6f} Ω")
    r_ca_longitud = resistencia_ca * longitud
    print(f"Resistencia en CA para {longitud} km: {r_ca_longitud:.6f} Ω")
    
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CÁLCULO MATRIZ IMPEDANCIAS
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print("\n2. MATRIZ DE IMPEDANCIAS")
    penetracion_terreno = raiz(resistividad/(pi*frecuencia*mu0))
    print(f"Penetración terreno: {penetracion_terreno:.6f} m")
    # Matriz de distancias
    puntos = [uno, dos, tres, cuatro, cinco, seis, tierra]
    n_puntos = len(puntos)
    matriz_distancias = [[0]*n_puntos for _ in range(n_puntos)]
    def _xy(pt):
        try:
            if isinstance(pt, (list, tuple)) and len(pt) >= 2:
                return float(pt[0]), float(pt[1])
            if isinstance(pt, str):
                s = pt.strip().replace(';', ',')
                try:
                    v = ast.literal_eval(s)
                    if isinstance(v, (list, tuple)) and len(v) >= 2:
                        return float(v[0]), float(v[1])
                except Exception:
                    parts = s.split(',')
                    if len(parts) >= 2:
                        return float(parts[0]), float(parts[1])
            if isinstance(pt, (int, float)):
                return float(pt), 0.0
        except Exception:
            pass
        raise ValueError(f"Formato de punto inválido: {pt}")
    for i in range(n_puntos):
        for j in range(n_puntos):
            x1, y1 = _xy(puntos[i])
            x2, y2 = _xy(puntos[j])
            distancia = raiz((x2-x1)**2 + (y2-y1)**2)
            matriz_distancias[i][j] = distancia
    print("\nMatriz de distancias (km):")
    df = pd.DataFrame(matriz_distancias).round(4)
    latex_matrix = r"D = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
        latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)

    # Matriz D_prima: distancia entre puntos y sus espejos respecto al suelo
    matriz_D_prima = [[0]*n_puntos for _ in range(n_puntos)]
    for i in range(n_puntos):
        for j in range(n_puntos):
            x1, y1 = _xy(puntos[i])
            # Espejo del punto j respecto al suelo
            x2, y2 = _xy(puntos[j])
            x2_espejo, y2_espejo = x2, -y2
            distancia = raiz((x2_espejo-x1)**2 + (y2_espejo-y1)**2)
            matriz_D_prima[i][j] = distancia
    print("\nMatriz D' (distancias entre puntos y espejos) (km):")
    df = pd.DataFrame(matriz_D_prima).round(4)
    latex_matrix = r"D' = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
        latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # kij = raiz(2)*D_prima_ij / penetracion_terreno
    matriz_kij = [[0]*n_puntos for _ in range(n_puntos)]
    for i in range(n_puntos):
        for j in range(n_puntos):
            matriz_kij[i][j] = raiz(2) * matriz_D_prima[i][j] / penetracion_terreno  
    print("\nMatriz k_ij:")
    df = pd.DataFrame(matriz_kij).round(4)
    latex_matrix = r"k_{ij} = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
        latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # P_ij = pi/8 - k_ij*cos(tetha_ij)/(3*raiz(2)) siendo cos(tetha_ij) = (h_i + h_j )/ D_prima_ij
    matriz_Pij = [[0]*n_puntos for _ in range(n_puntos)]
    for i in range(n_puntos):
        for j in range(n_puntos):
            h_i = _xy(puntos[i])[1]
            h_j = _xy(puntos[j])[1]
            D_prima_ij = matriz_D_prima[i][j]
            if D_prima_ij != 0:
                cos_tetha_ij = (h_i + h_j) / D_prima_ij
            else:
                cos_tetha_ij = 0
            k_ij = matriz_kij[i][j]
            P_ij = (pi / 8) - (k_ij * cos_tetha_ij) / (3 * raiz(2))
            matriz_Pij[i][j] = P_ij
    print("\nMatriz P_ij:")
    df = pd.DataFrame(matriz_Pij).round(4)
    latex_matrix = r"P_{ij} = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
        latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # Q_ij = 0.5*log neperiano(1.85138/k_ij) + k_ij*cos(tetha_ij)/(3*raiz(2))
    matriz_Qij = [[0]*n_puntos for _ in range(n_puntos)]
    for i in range(n_puntos):
        for j in range(n_puntos):
            k_ij = matriz_kij[i][j]
            h_i = _xy(puntos[i])[1]
            h_j = _xy(puntos[j])[1]
            D_prima_ij = matriz_D_prima[i][j]
            if D_prima_ij != 0:
                cos_tetha_ij = (h_i + h_j) / D_prima_ij
            else:
                cos_tetha_ij = 0
            if k_ij != 0:
                Q_ij = 0.5 * log(1.85138 / k_ij) + (k_ij * cos_tetha_ij) / (3 * raiz(2))              
            else:
                Q_ij = 0
            matriz_Qij[i][j] = Q_ij
    print("\nMatriz Q_ij:")
    df = pf.DtaFrame(matriz_Qij).round(4)
    latex_matrix = r"Q_{ij} = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
        latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # Resistencia de cada uno de los conductores (ultimo a tierra con dato de resistencia a tierra directamente)
    resistencias_conductores = [resistencia_ca, resistencia_ca, resistencia_ca, resistencia_ca, resistencia_ca, resistencia_ca, resistencia_tierra]
    print("\nResistencias de los conductores (Ω/km):")
    for i, R in enumerate(resistencias_conductores):
        print(f" {R:.6f} Ω/km")
    # Radio de cada conductor en mm
    radio_conductores_mm = [
        diametro / 2,
        diametro / 2,
        diametro / 2,
        diametro / 2,
        diametro / 2,
        diametro / 2,
        diametro_tierra / 2
    ]
    print("\nRadio de los conductores (mm):")
    df=pd.DataFrame(radio_conductores_mm).round(4)
    latex_matrix = r"r = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
        latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # Matriz de impedancias con:
    # Z_ii = R_i + j*mu0*frecuencia*((1/(4*n_i)+log(D_prima_ii/r_i)))+mu0*2*frecuencia*(P_ii + j*Q_ii) 
    # Z_ij = j*mu0*frecuencia*(log(D_prima_ij/D_ij)) + mu0*2*frecuencia*(P_ij + j*Q_ij) para i != j
    matriz_impedancias = [[0]*n_puntos for _ in range(n_puntos)]
    for i in range(n_puntos):
        for j in range(n_puntos):
            R_i = resistencias_conductores[i]
            D_ij = matriz_distancias[i][j]
            D_prima_ij = matriz_D_prima[i][j]
            P_ij = matriz_Pij[i][j]
            Q_ij = matriz_Qij[i][j]
            r_i = radio_conductores_mm[i] / 1000  # Convertir mm a m
            if i == j:
                n_i = ndeconductoresporfase
                Z_ii = R_i + (1j * mu0 * frecuencia * ((1 / (4 * n_i)) + log(D_prima_ij / r_i)) + mu0 * 2 * frecuencia * (P_ij + 1j * Q_ij))*1000
                matriz_impedancias[i][j] = Z_ii
            else:
                Z_ij = (1j * mu0 * frecuencia * (log(D_prima_ij / D_ij)) + mu0 * 2 * frecuencia * (P_ij + 1j * Q_ij))*1000
                matriz_impedancias[i][j] = Z_ij
    print("\nMatriz de impedancias (Ω/km):")
    for i, fila in enumerate(matriz_impedancias):
        print(f"", end="")
        for Z in fila:
            print(f"{Z:18.4f}", end="  ")
        print()
    df = pd.DataFrame(matriz_impedancias).round(4)
    latex_matrix = r"Z = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
        latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    """ Como la linea posee cable de tierra, es necesario realizar un análisis matricial para eliminarlos y obtener una matriz 3*3 que representa las impedancias por fase
    Z=[Zf, Zft; Ztf, Zt]
    Zfas = Zf - Zft * inv(Zt) * Ztf
    """
    Z = np.array(matriz_impedancias)
    Zf = Z[0:6, 0:6]
    Zt = Z[6:7, 6:7]
    Zft = Z[0:6, 6:7]
    Ztf = Z[6:7, 0:6]
    Zt_inv = np.linalg.inv(Zt)
    Zfas = Zf - Zft @ Zt_inv @ Ztf
    print("\nMatriz de impedancias por fase (Ω/km):")
    df = pd.DataFrame(Zfas).round(4)
    latex_matrix = r"Z = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.4f}" for v in row)
        latex_matrix += r" \\ "
        latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # Finalmente, multiplicar por la longitud de la línea para obtener las impedancias totales
    Zfas_total = Zfas * longitud
    print("\nMatriz de impedancias por fase total para la longitud dada (Ω):")
    df = pd.DataFrame(Zfas_total).round(4)
    latex_matrix = r"Z = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.3f}" for v in row)
        latex_matrix += r" \\ "
        latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # Impedancias de secuencia
    print("\n3. IMPEDANCIAS DE SECUENCIA")
    A = np.array([[1, 1, 1, 0, 0, 0],
                [complex(-0.5, -raiz(3)/2), complex(-0.5, raiz(3)/2), 1, 0, 0, 0],
                [complex(-0.5, raiz(3)/2), complex(-0.5, -raiz(3)/2), 1, 0, 0, 0],
                [0, 0, 0, 1, 1, 1],
                [0, 0, 0, complex(-0.5, -raiz(3)/2), complex(-0.5, raiz(3)/2), 1],
                [0, 0, 0, complex(-0.5, raiz(3)/2), complex(-0.5, -raiz(3)/2), 1]])
    A_inv = np.linalg.inv(A)
    Z_seq = A_inv @ Zfas @ A
    print("\nMatriz de impedancias de secuencia (Ω):")
    df = pd.DataFrame(Z_seq).round(4)
    latex_matrix = r"Z_{012} = \begin{bmatrix}"
    for row in df.values:
        latex_matrix += " & ".join(f"{v:.3f}" for v in row)
        latex_matrix += r" \\ "
        latex_matrix += r"\end{bmatrix}"
    st.latex(latex_matrix)
    # Matriz anterior de impedancias de secuencia pero con argumento y angulo
    print("\nMatriz de impedancias de secuencia con magnitud y ángulo (Ω):")
    for i, fila in enumerate(Z_seq):
        print(f" ", end="")
        for Z in fila:
            magnitud = abs(Z)
            angulo = math.degrees(math.atan2(Z.imag, Z.real))
            print(f"{magnitud:12.4f} ∠ {angulo:8.2f}°", end="  ")
        print()
    Z0 = Z_seq[2, 2]
    print(f"\nImpedancia homopolar de la línea (Z0): {Z0:.3f} Ω/km")
    Z1 = Z_seq[1, 1]
    print(f"Impedancia directa e inversa de la línea (Z1): {Z1:.3f} Ω/km")
    #Teniendo en cuenta la longitud de la línea
    Z_seq_total = Z_seq * longitud
    print("\nMatriz de impedancias de secuencia total para la longitud dada (Ω):")
    for i, fila in enumerate(Z_seq_total):
        print(f" ", end="")
        for Z in fila:
            print(f"{Z:18.3f}", end="  ")
        print()
    # Matriz anterior de impedancias de secuencia pero con argumento y angulo
    print("\nMatriz de impedancias de secuencia total con magnitud y ángulo para la longitud dada (Ω):")
    for i, fila in enumerate(Z_seq_total):
        print(f" ", end="")
        for Z in fila:
            magnitud = abs(Z)
            angulo = math.degrees(math.atan2(Z.imag, Z.real))
            print(f"{magnitud:12.3f} ∠ {angulo:8.2f}°", end="  ")
        print()
    Z0_total = Z_seq_total[2, 2]
    print(f"\nImpedancia homopolar de la línea total (Z0): {Z0_total:.3f} Ω")
    Z1_total = Z_seq_total[1, 1]
    print(f"Impedancia directa e inversa de la línea total (Z1): {Z1_total:.3f} Ω")
    # Resistencias R0/R1
    R0 = Z0.real
    R1 = Z1.real
    print(f"\nR0/R1: {R0/R1:.3f}")
    # Reactancias X0/X1
    X0 = Z0.imag
    X1 = Z1.imag
    print(f"X0/X1: {X0/X1:.3f}")

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # MATRIZ DE CAPACIDADES
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print("\n4. MATRIZ DE CAPACIDADES")
    # Coeficientes de potencial
    # C_ii = (1/(2*pi*epsilon0))*log(D_prima_ii/r_i)
    # C_ij = (1/(2*pi*epsilon0))*log(D_prima_ij/D_ij) para i != j
    matriz_capacidades = [[0]*n_puntos for _ in range(n_puntos)]
    for i in range(n_puntos):
        for j in range(n_puntos):
            D_ij = matriz_distancias[i][j]/1000   # Convertir m a km
            D_prima_ij = matriz_D_prima[i][j]/1000   # Convertir m a km
            r_i = radio_conductores_mm[i] / 1000000  # Convertir mm a km
            if i == j:
                C_ii = ((1 / (2 * pi * epsilon0)) * log(D_prima_ij / r_i))/1000000000
                matriz_capacidades[i][j] = C_ii
            else:
                C_ij = (1 / (2 * pi * epsilon0)) * log(D_prima_ij / D_ij)/1000000000
                matriz_capacidades[i][j] = C_ij
    print("\nMatriz de coeficientes de potencial (km/uF):")
    for i, fila in enumerate(matriz_capacidades):
        print(f" ", end="")
        for C in fila:
            print(f"{C:18.2f}", end="  ")
        print()
    # Matriz de capacitancias por fase (eliminando tierra)
    C = np.array(matriz_capacidades)
    Cf = C[0:6, 0:6]
    Ct = C[6:7, 6:7]
    Cft = C[0:6, 6:7]
    Ctf = C[6:7, 0:6]
    Ct_inv = np.linalg.inv(Ct)
    Cfas = Cf - Cft @ Ct_inv @ Ctf
    print("\nMatriz de coeficientes de potencial por fase (km/μF):")
    for i, fila in enumerate(Cfas):
        print(f" ", end="")
        for C in fila:
            print(f"{C:18.2f}", end="  ")
        print()

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # MATRIZ DE SUSCEPTANCIAS
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print("\n5. MATRIZ DE SUSCEPTANCIAS")
    # Matriz de susceptancias por km (B=j*2*pi*frecuencia*C^-1) solo parte imaginaria, no mostrando la parte real
    Cfas_inv = np.linalg.inv(Cfas)
    Bfas = 1j * 2 * pi * frecuencia * Cfas_inv
    print("\nMatriz de susceptancias por fase (μS/km):")
    for i, fila in enumerate(Bfas):
        print(f"", end="")
        for B in fila:
            print(f"{B:18.3f}", end="  ")
        print()
    # Matriz de susceptancias total para la longitud dada
    Bfas_total = Bfas * longitud
    print("\nMatriz de susceptancias por fase total para la longitud dada (μS):")
    for i, fila in enumerate(Bfas_total):
        print(f" ", end="")
        for B in fila:
            print(f"{B:18.2f}", end="  ")
        print()
    # Susceptancias de secuencia
    B_seq = A_inv @ Bfas @ A
    print("\nMatriz de susceptancias de secuencia (μS):")
    for i, fila in enumerate(B_seq):
        print(f"", end="")
        for B in fila:
            print(f"{B:18.3f}", end="  ")
        print()
    # Susceptancia homopolar de la linea (B0)
    B0 = B_seq[2, 2]
    print(f"\nSusceptancia homopolar de la línea (B0): {B0:.3f} μS/km")
    # Susceptancia directa e inversa de la linea (B1)
    B1 = B_seq[1, 1]
    print(f"Susceptancia directa e inversa de la línea (B1): {B1:.3f} μS/km")
    # Susceptancia homopolar de la linea total (B0)
    B0_total = B_seq[2, 2] * longitud
    print(f"\nSusceptancia homopolar de la línea total (B0): {B0_total:.3f} μS")
    # Susceptancia directa e inversa de la linea total (B1)
    B1_total = B_seq[1, 1] * longitud
    print(f"Susceptancia directa e inversa de la línea total (B1): {B1_total:.3f} μS")

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CÁLCULO IMPEDANCIA CARACTERÍSTICA Y CONSTANTE DE PROPAGACIÓN
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print("\n6. IMPEDANCIA CARACTERÍSTICA Y CONSTANTE DE PROPAGACIÓN")
    # Impedancia característica (Zc=raiz((R+jX)/(jB)))
    Zc = cmath.sqrt(Z1 / (1 * B1*1e-6))
    magnitud_Zc = abs(Zc)
    angulo_Zc = math.degrees(math.atan2(Zc.imag, Zc.real))
    print(f"\nImpedancia característica de secuencia directa (Zc1):{Zc:.3f} = {magnitud_Zc:.4f} ∠ {angulo_Zc:.2f}° Ω")
    # Constante de propagación (gamma=raiz((R+jX)*(jB)))
    gamma = cmath.sqrt(Z1 * (1 * B1*1e-6))
    magnitud_gamma = abs(gamma)
    angulo_gamma = math.degrees(math.atan2(gamma.imag, gamma.real))
    print(f"\nConstante de propagación de secuencia directa (γ1): {gamma:.4e} = {magnitud_gamma:.4e} ∠ {angulo_gamma:.2f}° 1/km")
    # Constante de propagación total para la longitud dada
    gamma_total = gamma * longitud
    magnitud_gamma_total = abs(gamma_total)
    angulo_gamma_total = math.degrees(math.atan2(gamma_total.imag, gamma_total.real))
    print(f"\nConstante de propagación para la longitud dada: {gamma_total:.4e} = {magnitud_gamma_total:.4e} ∠ {angulo_gamma_total:.2f}°")

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # POTENCIA CARACTERÍSTICA
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print("\n7. POTENCIA CARACTERÍSTICA")
    # Pc = tension_nominal**2 / Zc
    Pc = (tension_nominal) ** 2 / Zc
    magnitud_Pc = abs(Pc)
    angulo_Pc = math.degrees(math.atan2(Pc.imag, Pc.real))
    print(f"\nPotencia característica: {Pc:.4f} = {magnitud_Pc:.4f} ∠ {angulo_Pc:.2f}° kVA")

    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CAIDA DE TENSIÓN
    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print("\n8. CAIDA DE TENSIÓN")
    # ΔU=potencia_transportada(MW)*longitud*(R*cos(phi)+X*sin(phi))/((tension_nominal**2)*10*cos(phi))
    potencia_transportada_MW = potencia_transportada * cos_phi  # Convertir MVA a MW usando el factor de potencia
    ΔU = potencia_transportada_MW * 1000 * longitud * (Z1.real * cos_phi + Z1.imag * math.sin(math.acos(cos_phi))) / ((tension_nominal ** 2) * 10 * cos_phi)
    print(f"\nCaída de tensión ΔU: {ΔU:.4f} %")
    # Poner si la caida de tension es inferior al 5%: La caída de tensión es inferior al 5%; y si no lo es: La caída de tensión es superior al 5%
    if ΔU < 5:
        print("La caída de tensión es inferior al 5%")
    else:
        print("La caída de tensión es superior al 5%")

    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # DENSIDAD MÁXIMA DE CORRIENTE E INTENSIDAD MÁXIMA POR CABLE
    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print("\n9. DENSIDAD MÁXIMA DE CORRIENTE E INTENSIDAD MÁXIMA POR CABLE")
    # Valores densidad de corriente según la sección
    secciones = [15, 25, 35, 50, 70, 95, 125, 160, 200, 250, 300, 400, 500, 600]
    densidades = [6, 5, 4.55, 4.0, 3.55, 3.2, 2.9, 2.7, 2.5, 2.3, 2.15, 1.95, 1.8, 1.65]  # A/mm2
    def densidad_corriente(seccion):
        if seccion < secciones[0] or seccion > secciones[-1]:
            return None  # Fuera de rango
        for i in range(len(secciones) - 1):
            if secciones[i] <= seccion <= secciones[i + 1]:
                # Regla de tres
                densidad = densidades[i] + (densidades[i + 1] - densidades[i]) * (seccion - secciones[i]) / (secciones[i + 1] - secciones[i])
                return densidad
        return None
    densidad_maxima = densidad_corriente(seccion)
    print(f"\nDensidad máxima de corriente para la sección {seccion} mm2: {densidad_maxima:.4f} A/mm2")
    # Valores coeficiente reductor según composición material
    composiciones = ["30+7", "6+1", "26+7", "54+7", "45+7"]
    coef_reductores = [0.916, 0.937, 0.937, 0.95, 0.97]
    def coeficiente_reductor(composicion):
        if composicion in composiciones:
            index = composiciones.index(composicion)
            return coef_reductores[index]
        return None
    coef_reductor = coeficiente_reductor(composicion)
    print(f"Coeficiente reductor para la composición {composicion}: {coef_reductor:.4f}")

    densidad_max_con_reduccion = densidad_maxima * coef_reductor
    print(f"Densidad máxima de corriente con coeficiente reductor: {densidad_max_con_reduccion:.4f} A/mm2")
    # Intensidad máxima admisible
    intensidad_maxima_conductor = densidad_max_con_reduccion * seccion  # A
    print(f"Corriente máxima admisible del conductor: {intensidad_maxima_conductor:.2f} A")

    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # POTENCIA MÁXIMA ADMISIBLE POR INTENSIDAD
    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print("\n10. POTENCIA MÁXIMA ADMISIBLE POR INTENSIDAD")
    potencia_maxima_admisible = (intensidad_maxima_conductor * tension_nominal * raiz(3)) * cos_phi 
    print(f"\nPotencia máxima admisible por intensidad: {potencia_maxima_admisible/1000:.4f} MW")
    # Comparar potencia máxima admisible con potencia transportada
    if potencia_maxima_admisible/1000 >= potencia_transportada_MW:
        print("La potencia máxima admisible por intensidad es mayor que la potencia transportada.")
    else:
        print("La potencia máxima admisible por intensidad es menor que la potencia transportada.")

    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # POTENCIA DE TRANSPORTE EN FUNCIÓN DE CONDICIONES METEOROLÓGICAS
    #-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print("\n11. CÁLCULO POTENCIA DE TRANSPORTE EN FUNCIÓN DE CONDICIONES METEOROLÓGICAS")
    # Calor aportado por radiación solar para invierno y verano Qs = coeficiente_absorcion * radiacion * diametro/1000
    print("\n   11.1. CALOR APORTADO POR RADIACIÓN SOLAR")
    Qs_invierno = coeficiente_absorcion * radiacion_invierno * diametro/1000 # W/m
    Qs_verano = coeficiente_absorcion * radiacion_verano * diametro/1000 # W/m
    print(f"\nCalor aportado por radiación solar en invierno: {Qs_invierno:.3f} W/m")
    print(f"Calor aportado por radiación solar en verano: {Qs_verano:.3f} W/m")
    # Calor cedido por radiación Qr = pi * diametro/1000 * emisividad_conductor * cte_boltzmann * ((Tc+273.15)^4 - (Ta+273.15)^4)
    print("\n   11.2. CALOR CEDIDO POR RADIACIÓN SOLAR")
    cte_boltzmann = 5.6704e-8 # W/m2K4
    Qr_invierno = pi * (diametro/1000) * emisividad_conductor * cte_boltzmann * ((Tc + 273.15)**4 - (temperatura_invierno + 273.15)**4) # W/m
    Qr_verano = pi * (diametro/1000) * emisividad_conductor * cte_boltzmann * ((Tc + 273.15)**4 - (temperatura_verano + 273.15)**4) # W/m
    print(f"\nCalor cedido por radiación en invierno: {Qr_invierno:.3f} W/m")
    print(f"Calor cedido por radiación en verano: {Qr_verano:.3f} W/m")

    # Calor cedido por convección Qc = pi * conductividad_termica * (Tc - Ta) * Nu
    print("\n   11.3. CALOR CEDIDO POR CONVECCIÓN")
    # 1.  Convección natural: Nu = A * (Gr*Pr)^m
    print("\n       11.3.1. CONVECCIÓN NATURAL")
    """ Valores para A y m dependiendo del rango de Gr*Pr:
    0.1 <= Gr*Pr < 100: A=1.02, m=0.148
    100 <= Gr*Pr < 10000: A=0.85, m=0.188
    10000 <= Gr*Pr < 1e7: A=0.48, m=0.25
    1e7 <= Gr*Pr < 1e12: A=0.125, m=0.333
    """
    def coeficientes_conveccion_natural(Gr_Pr):
        if 0.1 <= Gr_Pr < 100:
            return 1.02, 0.148
        elif 100 <= Gr_Pr < 10000:
            return 0.85, 0.188
        elif 10000 <= Gr_Pr < 1e7:
            return 0.48, 0.25
        elif 1e7 <= Gr_Pr < 1e12:
            return 0.125, 0.333
        else:
            return None, None
    # Gr = (diametro/1000)^3 * (Tc - Ta) * g / ((Tav + 27.15) * (viscosidad_cinematica)^2)
    g = 9.81  # m/s2
    # viscosidad_cinematica = 1.32e-5 + 9.5e-8 * Tav
    def viscosidad_cinematica(Tav):
        return 1.32e-5 + 9.5e-8 * Tav
    viscosidad_cinematica_invierno = viscosidad_cinematica((Tc + temperatura_invierno) / 2)
    viscosidad_cinematica_verano = viscosidad_cinematica((Tc + temperatura_verano) / 2)
    print(f"\nViscosidad cinemática en invierno: {viscosidad_cinematica_invierno:.3e} m²/s")
    print(f"Viscosidad cinemática en verano: {viscosidad_cinematica_verano:.3e} m²/s")
    # densidad_relativa_aire = exp(-1.16e-4* h)
    def Gr(Ta, viscosidad_cinematica):
        Tav = (Tc + Ta) / 2
        return ((diametro/1000)**3 * (Tc - Ta) * g) / ((Tav + 273.15) * (viscosidad_cinematica**2))
    Gr_invierno = Gr(temperatura_invierno, viscosidad_cinematica_invierno)
    Gr_verano = Gr(temperatura_verano, viscosidad_cinematica_verano)
    print(f"\nNúmero de Grashof en invierno: {Gr_invierno:.2f}")
    print(f"Número de Grashof en verano: {Gr_verano:.2f}")
    # Pr = calor_esp_aire * viscosidad_dinamica / conductividad_termica_aire
    calor_esp_aire = 1005  # J/kgK
    # viscosidad_dinamica = (17.239 + 4.625e-2 * Tav - 2.03e-5 * Tav **2)e-6 kg/m s
    def viscosidad_dinamica(Tav):
        return (17.239 + 4.625e-2 * Tav - 2.03e-5 * Tav **2) * 1e-6  # kg/m s
    viscosidad_dinamica_invierno = viscosidad_dinamica((Tc + temperatura_invierno) / 2)
    viscosidad_dinamica_verano = viscosidad_dinamica((Tc + temperatura_verano) / 2)
    print(f"\nViscosidad dinámica en invierno: {viscosidad_dinamica_invierno:.3e} kg/m s")
    print(f"Viscosidad dinámica en verano: {viscosidad_dinamica_verano:.3e} kg/m s")
    # conductividad_termica_aire = 2.368e-2 + 7.23e-5 * Tav - 2.763e-8 * Tav**2 W/Km
    def conductividad_termica_aire(Tav):
        return 2.368e-2 + 7.23e-5 * Tav - 2.763e-8 * Tav**2  # W/Km
    conductividad_termica_aire_invierno = conductividad_termica_aire((Tc + temperatura_invierno) / 2)
    conductividad_termica_aire_verano = conductividad_termica_aire((Tc + temperatura_verano) / 2)
    print(f"\nConductividad térmica del aire en invierno: {conductividad_termica_aire_invierno:.3e} W/K·m")
    print(f"Conductividad térmica del aire en verano: {conductividad_termica_aire_verano:.3e} W/K·m")

    Pr_invierno = calor_esp_aire * viscosidad_dinamica_invierno / conductividad_termica_aire_invierno
    Pr_verano = calor_esp_aire * viscosidad_dinamica_verano / conductividad_termica_aire_verano
    print(f"\nNúmero de Prandtl en invierno: {Pr_invierno:.4f}")
    print(f"Número de Prandtl en verano: {Pr_verano:.4f}")
    Gr_Pr_invierno = Gr_invierno * Pr_invierno
    Gr_Pr_verano = Gr_verano * Pr_verano
    Nu_invierno = 0
    Nu_verano = 0
    A_invierno, m_invierno = coeficientes_conveccion_natural(Gr_Pr_invierno)
    if A_invierno is not None:
        Nu_invierno = A_invierno * (Gr_Pr_invierno ** m_invierno)
    A_verano, m_verano = coeficientes_conveccion_natural(Gr_Pr_verano)
    if A_verano is not None:
        Nu_verano = A_verano * (Gr_Pr_verano ** m_verano)
    print(f"\nNúmero de Nusselt en invierno: {Nu_invierno:.4f}")
    print(f"Número de Nusselt en verano: {Nu_verano:.4f}")

    # 2. Convección forzada: Nu = B1 * Re^n
    print("\n       11.3.2. CONVECCIÓN FORZADA")
    # Re = (diametro/1000) * velocidad_viento / viscosidad_cinematica
    # Rugosidad: Rf = diametro_alambre_ext / (2*(diametro-diametro_alambe_ext))
    Rf = diametro_alambre_ext / (2 * (diametro - diametro_alambre_ext))
    print(f"\nRugosidad Rf: {Rf:.3f}")

    def Re(viscosidad_cinematica):
        return (diametro/1000) * velocidad_viento / viscosidad_cinematica
    Re_invierno = Re(viscosidad_cinematica_invierno)
    Re_verano = Re(viscosidad_cinematica_verano)
    print(f"\nNúmero de Reynolds en invierno: {Re_invierno:.2f}")
    print(f"Número de Reynolds en verano: {Re_verano:.2f}")
    # B1 y n según Re
    def coeficientes_conveccion_forzada(Re, Rf):
        if 100 < Re < 2.65e3:
            return 0.641, 0.471
        elif Rf <= 0.05 and 2.65e3 <= Re < 5e4:
            return 0.178, 0.633
        elif Rf > 0.05 and 2.65e3 <= Re < 5e4:
            return 0.048, 0.8
        else:
            return None, None
    Nu_invierno_forzada = 0
    Nu_verano_forzada = 0
    B1_invierno, n_invierno = coeficientes_conveccion_forzada(Re_invierno, Rf)
    if B1_invierno is not None:
        Nu_invierno_forzada = B1_invierno * (Re_invierno ** n_invierno)
    B1_verano, n_verano = coeficientes_conveccion_forzada(Re_verano, Rf)
    if B1_verano is not None:
        Nu_verano_forzada = B1_verano * (Re_verano ** n_verano)
    print(f"\nNúmero de Nusselt por convección forzada en invierno: {Nu_invierno_forzada:.4f}")
    print(f"Número de Nusselt por convección forzada en verano: {Nu_verano_forzada:.4f}")
    Nu_45_invierno = Nu_invierno_forzada * (0.42 + 0.58*math.sin(45*pi/180)**0.9)
    Nu_45_verano = Nu_verano_forzada * (0.42 + 0.58*math.sin(45*pi/180)**0.9)
    print(f"\nNúmero de Nusselt corregido para ángulo de 45° en invierno: {Nu_45_invierno:.4f}")
    print(f"Número de Nusselt corregido para ángulo de 45° en verano: {Nu_45_verano:.4f}")
    # Qc = pi* conductividad_termica * (Tc-Ta) * Nu con Nu=max(Nu_natural, Nu_45_forzada)
    Qc_invierno = pi * conductividad_termica_aire_invierno * (Tc - temperatura_invierno) * max (Nu_invierno, Nu_45_invierno)
    Qc_verano = pi * conductividad_termica_aire_verano * (Tc - temperatura_verano) * max (Nu_verano, Nu_45_verano)
    print(f"\nCalor cedido por convección en invierno: {Qc_invierno:.3f} W/m")
    print(f"Calor cedido por convección en verano: {Qc_verano:.3f} W/m")

    # Resultados corriente máxima: I = raiz((Qr+Qc-Qs)/resistencia_ca/1000)
    print("\n   11.4. RESULTADOS CORRIENTE MÁXIMA")
    print(f"\nResistencia del conductor a la temperatura de cálculo: {resistencia_ca:.6f} Ω/km")
    I_max_invierno = raiz((Qr_invierno + Qc_invierno - Qs_invierno)/ (resistencia_ca * 1e-3))
    I_max_verano = raiz((Qr_verano + Qc_verano - Qs_verano) / (resistencia_ca * 1e-3))
    print(f"\nCorriente máxima admisible en invierno según condiciones meteorológicas: {I_max_invierno:.2f} A")
    print(f"Corriente máxima admisible en verano según condiciones meteorológicas: {I_max_verano:.2f} A")

    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # POTENCIA MÁXIMA DE TRANSPORTE SEGÚN CONDICIONES METEOROLÓGICAS
    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print("\n   11.5. POTENCIA MÁXIMA DE TRANSPORTE")
    potencia_maxima_invierno = (I_max_invierno * tension_nominal * raiz(3)) * cos_phi
    potencia_maxima_verano = (I_max_verano * tension_nominal * raiz(3)) * cos_phi
    print(f"\nPotencia máxima de transporte en invierno según condiciones meteorológicas: {potencia_maxima_invierno/1000:.2f} MW")
    print(f"Potencia máxima de transporte en verano según condiciones meteorológicas: {potencia_maxima_verano/1000:.2f} MW")

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # PÉRDIDAS DE POTENCIA
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print("\n12. PÉRDIDAS DE POTENCIA")
    perdidas_potencia = ((potencia_transportada_MW * resistencia_ca * longitud) / (tension_nominal**2 * cos_phi**2))*100
    print(f"\nPérdidas de potencia en la línea: {perdidas_potencia:.5f} %")
    # En valor absoluto
    perdidas_potencia_valor = (potencia_transportada_MW * resistencia_ca * longitud) / (tension_nominal**2 * cos_phi**2) * potencia_transportada_MW
    print(f"Pérdidas de potencia en la línea en valor absoluto: {perdidas_potencia_valor:.5f} MW")

    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CORTOCIRCUITO MÁXIMO
    #----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print("\n13. CORTOCIRCUITO MÁXIMO")
    # Contantes para cálculo de factor K en función del material del conductor
    materiales = ["Cobre", "Aluminio-Acero", "Acero"]
    c_conductores = [390, 910, 480]  # J/kgK
    densidades_conductores = [8900, 2700, 7850]  # kg/m3
    ctes_conductores = [56e6, 34.8e6, 7.25e6]  # 1/ohmio m
    alphas_conductores = [0.0039, 0.004, 0.0045]  # 1/ºC
    def propiedades_material(material):
        if material in materiales:
            index = materiales.index(material)
            return (c_conductores[index], densidades_conductores[index], ctes_conductores[index], alphas_conductores[index])
        return (None, None, None, None)
    c_conductor, densidad_conductor, cte_conductor, alpha_conductor = propiedades_material(material)
    print(f"c_conductor, densidad_conductor, cte_conductor, alpha_conductor", c_conductor, densidad_conductor, cte_conductor, alpha_conductor )
    # Temperatura máxima recomendada en función del material (ºC)
    if material == "Acero":
        temperatura_max_recomendada = 300
    else:
        temperatura_max_recomendada = 200
    # Factor K
    multiplicacion = log ((1+alpha_conductor*(temperatura_max_recomendada-20))/(1+alpha_conductor*(Tc-20)))
    K = (raiz (cte_conductor*c_conductor*densidad_conductor*multiplicacion/alpha_conductor))/1000000
    print(f"Factor K: {K:.4f} Araiz(s)/mm2")
    # Corriente máxima de cortocircuito
    Icc_max = K * seccion / raiz(tiempo_accionamiento_proteccion)
    print(f"Icc max: {Icc_max:.4f} A")

    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # EFECTO CORONA
    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print("\n14. EFECTO CORONA ")
    # Presión barométrica
    presion_barometrica = 1/(log(log(76)-altitud_media/18330))*100 
    print(f"\nPresión barométrica h = {presion_barometrica:.3f} cmHg")
    # Factor corrección densidad aire: delta
    def delta (temperatura):
        return 3.92*presion_barometrica/(273+temperatura)
    delta_invierno = delta (temperatura_invierno)
    delta_verano = delta (temperatura_verano)
    print(f"\nFactor corrección densidad aire invierno (δinv) = {delta_invierno:.4f}")
    print(f"Factor corrección densidad aire verano (δver) = {delta_verano:.4f}")


    DMG = 100*(matriz_distancias[0][1]*matriz_distancias[1][2]*matriz_distancias[2][0])**(1/6)
    print(f"DMG: {DMG}")
    DMG = 807
    print(f"DMG: {DMG}")

    # Tensión crítica dieléctrica
    def Uc (mt,delta):
        return raiz(3)*mc*mt*(30/raiz(2))*delta*(diametro/20)*log(DMG/(diametro/20))
    Uc_invierno = Uc(mt_invierno, delta_invierno)
    Uc_verano = Uc(mt_verano,delta_verano)
    print(f"\nTensión critica invierno {Uc_invierno:.4f} kV")
    print(f"Tensión crítica verano {Uc_verano:.4f} kV")

    # Pérdidas por efecto corona
    if Uc_invierno < tension_nominal:
        Perdidas_efecto_corona_invierno = 241/delta_invierno * (frecuencia + 25) * raiz (diametro/20/DMG) * (tension_nominal - Uc_invierno)/raiz(3) * 1e-5
        print(f"\nHay pérdidas por efecto corona en invierno: {Perdidas_efecto_corona_invierno:.4f} kW·km/fase")
    else: print(f"\nNo hay pérdidas por efecto corona en invierno")

    if Uc_verano < tension_nominal:
        Perdidas_efecto_corona_verano = 241/delta_verano * (frecuencia + 25) * raiz (diametro/20/DMG) * (tension_nominal - Uc_verano)/raiz(3) * 1e-5
        print(f"Hay pérdidas por efecto corona en verano: {Perdidas_efecto_corona_verano:.4f} kW·km/fase")
    else: print(f"No hay pérdidas por efecto corona en verano")


    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # CAMPO ELÉCTRICO
    #------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    print("\n15. CAMPO ELÉCTRICO ")
    # Matriz coeficientes de potencial
    print("\nMatriz de capacidades (km/uF):")
    for i, fila in enumerate(matriz_capacidades):
        print(f"", end="")
        for C in fila:
            print(f"{C:18.2f}", end="  ")
        print()       
    # Matriz coeficientes de potencial inversa
    matriz_capacidades_inversa = np.linalg.inv(matriz_capacidades) * 1000
    print("\nMatriz de capacidades inversa(nF/km):")
    for i, fila in enumerate(matriz_capacidades_inversa):
        print(f"", end="")
        for C in fila:
            print(f"{C:18.2f}", end="  ")
        print()
    # Vector de potenciales (tensiones nominales/raiz(3) (angulos 0 -120 y 120, tension tierra)) en numeros complejos y kV
    V_phase = tension_nominal / raiz(3)
    V_vector = np.array([V_phase * cmath.rect(1, math.radians(0)),
                        V_phase * cmath.rect(1, math.radians(-120)),
                        V_phase * cmath.rect(1, math.radians(120)),
                        V_phase * cmath.rect(1, math.radians(0)),
                        V_phase * cmath.rect(1, math.radians(-120)),
                        V_phase * cmath.rect(1, math.radians(120)),
                        0])
    print("\nVector de potenciales (kV):")
    for i, V in enumerate(V_vector):
        print(f"{V:.3f}")
    # Vector de cargas (matriz_capacidades_inversa @ V_vector)
    Q_vector = matriz_capacidades_inversa @ V_vector /1000
    print("\nVector de cargas (kV·nF/km):")
    for i, Q in enumerate(Q_vector):
        print(f"{Q:.3f} mC/km")




    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    # EJECUCIÓN PROGRAMA
    #---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    sys.stdout = sys.__stdout__
    texto_resultados = buffer.getvalue()
    html = """
    <html>
    <head>
    <meta charset="utf-8">
    <style>
    body { font-family: Arial; }
    .mayus { font-weight: bold; }
    </style>
    </head>
    <body>
    <pre>
    """
    def es_mayusculas(linea):
        letras = [c for c in linea if c.isalpha()]
        return len(letras)>0 and all(c.isupper() for c in letras)

    for linea in texto_resultados.split("\n"):
        # Detectar líneas completamente en mayúsculas
        if es_mayusculas(linea):
            html += f'<span class="mayus">{linea}</span>\n'
        else:
            html += linea + "\n"
    html += """
    </pre>
    </body>
    </html>
    """
    with open("resultados.html", "w", encoding="utf-8") as f:
        f.write(html)






    from docx import Document

    doc = Document()
    doc.add_heading("Resultados", 1)

    for linea in texto_resultados.split("\n"):
        doc.add_paragraph(linea)

    doc.save("resultados.docx")




#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# SIMPLE CIRCUITO DUPLEX
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
