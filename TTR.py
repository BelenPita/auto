
import math 
import cmath
import matplotlib.pyplot as plt
import numpy as np


#--------------------------------------------------------------------------------------------------------------------------------------------
# APARTADO 1 
#--------------------------------------------------------------------------------------------------------------------------------------------

# DATOS
U0 = 107.8 #kV
L = 20e-3 #He
C = 0.1e-6 #F
w = 100 * math.pi #rad/s
w1 = 1 / math.sqrt(L * C)

# APARTADO A: EXPRESIÓN ANALÍTICA TTR
t = np.linspace(0, 0.02, 1000)
U_s = U0 * np.cos(w * t)
U_TTR = U0 * np.cos(w * t) - U0 * np.cos(w1 * t)

plt.plot(t, U_s, label='U_s')
plt.plot(t, U_TTR, label='U_TTR')
plt.xlim(0, 0.02)
plt.xlabel('Tiempo (s)')
plt.ylabel('Tensión (kV)')
plt.title('Tensión TTR')
plt.legend()
plt.grid(True)

# APARTADO B: RESISTENCIA DE PREINSERCIÓN y TTR
R_c = 0.5 * math.sqrt(L / C)
print(f'Resistencia de preinserción R_c = {R_c:.2f} Ohmios')

Y = w*C*1j - 1j/(w*L) + 1/R_c
U = - U0/(w*L) * 1j * (1/Y)
r = -1 / (2*R_c*C)
A = - abs(U)
B = - A * r + abs(U) * w * math.sin(np.angle(U))

t = np.arange(0, 0.02, 1e-6)
U_TTR2 = abs(U) * np.cos(w*t - np.angle(U)) + (A + B * t) * np.exp(r*t)
U_s = U0 * np.cos(w * t)

plt.figure()
plt.plot(t, U_TTR2, label='U_TTR')
plt.plot(t, U_s, label='U_s')
plt.xlim(0, 0.02)
plt.xlabel('Tiempo (s)')
plt.ylabel('Tensión (kV)')
plt.title('U_TTR con resistencia de preinserción')
plt.legend()
plt.grid(True)
plt.show()

# APARTADO C: CON 2R_c, 4R_c y 10R_c
# 2Rc
R_s = 2 * R_c
Y = w*C*1j - 1j/(w*L) + 1/R_s
U = - U0/(w*L) * 1j * (1/Y)

alpha = - 1 / (2*R_s*C)
beta = math.sqrt(4/(L*C) - 1/(R_s**2 * C**2)) / 2
A = - abs(U) * np.cos(np.angle(U))
B = - abs(U) * w * (math.sin(np.angle(U)) + alpha*A) / beta

t = np.arange(0, 0.02, 1e-6)
U_TTR3 = abs(U) * np.cos(w*t - np.angle(U)) + (A*np.cos(beta*t) + B * np.sin(beta*t)) * np.exp(alpha*t)
U_s = U0 * np.cos(w * t)

plt.figure()
plt.plot(t, U_TTR3, label='U_TTR')
plt.plot(t, U_s, label='U_s')
plt.xlim(0, 0.02)
plt.ylim(-200, 200)
plt.xlabel('Tiempo (s)')
plt.ylabel('Tensión (kV)')
plt.title('U_TTR con 2Rc')
plt.legend()
plt.grid(True)

# 4Rc
R_s = 4 * R_c
Y = w*C*1j - 1j/(w*L) + 1/R_s
U = - U0/(w*L) * 1j * (1/Y)

alpha = - 1 / (2*R_s*C)
beta = math.sqrt(4/(L*C) - 1/(R_s**2 * C**2)) / 2
A = - abs(U) * np.cos(np.angle(U))
B = - abs(U) * w * (math.sin(np.angle(U)) + alpha*A) / beta

t = np.arange(0, 0.02, 1e-6)
U_TTR4 = abs(U) * np.cos(w*t - np.angle(U)) + (A*np.cos(beta*t) + B * np.sin(beta*t)) * np.exp(alpha*t)
U_s = U0 * np.cos(w * t)

plt.figure()
plt.plot(t, U_TTR4, label='U_TTR')
plt.plot(t, U_s, label='U_s')
plt.xlim(0, 0.02)
plt.ylim(-200, 200)
plt.xlabel('Tiempo (s)')
plt.ylabel('Tensión (kV)')
plt.title('U_TTR con 4Rc')
plt.legend()
plt.grid(True)

# 10Rc
R_s = 10 * R_c
Y = w*C*1j - 1j/(w*L) + 1/R_s
U = - U0/(w*L) * 1j * (1/Y)

alpha = - 1 / (2*R_s*C)
beta = math.sqrt(4/(L*C) - 1/(R_s**2 * C**2)) / 2
A = - abs(U) * np.cos(np.angle(U))
B = - abs(U) * w * (math.sin(np.angle(U)) + alpha*A) / beta

t = np.arange(0, 0.02, 1e-6)
U_TTR5 = abs(U) * np.cos(w*t - np.angle(U)) + (A*np.cos(beta*t) + B * np.sin(beta*t)) * np.exp(alpha*t)
U_s = U0 * np.cos(w * t)

plt.figure()
plt.plot(t, U_TTR5, label='U_TTR')
plt.plot(t, U_s, label='U_s')
plt.xlim(0, 0.02)
plt.ylim(-200, 200)
plt.xlabel('Tiempo (s)')
plt.ylabel('Tensión (kV)')
plt.title('U_TTR con 10Rc')
plt.legend()
plt.grid(True)
plt.show()

# APARTADO D: TTR CON BOBINA REAL DE 1 OHMIO. FACTOR DE AMPLITUD
R_l = 1 #Ohmio

Icc = U0 / math.sqrt((w*L)**2 + R_l**2) * (np.cos(-math.atan(w*L/R_l)) + np.sin(-math.atan(w*L/R_l))*1j)
U = ((w*L*1j + R_l) * (-1j / (w*C))) / (w*L*1j + R_l - 1j/(w*C)) * Icc



alpha = - R_l / (2*L)
beta = math.sqrt(4/(L*C) - (R_l**2/L**2)) / 2
A = - abs(U) * np.cos(np.angle(U))
B = - (abs(U) * w * math.sin(np.angle(U)) + alpha*A) / beta

t = np.arange(0, 0.02, 1e-6)
U_TTR6 = abs(U) * np.cos(w*t - np.angle(U)) + (A*np.cos(beta*t) + B * np.sin(beta*t)) * np.exp(alpha*t)
U_s = U0 * np.cos(w * t)

factor_amplitud = np.max(abs(U_TTR6)) / U0
print(f"Factor de amplitud = {factor_amplitud:.4f}")

plt.figure()
plt.plot(t, U_TTR6, label='U_TTR')
plt.plot(t, U_s, label='U_s')
plt.xlim(0, 0.02)
plt.ylim(-250, 250)
plt.xlabel('Tiempo (s)')
plt.ylabel('Tensión (kV)')
plt.title('U_TTR con bobina real de 1 Ohmio')
plt.legend()
plt.grid(True)
plt.show()

# ---------------------------------------------------------------------------------------------------------------------------------------------
# APARTADO 2
# ---------------------------------------------------------------------------------------------------------------------------------------------

# DATOS
X1 = 5.8572 #Ohmios/km
B1 = 0.19143e-4 #S/km
Z1 = 553.23 #Ohmios
v = 0.29666e9 #m/s

# TTR SI FALTA TRIFÁSICA A 2 km, 10 km, 20 km y 90 km
L1 = X1 / w
C1 = B1 / w

# 2 km
l = 2
tau = l/v/1000

U10 = U0 * (L1 * l) / (L + L1 *l)
U1 = (U0 - U10) * (np.cos(w*t) - np.cos(w1*t))

x = np.arange(0, 0.02, 1e-6)
y = [-2*U10, 0]
i = 0

while len(y) < len(x):
    if i % 2 == 0:
        y.append(0)
    else:
        y.append(-2*U10)
    i += 1
U2 = np.array(y)


U_TTR7 = U1 - U2

plt.figure()
plt.plot(t, U_TTR7, label='U_TTR')
plt.plot(t, y, label='U_s')
plt.xlim(0, 0.02)
plt.ylim(-250, 250)
plt.xlabel('Tiempo (s)')
plt.ylabel('Tensión (kV)')
plt.title('U_TTR a 2 km')
plt.legend()
plt.grid(True)


# 2 km
l = 10

U10 = U0 * (L1 * l) / (L + L1 *l)
U1 = (U0 - U10) * (np.cos(w*t) - np.cos(w1*t))
U2 = - 2 * U10

U_TTR8 = U1 - U2

plt.figure()
plt.plot(t, U_TTR8, label='U_TTR')
plt.plot(t, U_s, label='U_s')
plt.xlim(0, 0.02)
plt.ylim(-250, 250)
plt.xlabel('Tiempo (s)')
plt.ylabel('Tensión (kV)')
plt.title('U_TTR a 10 km')
plt.legend()
plt.grid(True)


# 20 km
l = 20

U10 = U0 * (L1 * l) / (L + L1 *l)
U1 = (U0 - U10) * (np.cos(w*t) - np.cos(w1*t))
U2 = - 2 * U10

U_TTR9 = U1 - U2

plt.figure()
plt.plot(t, U_TTR9, label='U_TTR')
plt.plot(t, U_s, label='U_s')
plt.xlim(0, 0.02)
plt.ylim(-250, 250)
plt.xlabel('Tiempo (s)')
plt.ylabel('Tensión (kV)')
plt.title('U_TTR a 20 km')
plt.legend()
plt.grid(True)

# 90 km
l = 90

U10 = U0 * (L1 * l) / (L + L1 *l)
U1 = (U0 - U10) * (np.cos(w*t) - np.cos(w1*t))
U2 = - 2 * U10

U_TTR10 = U1 - U2

plt.figure()
plt.plot(t, U_TTR10, label='U_TTR')
plt.plot(t, U_s, label='U_s')
plt.xlim(0, 0.02)
plt.ylim(-250, 250)
plt.xlabel('Tiempo (s)')
plt.ylabel('Tensión (kV)')
plt.title('U_TTR a 90 km')
plt.legend()
plt.grid(True)
plt.show()


#--------------------------------------------------------------------------------------------------------------------------------------------
# APARTADO 3: falta trifásica punto LX
#--------------------------------------------------------------------------------------------------------------------------------------------
# APARTADO A: TTR segunda lína 90 km
R_s = Z1
Y = w*C*1j - 1j/(w*L) + 1/R_s
U = - U0/(w*L) * 1j * (1/Y)

alpha = - 1 / (2*R_s*C)
beta = math.sqrt(4/(L*C) - 1/(R_s**2 * C**2)) / 2
A = - abs(U) * np.cos(np.angle(U))
B = - abs(U) * w * (math.sin(np.angle(U)) + alpha*A) / beta

t = np.arange(0, 0.02, 1e-6)
U_TTR11 = abs(U) * np.cos(w*t - np.angle(U)) + (A*np.cos(beta*t) + B * np.sin(beta*t)) * np.exp(alpha*t)
U_s = U0 * np.cos(w * t)

plt.figure()
plt.plot(t, U_TTR11, label='U_TTR')
plt.plot(t, U_s, label='U_s')
plt.xlim(0, 0.02)
plt.ylim(-500, 500)
plt.xlabel('Tiempo (s)')
plt.ylabel('Tensión (kV)')
plt.title('U_TTR con línea de 90 km')
plt.legend()
plt.grid(True)
plt.show()

# APARTADO B: TTR tramo lína 120 m entre barras e interruptor
