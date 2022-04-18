from matplotlib.pyplot import plot
import scipy.io
import control

import scipy as sp
import numpy as np
import scipy.linalg as la

from math import *
from cmath import phase
from control.matlab import *


def d2c(sys):
    sys = tf2ss(sys)

    a = sys.A
    b = sys.B
    c = sys.C
    d = sys.D
    Ts = sys.dt
    n = np.shape(a)[0]
    nb = np.shape(b)[1]

    if n == 1 and a[0, 0] == 1:
        A = 0
        B = b/sys.dt
        C = c
        D = d
    else:
        tmp1 = np.hstack((a, b))
        tmp2 = np.hstack((np.zeros((nb, n)), np.eye(nb)))
        tmp = np.vstack((tmp1, tmp2))
        s = la.logm(tmp)
        s = s/Ts
        if la.norm(np. imag(s), np.inf) > np.sqrt(sp.finfo(float).eps):
            print('Warning: accuracy may be poor')
        s = np.real(s)
        A = s[0:n, 0:n]
        B = s[0:n, n:n+nb]
        C = c
        D = d

    sysc = StateSpace(A, B, C, D)
    return ss2tf(sysc)


def mean_square(x1, y1):
    u = x1
    y = y1

    x = len(u)
    N = x
    M = N - 1

    T = 0.1

    F = np.block([y[0:M], u[0:M]])
    Y = np.block(y[1:N])

    theta = np.linalg.inv(np.transpose(F) @ F) @ np.transpose(F) @ Y
    a1 = theta[0][0]
    b1 = theta[1][0]
    sysz = tf([b1], [1, -a1], T)

    mesh_open_tf = d2c(sysz)
    return a1, b1, mesh_open_tf


def rf_pid_sintonization(num, den, mp, ta):

    Num = num  # Numerador da F.T.
    Den = den  # Denominador da F.T.

    Mp = mp  # Máximo Pico de 15#
    Ta = ta  # Tempo de Acomodação

    Ess = 0
    Ki = 0
    Kp = 0
    Kd = 0
    # ________________ Cálculo Ki_______________________
    if Ess != 0:
        Kv = 1/Ess
        Ki = Kv/Num(1)
    else:
        Kv = 0
    # __________________________________________

    Mp1 = Mp/100
    Mp2 = log(Mp1)/(-pi)

    Qsi = sqrt((Mp2 ** 2)/(1+Mp2 ** 2))  # Cálculo de Qsi

    MF = degrees((asin(Qsi)))*2

    Wn = 4/(Qsi*Ta)  # Frequência Natural não amortecida
    Wcg = 1j * Wn  # Frequência de Cruzamento de Ganho

    G_jwcg = Num/((Den[0]*Wcg) + 1)

    Mod_Gjwcg = abs(G_jwcg)

    Anglo_Gjwcg = (phase(G_jwcg)*180)/pi

    Theta = -180 + MF - Anglo_Gjwcg

    if Kv != 0:
        Kp = cos(radians(Theta)) / Mod_Gjwcg
        Kd = (Ki/Wn ** 2) + (sin(Theta)/(Mod_Gjwcg*Wn))

    else:
        Kp = cos(radians(Theta)) / Mod_Gjwcg
        Ki = -(sin(radians(Theta))*Wn ** 2)/(Mod_Gjwcg * Wn)
        Kd = 0

    return Kp, Ki, Kd


def output_close_mesh_pi(a1, b1, kp, ki, sp):
    # Aplicando o algoritmo para verificar a saída do sistema em malha fechada com controlador PI
    a1 = a1
    b1 = b1
    T = 0.1
    M = 250

    # Definindo pv como sendo a saída do sistema
    pv = [0] * (M-1)

    # Definindo cont como a saída do controlador PI
    cont = [0] * (M-1)

    # Definindo a ação proporcional do sitema
    P = [0] * (M-1)

    # Definindo a ação integral do sistema
    I = [0] * (M-1)
    sp = sp

    # Definindo o ganho proporcional do sistema
    Kp = kp

    # Definindo o ganho integral do sistema
    Ki = ki

    erro = [0] * (M-1)
    t = [0] * (M-1)
    # Implementando a malha fechada do sistema com controlador PI
    # Por meio da utilização da equação a diferenças do sistema

    for i in range(1, M-1, 1):
        pv[i] = a1*pv[i-1]+b1*cont[i-1]  # Saída instantanea do sistema
        erro[i] = sp - pv[i]  # Erro instataneo do sistema
        P[i] = Kp*erro[i]  # Ação proporcioanl
        I[i] = I[i-1]+Ki*erro[i]*T  # Ação Integral
        cont[i] = P[i]+I[i]  # Ação do controlador PI saida
        t[i] = T*i  # variável de tempo

    # Plotando os gráficos de entrada e saída do sistema
    return t, pv


def read_matlab_data(file):
    """
    returns x1, y1, T
    """
    mat = scipy.io.loadmat(file)

    return mat['x1'], mat['y1'], mat['T'][0]


def system_info(tf, t):
    S = stepinfo(tf, t)
    return S['Overshoot'], S['SettlingTime'], S['RiseTime']
