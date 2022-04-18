from cProfile import label
from ctypes import alignment
import matplotlib.pyplot as plt

from tkinter import *
from tkinter import messagebox
from py_control import *

root = Tk()
root.title("Controlador PID")
root.geometry("300x200")

entries = [Entry(root) for _ in range(5)]
system_info_var = StringVar(value=0)


def on_press():
    inputs = []
    for entry in entries:
        if not (data := entry.get()):
            messagebox.showinfo("Erro", "Preencha todos os campos!")
            break
        inputs.append(data)
    plot(inputs)


def plot(data: list[str]):
    x1, y1, T = read_matlab_data("dadosgrupo7.mat")

    sp = int(data[0])  # set-point
    mp = int(data[1])  # maximo pico
    ta = int(data[2])  # tempo de acomodacao
    kp_acresc = float(data[3])  # kp acrescimento
    ki_acresc = float(data[4])  # ki acrescimento

    a1, b1, mesh_open_tf = mean_square(x1, y1)
    print("malha_aberta: ", mesh_open_tf)

    mesh_close_tf = feedback(mesh_open_tf, 1)
    print("malha_fechada: ", mesh_close_tf)

    den_0 = mesh_open_tf.den[0][0][0]
    den_1 = mesh_open_tf.den[0][0][1]
    num = mesh_open_tf.num[0][0][0]

    canonical_mesh_open_tf = tf([num / den_1], [den_0 / den_1, 1])
    print("malha_fechada_canonica: ", canonical_mesh_open_tf)

    Num = canonical_mesh_open_tf.num[0][0][0]
    Den = canonical_mesh_open_tf.den[0][0]
    kp, ki, kd = rf_pid_sintonization(
        Num,
        Den,
        mp,
        ta)

    kp = kp+kp_acresc
    ki = ki+ki_acresc

    t, pv = output_close_mesh_pi(a1, b1, kp, ki, sp)

    # saida real
    C_s = tf([kd, kp, ki], [1, 0])
    sys = tf(Num, Den)
    D = series(C_s, sys)
    E = feedback(D, 1)
    y_real, T_real = step(E*sp, t)
    # printa informacoes do sistema: overshoot, tempo de acomodacao e tempo de subida
    system_info_var.set(
        f"{system_info(E, t)[0]:.2f} (%), {system_info(E, t)[1]:.2f} (s)")

    # plotando graficos
    _, axs = plt.subplots(2, 2)
    # malha aberta
    axs[0, 0].plot(T, x1, 'b', label='x1')
    axs[0, 0].plot(T, y1, 'r', label='y1')
    axs[0, 0].set_ylabel('altura')
    axs[0, 0].legend()
    axs[0, 0].grid()
    axs[0, 0].set_title('Malha aberta')

    # malha fechada
    yout, T = step(mesh_close_tf*sp, t)
    axs[0, 1].plot(T, yout)
    axs[0, 1].set_title(f'Malha fechada com degrau de {sp}')
    axs[0, 1].set_ylabel('altura')
    axs[0, 1].grid()

    # saida do controlador
    axs[1, 0].plot(t, pv)
    axs[1, 0].set_title('Saida do controlador')
    axs[1, 0].set_ylabel('altura')
    axs[1, 0].grid()

    # saida real vs saida controlada
    axs[1, 1].plot(t, pv, 'b', label='controlada')
    axs[1, 1].plot(T_real, y_real, 'r', label='real')
    axs[1, 1].set_title('Saida real vs saida controlada')
    axs[1, 1].set_ylabel('altura')
    axs[1, 1].legend()
    axs[1, 1].grid()
    plt.show()


input_names = [
    "Set point", "Máximo pico (%)", "Tempo de acomodação (s)", "kp acrescimo", "ki acrescimo"]
for y, entry in enumerate(entries, start=1):
    entry.grid(row=y, column=1)
    Label(text=input_names[y-1]).grid(row=y, column=0)

Button(root, text="Calcular", command=on_press).grid(row=y+1, column=0)
Label(root, text="Informações do sistema:").grid(row=y+2, column=0)
Label(root, text="Overshot, Tempo de Acomodação").grid(row=y+3, column=0)
Label(root, textvariable=system_info_var).grid(row=y+3, column=1)
root.mainloop()
