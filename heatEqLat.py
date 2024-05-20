import numpy as np
import matplotlib.pyplot as plt


# Задача
a = 110
length = np.pi #mm
time = 10 #seconds
nodes = 101

# Инициализация
dt = 1/10000
dx = np.sqrt(2 * a * dt)
nodes = int(np.round(length/dx))

t_nodes = int(time/dt)


# Начальные условия
d = 10
lam = 200
psi = 10
xx = np.linspace(-1*np.pi /2.0, np.pi /2.0, nodes+1)

p = 400
# u = np.exp(np.abs(xx))*np.sin(np.abs(xx))/np.exp(np.pi /2.0)*1000
u = np.sin(np.abs(xx))*lam*psi

def F(u):
    if u > p:
        return 0
    return 500*u

# Визуализация с графикой
fig, axis = plt.subplots()

pcm = axis.pcolormesh([u], cmap='rainbow', vmin=0, vmax=2000)
plt.colorbar(pcm, ax=axis)
axis.set_ylim([-2, 3])

# Симуляция
t = 1

w = []

Tlat = 1
latencyStep = int(np.round(Tlat/dt))

T = latencyStep

M = 0

while t < t_nodes  :

    w.append(u)

    M = 1
    T = latencyStep
    if (t-latencyStep) < 0:
        M = 0
        T = 0

    for i in range(0, nodes+1):

        if i == 0 or i == nodes:
            if i == 0:
                k1 = -1 * w[t-1][i] + d * (2 * w[t-1][i + 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i]) * M
                k2 = -1 * (w[t-1][i] + k1 * dt / 2.0) + d * (2 * w[t-1][i + 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i] + k1 * dt / 2.0) * M
                k3 = -1 * (w[t-1][i] + k2 * dt / 2.0) + d * (2 * w[t-1][i + 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i] + k2 * dt / 2.0) * M
                k4 = -1 * (w[t-1][i] + k3 * dt) + d * (2 * w[t-1][i + 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i] + k3 * dt) * M
                u[i] = w[t-1][i] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
            else:
                k1 = -1 * w[t-1][i] + d * (2 * w[t-1][i - 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i]) * M
                k2 = -1 * (w[t-1][i] + k1 * dt / 2.0) + d * (2 * w[t-1][i - 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i] + k1 * dt ) * M
                k3 = -1 * (w[t-1][i] + k2 * dt / 2.0) + d * (2 * w[t-1][i - 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i] + k2 * dt) * M
                k4 = -1 * (w[t-1][i] + k3 * dt) + d * (2 * w[t-1][i - 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i] + k3 * dt) * M
                u[i] = w[t-1][i] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        else:
            k1 = -1 * w[t-1][i] + d * (w[t-1][i - 1] + w[t-1][i + 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i]) * M
            k2 = -1 * (w[t-1][i] + k1 * dt / 2.0) + d * (w[t-1][i - 1] + w[t-1][i + 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i]+ k1 * dt / 2.0) * M
            k3 = -1 * (w[t-1][i] + k2 * dt / 2.0) + d * (w[t-1][i - 1] + w[t-1][i + 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i]+ k2 * dt / 2.0) * M
            k4 = -1 * (w[t-1][i] + k3 * dt) + d * (w[t-1][i - 1] + w[t-1][i + 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i] + k3 * dt) * M
            u[i] = w[t-1][i] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0

    t += 1

    print("t: {:.3f} [s], Average temperature: {:.2f} Celcius".format(dt*t, np.average(u)))

    # Обновление графика

    if t % 50 == 0:
        pcm.set_array([u])
        if t > 1:
            axis.set_title("Distribution at t: {:.3f} [s]. M: {:}, F: {:}".format(dt*t,M, M*F(w[t-T-2][i])))
        plt.pause(0.0001)


plt.show()











