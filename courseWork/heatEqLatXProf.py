import numpy as np
import matplotlib.pyplot as plt


# Задача
a = 110
length = np.pi #mm
time = 16 #seconds
nodes = 101

# Инициализация
dt = 1/5000
dx = np.sqrt(2 * a * dt)
nodes = 101

t_nodes = int(time/dt)


# Начальные условия
d = 2
lam = 2.71828**3
psi = 10
xx = np.linspace(-1*np.pi /2.0, np.pi /2.0, nodes+1)

p = 25
# u = np.exp(np.abs(xx))*np.sin(np.abs(xx))/np.exp(np.pi /2.0)*1000
u = np.zeros(xx.shape)+24

def F(x):
    if x >= p:
        return 0
    return x

# Визуализация с графикой
fig, axis = plt.subplots(1,2)

# pcm = axis.pcolormesh([u], cmap='rainbow', vmin=0, vmax=2000)
# plt.colorbar(pcm, ax=axis)
# axis.set_ylim([-2, 3])


axis[0].plot(xx + np.pi,u)
pcm = axis[1].pcolormesh([u], cmap='rainbow', vmin=0, vmax=200)
plt.colorbar(pcm, ax=axis[1])
axis[1].set_ylim([-2, 3])


# Симуляция
t = 1

w = []
wX = []

Tlat = 2
latencyStep = int(np.round(Tlat/dt))

T = latencyStep

M = 0

while t < t_nodes  :
    a = np.copy(u)
    w.append(a)

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
                k2 = -1 * (w[t-1][i] + k1 * dt / 2.0) + d * (2 * w[t-1][i - 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i] + k1 * dt / 2.0) * M
                k3 = -1 * (w[t-1][i] + k2 * dt / 2.0) + d * (2 * w[t-1][i - 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i] + k2 * dt/ 2.0) * M
                k4 = -1 * (w[t-1][i] + k3 * dt) + d * (2 * w[t-1][i - 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i] + k3 * dt) * M
                u[i] = w[t-1][i] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        else:
            k1 = -1 * w[t-1][i] + d * (w[t-1][i - 1] + w[t-1][i + 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i]) * M
            k2 = -1 * (w[t-1][i] + k1 * dt / 2.0) + d * (w[t-1][i - 1] + w[t-1][i + 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i]+ k1 * dt / 2.0) * M
            k3 = -1 * (w[t-1][i] + k2 * dt / 2.0) + d * (w[t-1][i - 1] + w[t-1][i + 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i]+ k2 * dt / 2.0) * M
            k4 = -1 * (w[t-1][i] + k3 * dt) + d * (w[t-1][i - 1] + w[t-1][i + 1] - 2 * w[t-1][i]) / dx ** 2 + lam*F(w[t-T-1][i] + k3 * dt) * M
            u[i] = np.copy(w[t-1][i] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0)

    wX.append(w[t-1][10])
    t += 1


    print("t: {:.3f} [s], Average temperature: {:.2f} Celcius".format(dt*t, np.average(u)))

    # Обновление графика

    if t % 500 == 0:
        axis[0].cla()
        axis[1].cla()

        axis[0].scatter(np.linspace(0,dt*t,t-1), wX)

        pcm = axis[1].pcolormesh([u], cmap='rainbow', vmin=0, vmax=100)
        axis[1].set_ylim([-2, 3])

        axis[0].set_ylim([0, 100])

        if t > 1:
            plt.title("Distribution at t: {:.3f} [s]. M: {:}, F: {:} \n Average temperature: {:.2f} Celcius \n Average temperature latency: {:.2f} Celcius".format(dt*t,M, M*F(w[t-T-2][i]), np.average(u), np.average(w[t-T-2])))
        plt.pause(0.0001)


plt.show()
