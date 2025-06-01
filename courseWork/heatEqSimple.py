import numpy as np
import matplotlib.pyplot as plt

# Задача
a = 10
length = np.pi #mm
time = 5 #seconds
nodes = 101

# Инициализация
dt = 1/10000
dx = np.sqrt(2 * a * dt)
nodes = 101

t_nodes = int(time/dt)

# u = np.zeros(nodes+1) + 50# Plate is initially as 20 degres C
xx = np.linspace(-1*np.pi /2.0, np.pi /2.0, nodes+1)
# Начальные условия
d = 5
lam = 500
psi = 2
u = np.sin(np.abs(xx))*100


# Визуализация с графикой
fig, axis = plt.subplots()

plt.plot(xx + np.pi,u)

# Симуляция
counter = 0

while counter < time  :

    w = u.copy()

    for i in range(0, nodes+1):

        if i == 0 or i == nodes:
            if i == 0:
                k1 = -1 * w[i] + d * (2 * w[i + 1] - 2 * w[i]) / dx ** 2
                k2 = -1 * (w[i] + k1 * dt / 2.0) + d * (2 * w[i + 1] - 2 * w[i]) / dx ** 2
                k3 = -1 * (w[i] + k2 * dt / 2.0) + d * (2 * w[i + 1] - 2 * w[i]) / dx ** 2
                k4 = -1 * (w[i] + k3 * dt) + d * (2 * w[i + 1] - 2 * w[i]) / dx ** 2
                u[i] = w[i] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
            else:
                k1 = -1 * w[i] + d * (2 * w[i - 1] - 2 * w[i]) / dx ** 2
                k2 = -1 * (w[i] + k1 * dt / 2.0) + d * (2 * w[i - 1] - 2 * w[i]) / dx ** 2
                k3 = -1 * (w[i] + k2 * dt / 2.0) + d * (2 * w[i - 1] - 2 * w[i]) / dx ** 2
                k4 = -1 * (w[i] + k3 * dt) + d * (2 * w[i - 1] - 2 * w[i]) / dx ** 2
                u[i] = w[i] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        else:
            k1 = -1 * w[i] + d * (w[i - 1] + w[i + 1] - 2 * w[i]) / dx ** 2
            k2 = -1 * (w[i] + k1 * dt / 2.0) + d * (w[i - 1] + w[i + 1] - 2 * w[i]) / dx ** 2
            k3 = -1 * (w[i] + k2 * dt / 2.0) + d * (w[i - 1] + w[i + 1] - 2 * w[i]) / dx ** 2
            k4 = -1 * (w[i] + k3 * dt) + d * (w[i - 1] + w[i + 1] - 2 * w[i]) / dx ** 2
            u[i] = w[i] + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0

    counter += dt

    print("t: {:.3f} [s], Average temperature: {:.2f} Celcius".format(counter, np.average(u)))

    # Обновление графика
    if int(counter/dt) % 10 == 0:
        plt.clf()
        plt.plot(xx + np.pi/2, u)
        plt.ylim(0, 100)
        plt.xlim(0, np.pi)

        axis.set_title("Distribution at t: {:.3f} [s].".format(counter))
        plt.pause(0.0001)


plt.show()