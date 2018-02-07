#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author:Hirotaka Kondo

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import csv

x_min = -10.0
x_max = 30.0
y_min = -10.0
y_max = 10.0
n_last = 5000
x_np = np.linspace(x_min, x_max, 401)
y_np = np.linspace(y_max, y_min, 201)
p_np = []
x, y = np.meshgrid(x_np, y_np)


def read_p_csv(file_name):
    p_lis = []
    with open(file_name, 'r') as f:
        reader = csv.reader(f)
        for row in reader:  # y座標を固定したときの各xの圧力pのリストがrow
            p_row = []
            for p in row:
                p_row.append(round(float(p), 3))
            p_lis.append(p_row)
    p_np_part = np.array(p_lis)
    return p_np_part


def main1():
    """
    x,yのgrid格子をベースに,時間tが変化していったときの圧力分布を
    gifとして出力する.
    :return:
    """
    fig = plt.figure()
    ims = []
    for i in range(0, 200):
        print(i)
        p_np_part = read_p_csv("../data_pressure_kakutyu/" + str(2000 + i * 40) + "data.csv")
        p_np.append(p_np_part)
        im = plt.imshow(p_np[-1], extent=[np.min(x), np.max(x),
                                          np.min(y), np.max(y)], vmin=-0.1, vmax=0.1)  # vmin,vmaxで等高線のカラーリングを調整
        ims.append([im])

    anim = animation.ArtistAnimation(fig, ims, interval=5, blit=False)
    anim.save('karman1.gif', writer="imagemagick")
    plt.show()


def main2():
    """
    t=7000での等圧線をplot
    :return:
    """
    p = read_p_csv("../data/7000data.csv")
    plt.imshow(p, extent=[np.min(x), np.max(x), np.min(y), np.max(y)], vmin=-0.1,
               vmax=0.1)  # vmin,vmaxで等高線のカラーリングを調整
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.show()


def main3():
    cd, cl, cp1, cp2 = [], [], [], []
    t = np.linspace(1, 10000, 10000)
    with open("cdclcp.csv", 'r') as f:
        reader = csv.reader(f)
        for row in reader:  # y座標を固定したときの各xの圧力pのリストがrow
            cd.append(float(row[0]))
            cl.append(float(row[1]))
            cp1.append(float(row[2]))
            cp2.append(float(row[3]))
    plt.plot(t, cd, label="Cd")
    plt.plot(t, cl, label="Cl")
    plt.plot(t, cp1, label="Cp1")
    plt.plot(t, cp2, label="Cp2")
    plt.xlabel("step")
    plt.xlim([0, 5000])
    plt.ylim([-1, 3])
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main1()
