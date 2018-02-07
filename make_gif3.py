#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# author:Hirotaka Kondo

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import csv
import math

mx = 401
my = 201
r = 60  # 全体領域の円の半径
del_r = 1  # 円柱の半径
log_r = np.log(r)

x, y = [], []
for i in range(my):
    x_row = []
    for j in range(mx):
        x_row.append(np.exp(j / mx * np.log(r)) * np.cos(math.pi - i / my * 2 * math.pi))
    x.append(x_row)
    # x.append(
    #   np.linspace(del_r * np.cos(math.pi - i / my * 2 * math.pi), r * np.cos(math.pi - i / my * 2 * math.pi), mx))

for i in range(my):
    y_row = []
    for j in range(mx):
        y_row.append(np.exp(j / mx * np.log(r)) * np.sin(math.pi - i / my * 2 * math.pi))
    y.append(y_row)
    # y.append(
    #    np.linspace(del_r * np.sin(math.pi - i / my * 2 * math.pi), r * np.sin(math.pi - i / my * 2 * math.pi), mx))

p_np = []


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
    for i in range(0, 50):
        print(i)
        p_np_part = read_p_csv("../data50/" + str(10000 + i * 200) + "data.csv")
        p_np.append(p_np_part)
        im = plt.imshow(p_np[-1], extent=[np.min(x), np.max(x),
                                          np.min(y), np.max(y)], vmin=-0.1, vmax=0.1)  # vmin,vmaxで等高線のカラーリングを調整
        ims.append([im])

    anim = animation.ArtistAnimation(fig, ims, interval=5, blit=False)
    anim.save('karman.gif', writer="imagemagick")
    plt.show()


def main2():
    """
    t=7000での等圧線をplot
    :return:
    """
    p = read_p_csv("../data80/100000data.csv")
    plt.contour(x, y, p, 500)
    #plt.imshow(p, extent=[np.min(x), np.max(x), np.min(y), np.max(y)], vmin=-0.1,
    #        vmax=0.1)  # vmin,vmaxで等高線のカラーリングを調整
    plt.gca().set_aspect('equal')
    plt.colorbar()
    plt.show()

if __name__ == "__main__":
    main2()
