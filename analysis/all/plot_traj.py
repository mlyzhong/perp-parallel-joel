# plot trajectory function
# produces a figure that looks
# MZ 8 June, 2017


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec



style = 'coral'
figsize = (5, 6)


plt.ion()


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=14)


def plot_traj(ts, xs, ys, zs, tmin, tmax):
    indices = (ts >= tmin) * (ts <= tmax)
    t = ts[indices]
    x = 1000*xs[indices] # from m to mm
    y = 1000*ys[indices] # from m to mm
    z = 1000*zs[indices] # from m to mm

    plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(5, 4)

    plt.subplot(gs[:3, :4]) # t vs. z
    plt.plot(t, z, style)

    plt.xlabel(r'$\mathrm{Time\ (ms)}$')
    plt.ylabel(r"$z\mathrm{\ (mm)}$")
    plt.xlim([t.min(), t.max()])
    plt.text(2, 40, '(a)')

    # z vs. r
    plt.subplot(gs[3:, :2])
    plt.plot(z, np.sqrt(x*x + y*y), style)
    plt.xlabel(r"$z\mathrm{\ (mm)}$")
    plt.ylabel(r"$r\mathrm{\ (mm)}$")
    plt.text(-40, 3, '(b)')
    plt.ylim([0, 22])
    plt.yticks([0, 5, 10, 15, 20])
    plt.xticks([-40, 0, 40])

    # x vs. y
    ax = plt.subplot(gs[3:, 2:4])
    ax.yaxis.tick_right()
    ax.yaxis.set_ticks_position('both')
    ax.yaxis.set_label_position("right")
    plt.plot(x, y, style)
    plt.xlabel(r"$x\mathrm{\ (mm)}$")
    plt.ylabel(r"$y\mathrm{\ (mm)}$")
    plt.text(-11, -17, '(c)')
    plt.xlim([-22, 22])
    plt.ylim([-22, 22])
    plt.xticks([-20, 0, 20])
    plt.yticks([-20, -10, 0, 10, 20])

    plt.savefig("fig.pdf")
    return
