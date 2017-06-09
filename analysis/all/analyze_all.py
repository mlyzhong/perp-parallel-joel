# Analyzes a data file produced by perp-parallel-joel, test_type: all
# format of datafile: csv:
#   t,x,y,z,perp_x,perp_y,perp_z,perp_vx,perp_vy,perp_vz,dlog|perp|/dt
# included: script to produce plot of trajectory
#
# MZ 8 June, 2017

import sys
#import plot_traj
import numpy as np
import matplotlib.pyplot as plt

plt.ioff()


filename = "1sec3.5e-5.txt"
ignore_beginning = 3 # lines to ignore in the beginning
row_size = 11 # number of elements per row

def load_file(filename):
    "opens file, returns array of formatted csv file; and other outputs"
    f = open(filename)

    # skips first couple of rows
    for _ in range(ignore_beginning):
        f.readline()

    data = []
    other = []
    line = "_"
    while line != "":
        line = f.readline()
        try:
            row = eval(line)
            if len(row) == 11:
                data.append(row)
                if (row[0] % 1) < 0.001:
                    print(row[0])
            else:
                other.append(row)
        except:
            continue

    print("done")
    data = np.array(data)
    return data, other


def comparative_histograms(dist, indices, bins = 100, xlabel = ""):
    "plots normalized histogram of dist (green) and dist[indices] (red) overlayed over each other"
    n, bins, patches = plt.hist(dist, bins = bins, normed = True)
    plt.setp(patches, 'facecolor', 'g', 'alpha', 0.5)
    n, bins, patches2 = plt.hist(dist[indices], bins = bins, normed = True)
    plt.setp(patches2, 'facecolor', 'r', 'alpha', 0.5)
    plt.xlabel(xlabel)
    plt.show()

def comparative_scatter(x, y, indices, xlabel = "", ylabel = ""):
    "plots two scatter plots: x vs. y, and x[indices] vs. y[indices] (red) "
    plt.scatter(x, y, s = 0.1)
    plt.scatter(x[indices], y[indices], s = 0.8, color='r')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim([x.min(), x.max()])
    plt.ylim([y.min(), y.max()])
    plt.show()


def comparative_scatter_3(x, y, indices1, indices2, xlabel = "", ylabel = ""):
    """
    plots scatter plots:
    x vs. y (black)
    x[indices1] ... (red)
    x[indices2] ... (blue)
    x[indices1 and indices2] ...
    """
    indices_and = indices1 * indices2
    plt.scatter(x, y, s = 0.8)
    plt.scatter(x[indices1], y[indices1], s = 0.8, color='r')
    plt.scatter(x[indices2], y[indices2], s = 0.8, color='b')
    plt.scatter(x[indices_and], y[indices_and], s = 0.8, color='m')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xlim([x.min(), x.max()])
    plt.ylim([y.min(), y.max()])
    plt.show()


def analyze_lyapunov(data):
    """
    Analyzes instantaneous lyapunov
    We'll consider where, specifically, d(log|pert|)/dt > 2000 s^{-1}
    """
    dlnpertdt = data[:, -1]
    x = data[:, 1]
    y = data[:, 2]
    z = data[:, 3]
    dx = data[:, 4]
    dy = data[:, 5]
    dz = data[:, 6]
    dvx = data[:, 7]
    dvy = data[:, 8]
    dvz = data[:, 9]
    r = np.sqrt(x*x + y*y)
    dr = np.sqrt(dx*dx + dy*dy)
    dv = np.sqrt(dvx*dvx + dvy*dvy + dvz*dvz)

    """
    # study 1
    # shows that most growth occurs at small r**2
    # plt.hist(data[:,-1], bins = 100)
    indices = data[:, -1] > 2000
    plt.hist(data[:, -1], bins=100) # instantaneous lyapunov histogram
    plt.show()
    comparative_histograms(x, indices)
    comparative_histograms(y, indices)
    comparative_histograms(z, indices)
    comparative_histograms(r, indices)
    """


    # shows that r, dlogpert/dt are correlated
    plt.scatter(dlnpertdt, r, s = 0.1)
    plt.xlabel("dlog|pert|/dt (sec^{-1})")
    plt.ylabel("r (m)")



indices_r = r < 0.002
indices_pert = data[:, -1] > 2000

    comparative_histograms(r, indices_pert, xlabel = "r (m)")
    comparative_histograms(dlnpertdt, indices_r, xlabel = "d log|pert| / dt (s^{-1})")




    # study 2
    # examine points in trajectory in which r is close to center axis
    #plt.hist(r, bins=100)
    #plt.show()

    plt.show()
    plt.cla()

    comparative_scatter_3(dr, dv, indices_r, indices_pert, "dr", "dv")
    comparative_histograms(data[:, -1], indices)
    comparative_scatter(dx, dy, indices, "dx", "dy")
    comparative_scatter(dx, dz, indices, "dx", "dz")
    comparative_scatter(dvx, dvy, indices, "dvx", "dvy")
    comparative_scatter(dvx, dvz, indices, "dvx", "dvz")
    comparative_scatter(dx, dvx, indices, "dx", "dvx")
    comparative_scatter(dx, dvy, indices, "dx", "dvy")
    comparative_scatter(dy, dvy, indices, "dy", "dvy")
    comparative_scatter(dz, dvz, indices, "dz", "dvz")


    """
    comparative_histograms(data[:, -1], indices)
    comparative_histograms(dx, indices)
    comparative_histograms(dy, indices)
    comparative_histograms(dz, indices)
    """



if __name__ == "__main__":
    try:
        filename = sys.argv[1]
    except:
        print("You need to specify the filename of what's to be analyzed!")
    data, other = load_file(filename)
    #analyze_lyapunov(data)

    dlnpertdt = data[:, -1]
    x = data[:, 1]
    y = data[:, 2]
    z = data[:, 3]
    vz = np.append(0, np.diff(z) / 3.5E-5)
    dx = data[:, 4]
    dy = data[:, 5]
    dz = data[:, 6]
    dvx = data[:, 7]
    dvy = data[:, 8]
    dvz = data[:, 9]
    r = np.sqrt(x*x + y*y)
    dr = np.sqrt(dx*dx + dy*dy)
    dv = np.sqrt(dvx*dvx + dvy*dvy + dvz*dvz)
