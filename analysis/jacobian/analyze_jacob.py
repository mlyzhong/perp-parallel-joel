# Analyzes the jacobian matrix of
#  ./perp-parallel-joel oct all 1 >> 10sJacobian3.5e-6.txt
# format of datafile: csv:
#   t,x,y,z,L_1,L_2,L_3,L_4,L_5,L_6
# included: script to produce plot of trajectory
#
# MZ 20 June, 2017

import sys
#import plot_traj
import numpy as np
import matplotlib.pyplot as plt

plt.ioff()


filename = "10sJacobian3.5e-6.txt"
ignore_beginning = 3 # lines to ignore in the beginning
row_size = 10 # number of elements per row
dtim = 3.5E-6

def scatter(x,y, sample_size = 50000, s = 1, color = 'k'):
    "to save computation time, displays a random distribution of the scatterplot, instead of the whole damn thing"
    if len(x) < sample_size:
        plt.scatter(x, y, s = s, color = color)
        return
    else:
        indices = np.random.choice(len(x), size = sample_size, replace = False)
        plt.scatter(x[indices], y[indices], s = s, color = color)
    return

def load_file(filename):
    "opens file, returns array of formatted csv file; and other outputs"
    f = open(filename)

    # skips first couple of rows
    for _ in range(ignore_beginning):
        f.readline()

    data = []
    other = []
    line = "_"
    row = [0]
    while (line != ""): # and (row[0] < 1.0):
        line = f.readline()
        try:
            row = eval(line)
            if len(row) == row_size:
                data.append(row)
                if (row[0] % 1) < 0.01:
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
    scatter(x, y, s = 0.1)
    scatter(x[indices], y[indices], s = 0.1, color='r')
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

    x = data[:, 1]
    y = data[:, 2]
    z = data[:, 3]
    vx = np.append(0, np.diff(x) / dtim)
    vy = np.append(0, np.diff(y) / dtim)
    vz = np.append(0, np.diff(z) / dtim)
    ls = data[:, 4:]


    r = np.sqrt(x*x + y*y)
    vr = np.append(0, np.diff(r) / dtim)

# facts of the matter

# l1[0], l1[1] is corrupted
ls[0, :] = 0
ls[1, :] = 0

# see where sum is far from 1 -> buggy result
zeros = np.sum(ls, axis = 1)
# 3 outliers of close to zeros ... at [1029994], [1144035], and [1208838]
# from np.where(zeros == zeros.min())
# should I actually remove these though?
ls[1029994, :] = 0
ls[1144035, :] = 0
ls[1208838, :] = 0

l1 = data[:, 4]
l2 = data[:, 5]
l3 = data[:, 6]
l4 = data[:, 7]
l5 = data[:, 8]
l6 = data[:, 9]


# displays histogram of LLE
counts, bin_edges, _ = plt.hist(l1, bins = 10000, normed = True)
plt.xlabel("local LLE (1/s)")
plt.show()

# cdf it out!
cdf = (bin_edges[2] - bin_edges[1])*np.cumsum(counts)
plt.plot(bin_edges[1:], cdf)

# first, we will analyze where a very small LLE is, for the massive spike around l1 = 0

comparative_histograms(r, l1 < 0.5)
# seems as if small LLE means far away from r ...

# 3162060 points in total


indices = l1 < 0.5

comparative_scatter(r, z, l1 > 0.5)
