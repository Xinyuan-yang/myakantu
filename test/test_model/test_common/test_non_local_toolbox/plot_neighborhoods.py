#!/usr/bin/env python
__copyright__ = (
    "Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)"
    "Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)"
)
__license__ = "LGPLv3"


import sys
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.backends.backend_pdf import PdfPages


# ------------------------------------------------------------------------------
def read_file(filename):
    """Read a reference file."""
    quad_x_coords = list()
    quad_y_coords = list()
    with open(filename) as fh:
        for line in fh:
            if line.strip().startswith("#"):
                line = fh.next()
                x, y = [float(s) for s in line[1:-3].split(", ")]
                quad_x_coords.append(float(x))
                quad_y_coords.append(float(y))
    return quad_x_coords, quad_y_coords


# ------------------------------------------------------------------------------
def read_neighborhoods(filename):
    """Read a neighborhood file or verification."""
    neighborhoods = dict()
    with open(filename) as fh:
        first_line = fh.readline()
        key = int(first_line[0:-1].split(" ")[-1])
        current_list = list()

        for line in fh:
            if line.strip().startswith("#"):
                neighborhoods[key] = current_list
                key = int(line[0:-1].split(" ")[-1])
                current_list = list()
            else:
                coord = [float(s) for s in line[1:-3].split(", ")]
                current_list.append(coord)
        neighborhoods[key] = current_list
    return neighborhoods


# ------------------------------------------------------------------------------
def plot_neighborhood(quad_x_coords, quad_y_coords, neighborhood, r):
    """Plot a neighborhood."""
    ax = plt.subplot(111)
    plt.axis([-1, 1, -1, 1])
    # plot all quads
    ax.plot(
        quad_x_coords,
        quad_y_coords,
        label="free",
        linestyle="None",
        marker="+",
        markeredgewidth=0.2,
        markersize=4,
    )

    # plot center of neighborhood and circle
    x = neighborhood[0][0]
    y = neighborhood[0][1]
    ax.plot(
        x,
        y,
        label="free",
        marker="+",
        markeredgewidth=0.2,
        markersize=4,
        color="r",
    )
    circ = plt.Circle((x, y), radius=r, color="g", fill=False)
    ax.add_patch(circ)
    formatter = ticker.ScalarFormatter()
    formatter.set_powerlimits((-2, 2))
    ax.xaxis.set_major_formatter(formatter)

    # plot all the neighbors
    for i in range(1, len(neighborhood)):
        ax.plot(
            neighborhood[i][0],
            neighborhood[i][1],
            label="free",
            linestyle="None",
            marker="x",
            markeredgewidth=0.2,
            markersize=4,
            color="y",
        )

    ax.grid(b=True, which="major", color="0.75", linestyle="-", linewidth=0.5)
    ax.grid(b=True, which="minor", color="0.75", linestyle="-", linewidth=0.5)


# ------------------------------------------------------------------------------
def main():
    """Run the main function."""
    non_local_radius = 0.5
    rc("text", usetex=True)
    rc("font", family="serif", size=8, serif="Times")

    neighborhood_file = sys.argv[1]

    quad_x_coords, quad_y_coords = read_file(neighborhood_file)
    nb_quads = len(quad_x_coords)

    neighborhoods = read_neighborhoods(neighborhood_file)

    with PdfPages("resulting_neighborhoods.pdf") as pdf:
        for i in range(nb_quads):
            plt.figure(1, figsize=(7 / 2.54, 7 / 2.54))
            plot_neighborhood(
                quad_x_coords,
                quad_y_coords,
                neighborhoods[i],
                non_local_radius,
            )
            pdf.savefig()
            plt.close()


if __name__ == "__main__":
    main()
