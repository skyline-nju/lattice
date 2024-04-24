import numpy as np
import matplotlib.pyplot as plt
import struct
import platform
import matplotlib.tri as mtri
import os
import glob
from matplotlib.colors import hsv_to_rgb


root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def map_v_to_rgb(theta, module, m_max=None):
    """
    Transform orientation and magnitude of velocity into rgb.

    Parameters:
    --------
    theta: array_like
        Orietation of velocity field.
    module: array_like
        Magnitude of velocity field.
    m_max: float, optional
        Max magnitude to show.

    Returns:
    --------
    RGB: array_like
        RGB corresponding to velocity fields.
    """
    H = theta / 360
    V = module
    if m_max is not None:
        V[V > m_max] = m_max
        V /= m_max
    S = np.ones_like(H)
    HSV = np.dstack((H, S, V))
    RGB = hsv_to_rgb(HSV)
    return RGB


def get_para(fin):
    basename = os.path.basename(fin)
    s = basename.rstrip(".bin").split("_")
    Lx = int(s[0].lstrip("L"))
    Ly = int(s[1])
    Dr = float(s[2].lstrip("Dr"))
    Dt = float(s[3].lstrip("Dt"))
    eta = float(s[4].lstrip("e"))
    rho0 = float(s[5].lstrip("r"))
    phi = float(s[6].lstrip("p"))
    v0 = float(s[7].lstrip("v"))
    seed = int(s[8].lstrip("s"))
    dt = int(s[9].lstrip("dt"))
    t0 = int(s[10].lstrip("t"))
    return Lx, Ly, Dr, Dt, eta, v0, rho0, phi, seed, dt, t0


def plot_tri_snap(tri, rho, mx, my, rho0, Lx, Ly, t, fout=None):
    px, py = np.zeros((2, Ly, Lx))
    mask = rho != 0
    px[mask] = mx[mask] / rho[mask]
    py[mask] = my[mask] / rho[mask]

    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5), constrained_layout=True, sharex=True, sharey=True)
    im1 = axes[0].tripcolor(tri, rho.flat, vmin=0, vmax=3 * rho0)
    im2 = axes[1].tripcolor(tri, px[:, :].flat, vmin=-1, vmax=1, cmap="bwr")
    im3 = axes[2].tripcolor(tri, py[:, :].flat, vmin=-1, vmax=1, cmap="bwr")
    im = [im1, im2, im3]

    titles = [r"$\rho(\mathbf{r})$", r"$p_x(\mathbf{r})$", r"$p_y(\mathbf{r})$"]
    extend = ["max", None, None]
    for i, ax in enumerate(axes):
        ax.axis("scaled")
        ax.set_xlim(0, Lx)
        ax.set_ylim(0, Ly * np.sqrt(3)/2)
        ax.set_title(titles[i], fontsize="x-large")
        fig.colorbar(im[i], ax=axes[i], orientation="horizontal", extend=extend[i])
    fig.suptitle(r"$t=%g$" % t, fontsize="x-large")

    if fout is not None:
        plt.savefig(fout)
    else:
        plt.show()
    plt.close()


def add_colorbar(ax, mmin, mmax, theta_min=0, theta_max=360, orientation="h"):
    """ Add colorbar for the RGB image plotted by plt.imshow() """
    V, H = np.mgrid[0:1:50j, 0:1:180j]
    if orientation == "v":
        V = V.T
        H = H.T
        box = [mmin, mmax, theta_min, theta_max]
    else:
        box = [theta_min, theta_max, mmin, mmax]
    S = np.ones_like(V)
    HSV = np.dstack((H, S, V))
    RGB = hsv_to_rgb(HSV)
    ax.imshow(RGB, origin='lower', extent=box, aspect='auto')
    theta_ticks = [0, 45, 90, 135, 180, 225, 270, 315, 360]

    if orientation == "h":
        ax.set_xticks(theta_ticks)
        ax.set_xticklabels([r"$%d\degree$" % i for i in theta_ticks])
        ax.set_ylabel(r'$|\mathbf{m}|/\rho_0$', fontsize="large")
        ax.set_xlabel(r"$\theta$", fontsize="large")
    else:
        ax.yaxis.set_label_position('right')
        ax.yaxis.set_ticks_position("right")
        ax.set_yticks(theta_ticks)
        ax.set_yticklabels([r"$%d\degree$" % i for i in theta_ticks])
        ax.set_ylabel(r'orientation $\theta$', fontsize="large")
        ax.set_title(r"$|\mathbf{m}|/\rho_0$", fontsize="large")


def plot_rect_snap(tri, rho, mx, my, xi, yi, rho0, Lx, Ly, t, fout=None):
    interp_lin_rho = mtri.LinearTriInterpolator(tri, rho.flat)
    interp_lin_mx = mtri.LinearTriInterpolator(tri, mx.flat)
    interp_lin_my = mtri.LinearTriInterpolator(tri, my.flat)
    rho_lin = interp_lin_rho(xi, yi)
    mx_lin = interp_lin_mx(xi, yi)
    my_lin = interp_lin_my(xi, yi)

    theta = np.arctan2(my_lin, mx_lin)
    theta[theta < 0] += np.pi * 2
    theta *= 180 / np.pi
    module = np.sqrt(mx_lin ** 2 + my_lin ** 2)
    RGB = map_v_to_rgb(theta, module, m_max=rho0)
    box = [0, Lx, 0, Ly / 2 * np.sqrt(3)]

    fig, axes = plt.subplots(1, 2, figsize=(8, 4.5), sharex=True, sharey=True)
    im1 = axes[0].imshow(rho_lin, extent=box, origin="lower", vmin=0, vmax=3 * rho0)
    axes[1].imshow(RGB, extent=box, origin="lower")

    plt.tight_layout(rect=[-0.01, -0.01, 1.01, 0.99])
    bbox1 = axes[0].get_position().get_points().flatten()
    bbox2 = axes[1].get_position().get_points().flatten()
    fig.subplots_adjust(bottom=0.25)
    bbox1[1], bbox1[3] = 0.1, 0.08
    bbox1[0] += 0.025
    bbox1[2] = bbox1[2] - bbox1[0] - 0.025
    bbox2[1], bbox2[3] = 0.1, 0.1
    bbox2[0] += 0.025
    bbox2[2] = bbox2[2] - bbox2[0] - 0.025
    cb_ax1 = fig.add_axes(bbox1)
    cb_ax2 = fig.add_axes(bbox2)
    cb1 = fig.colorbar(im1, ax=axes[0], cax=cb_ax1, orientation="horizontal", extend="max")
    add_colorbar(cb_ax2, 0, 1, 0, 360, "h")
    cb1.set_label(r"$\rho$", fontsize="large")

    fig.suptitle(r"$t=%g$" % t, fontsize="x-large", y=0.998)

    if fout is not None:
        plt.savefig(fout)
    else:
        plt.show()
    plt.close()


def show_q6(fin=None, savefig=False):
    if fin is None:
        Lx = 256
        Ly = 256
        
        rho0 = 10
        phi = 10

        Dt = 0.07
        Dr = 0.05
        v0 = 1.76
        eta = -2
        
        dt = 5000
        t_beg = 0
        seed = 1000
        if (platform.system() == "Linux"):
            folder = f"{root_sohrab}/lat_QS/{Lx:d}_{Ly:d}_q6"
        elif (platform.system() == "Windows"):
            folder = "data"
        fin = f"{folder}/L{Lx:d}_{Ly:d}_Dr{Dr:.3f}_Dt{Dt:.3f}_e{eta:g}_r{rho0:g}_p{phi:g}_v{v0:g}_s{seed:d}_dt{dt:d}_t{t_beg:d}.bin"
    else:
        Lx, Ly, Dr, Dt, eta, v0, rho0, phi, seed, dt, t_beg = get_para(fin)
    print(fin)
    frame_size = Lx * Ly * 6
    with open(fin, "rb") as f:
        f.seek(0, 2)
        filesize = f.tell()

        n_frames = filesize//frame_size
        print("find", n_frames, "frames")
        beg = 0
        if savefig:
            outfolder = f"fig/Dr{Dr:.3f}_Dt{Dt:.3f}_e{eta:g}_r{rho0:g}_p{phi:g}_v{v0:g}_s{seed}"
            if not os.path.exists(outfolder):
                os.mkdir(outfolder)
        else:
            if n_frames > 10:
                beg = n_frames - 10
        f.seek(frame_size * beg)

        x = np.arange(Lx) + 0.5
        y = (np.arange(Ly) + 0.5) * np.sqrt(3) / 2
        xx, yy = np.meshgrid(x, y)
        for i in range(y.size):
            xx[i] -= i * 0.5
            mask = xx[i] < 0
            xx[i][mask] += Lx
        tri = mtri.Triangulation(xx.flat, yy.flat)

        xi, yi = np.meshgrid(
            np.linspace(0.5, Lx, Lx, endpoint=False),
            np.linspace(0.5 * np.sqrt(3)/2, Ly * np.sqrt(3)/2, int(Ly * np.sqrt(3)/2), endpoint=False)
        )

        while f.tell() < filesize:
            i_frame = f.tell() // frame_size
            if savefig:
                fout = f"{outfolder}/{i_frame:05d}.jpg"
                if os.path.exists(fout):
                    f.seek(frame_size, 1)
                    continue
            else:
                fout = None
            buf = f.read(frame_size)
            data = np.array(struct.unpack("%dB" % (Lx * Ly * 6), buf)).reshape(Ly, Lx, 6)
            rho = np.sum(data, axis=2)
            mx = data[:, : , 0] + 0.5 * data[:, :, 1] - 0.5 * data[:, :, 2] - data[:, :, 3] - 0.5 * data[:, :, 4] + 0.5 * data[:, :, 5]
            my = np.sqrt(3) / 2 * (data[:, :, 1] + data[:, :, 2] - data[:, :, 4] - data[:, :, 5])
            print("frame", f.tell()//frame_size, ", max particle number:", data.max())
 
            t = (i_frame + 1) * dt + t_beg
            # plot_tri_snap(tri, rho, mx, my, rho0, Lx, Ly, t, fout)
            plot_rect_snap(tri, rho, mx, my, xi, yi, rho0, Lx, Ly, t, fout)

def update_figs():
    folder = f"{root_sohrab}/lat_QS/256_256_q6"
    files = glob.glob(f"{folder}/*")
    for fin in files:
        show_q6(fin, only_rho=False, savefig=True)


if __name__ == "__main__":
    show_q6(savefig=False)
    # update_figs()