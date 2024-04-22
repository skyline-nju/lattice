import numpy as np
import matplotlib.pyplot as plt
import struct
import platform
import matplotlib.tri as mtri
import os
import glob


root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def get_para(fin):
    basename = os.path.basename(fin)
    s = basename.rstrip(".bin").split("_")
    Lx = int(s[0].lstrip("L"))
    Ly = int(s[1])
    Dr = float(s[2].lstrip("Dr"))
    Dt = float(s[3].lstrip("Dt"))
    etaAA = float(s[4].lstrip("e"))
    etaBB = float(s[5])
    etaAB = float(s[6].lstrip("J"))
    etaBA = float(s[7])
    v0 = float(s[8].lstrip("v"))
    rho0 = float(s[9].lstrip("r"))
    phiA = float(s[10].lstrip("p"))
    phiB = float(s[11])
    seed = int(s[12].lstrip("s"))
    dt = int(s[13].lstrip("dt"))
    t0 = int(s[14].lstrip("t"))
    return Lx, Ly, Dr, Dt, etaAA, etaBB, etaAB, etaBA, v0, rho0, phiA, phiB, seed, dt, t0

def show_q6(fin=None, only_rho=False, savefig=False):
    if fin is None:
        Lx = 256
        Ly = 256
        
        rho0 = 10
        phiA = 0
        phiB = 10

        Dt = 0.1
        Dr = 0.15
        v0 = 1
        etaAA = -1
        etaBB = -2
        etaAB = 1
        etaBA = -etaAB
        
        dt = 10000
        t_beg = 0
        seed = 1000
        if (platform.system() == "Linux"):
            if phiA != 0:
                folder = f"{root_sohrab}/lat_NRQS/{Lx:d}_{Ly:d}_q6"
            else:
                folder = f"{root_sohrab}/lat_NRQS/{Lx:d}_{Ly:d}_pA0_q6"
        elif (platform.system() == "Windows"):
            folder = "data"
        if phiA == rho0 and phiB == rho0:
            fin = "%s/L%d_%d_Dr%g_Dt%g_e%g_%g_J%g_%g_v%g_r%g_s%d_dt%d_t%d.bin" % (
                folder, Lx, Ly, Dr, Dt, etaAA, etaBB, etaAB, etaBA, v0, rho0, seed, dt, t_beg)
        else:
            fin = "%s/L%d_%d_Dr%g_Dt%g_e%g_%g_J%g_%g_v%g_r%g_p%g_%g_s%d_dt%d_t%d.bin" % (
                folder, Lx, Ly, Dr, Dt, etaAA, etaBB, etaAB, etaBA, v0, rho0, phiA, phiB, seed, dt, t_beg)
    else:
        Lx, Ly, Dr, Dt, etaAA, etaBB, etaAB, etaBA, v0, rho0, phiA, phiB, seed, dt, t_beg = get_para(fin)
    print(fin)
    frame_size = Lx * Ly * 12 * 2
    with open(fin, "rb") as f:
        f.seek(0, 2)
        filesize = f.tell()

        n_frames = filesize//frame_size
        print("find", n_frames, "frames")
        beg = 0
        if savefig:
            if phiA == 0:
                outfolder = f"fig/Dr{Dr:.3f}_Dt{Dt:.3f}_e{etaBB:g}_r{rho0:g}_p{phiB:g}_s{seed}"
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

        while f.tell() < filesize:
            i_frame = f.tell() // frame_size
            if savefig:
                fout = f"{outfolder}/{i_frame:05d}.jpg"
                if os.path.exists(fout):
                    f.seek(frame_size, 1)
                    continue
            buf = f.read(frame_size)
            data = np.array(struct.unpack("%dH" % (Lx * Ly * 12), buf)).reshape(Ly, Lx, 2, 6)
            rho = np.sum(data, axis=3)
            mx = data[:, :, :, 0] + 0.5 * data[:, :, :, 1] - 0.5 * data[:, :, :, 2] - data[:, :, :, 3] - 0.5 * data[:, :, :, 4] + 0.5 * data[:, :, :, 5]
            my = np.sqrt(3) / 2 * (data[:, :, :, 1] + data[:, :, :, 2] - data[:, :, :, 4] - data[:, :, :, 5])
            print("frame", f.tell()//frame_size, "max particle number:", data.max())
 
            if not only_rho:
                if phiA != 0:
                    fig, axes = plt.subplots(2, 3, figsize=(9, 5), constrained_layout=True, sharex=True, sharey=True)
                    axes[0, 0].tripcolor(tri, rho[:, :, 0].flat, vmin=0, vmax=3 * rho0)
                    axes[1, 0].tripcolor(tri, rho[:, :, 1].flat, vmin=0, vmax=3 * rho0)

                    px = np.zeros((Ly, Lx, 2))
                    py = np.zeros((Ly, Lx, 2))
                    mask = rho != 0
                    px[mask] = mx[mask] / rho[mask]
                    py[mask] = my[mask] / rho[mask]

                    axes[0, 1].tripcolor(tri, px[:, :, 0].flat, vmin=-1, vmax=1, cmap="bwr")
                    axes[1, 1].tripcolor(tri, px[:, :, 1].flat, vmin=-1, vmax=1, cmap="bwr")
                    axes[0, 2].tripcolor(tri, py[:, :, 0].flat, vmin=-1, vmax=1, cmap="bwr")
                    axes[1, 2].tripcolor(tri, py[:, :, 1].flat, vmin=-1, vmax=1, cmap="bwr")

                    for ax in axes.flat:
                        ax.axis("scaled")
                        ax.set_xlim(0, Lx)
                        ax.set_ylim(0, Ly * np.sqrt(3)/2)
                    if savefig:
                        plt.savefig(fout)
                    else:
                        plt.show()
                    plt.close()
                else:
                    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5), constrained_layout=True, sharex=True, sharey=True)
                    im1 = axes[0].tripcolor(tri, rho[:, :, 1].flat, vmin=0, vmax=3 * rho0)

                    px = np.zeros((Ly, Lx, 2))
                    py = np.zeros((Ly, Lx, 2))
                    mask = rho != 0
                    px[mask] = mx[mask] / rho[mask]
                    py[mask] = my[mask] / rho[mask]

                    im2 = axes[1].tripcolor(tri, px[:, :, 1].flat, vmin=-1, vmax=1, cmap="bwr")
                    im3 = axes[2].tripcolor(tri, py[:, :, 1].flat, vmin=-1, vmax=1, cmap="bwr")

                    for ax in axes:
                        ax.axis("scaled")
                        ax.set_xlim(0, Lx)
                        ax.set_ylim(0, Ly * np.sqrt(3)/2)
                    axes[0].set_title(r"$\rho(\mathbf{r})$", fontsize="x-large")
                    axes[1].set_title(r"$p_x(\mathbf{r})$", fontsize="x-large")
                    axes[2].set_title(r"$p_y(\mathbf{r})$", fontsize="x-large")
                    fig.colorbar(im1, ax=axes[0], orientation="horizontal", extend="max")
                    fig.colorbar(im2, ax=axes[1], orientation="horizontal")
                    fig.colorbar(im3, ax=axes[2], orientation="horizontal")
                    fig.suptitle(r"$t=%g$" % ((i_frame+1) * dt + t_beg), fontsize="x-large")

                    if savefig:
                        plt.savefig(fout)
                    else:
                        plt.show()
                    plt.close()
            else:
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5), constrained_layout=True, sharex=True, sharey=True)
                im1 = ax1.tripcolor(tri, rho[:, :, 0].flat, vmax=3 * rho0)
                im2 = ax2.tripcolor(tri, rho[:, :, 1].flat, vmax=3 * rho0)
                ax1.axis("scaled")
                ax2.axis("scaled")
                ax1.set_xlim(0, Lx)
                ax1.set_ylim(0, Ly * np.sqrt(3)/2)

                fig.colorbar(im1, ax=ax1, orientation="horizontal", extend="max")
                fig.colorbar(im2, ax=ax2, orientation="horizontal", extend="max")
                ax1.set_title(r"$\rho_A(\mathbf{r})$", fontsize="x-large")
                ax2.set_title(r"$\rho_B(\mathbf{r})$", fontsize="x-large")

                if savefig:
                    plt.savefig(fout)
                else:
                    plt.show()
                plt.close()


def update_figs():
    folder = f"{root_sohrab}/lat_NRQS/256_256_pA0_q6"
    files = glob.glob(f"{folder}/*")
    for fin in files:
        show_q6(fin, only_rho=False, savefig=True)


if __name__ == "__main__":
    # show_q6(only_rho=False, savefig=True)
    update_figs()