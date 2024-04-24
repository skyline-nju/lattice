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
    eta = float(s[4].lstrip("e"))
    rho0 = float(s[5].lstrip("r"))
    phi = float(s[6].lstrip("p"))
    v0 = float(s[7].lstrip("v"))
    seed = int(s[8].lstrip("s"))
    dt = int(s[9].lstrip("dt"))
    t0 = int(s[10].lstrip("t"))
    return Lx, Ly, Dr, Dt, eta, v0, rho0, phi, seed, dt, t0

def show_q6(fin=None, savefig=False):
    if fin is None:
        Lx = 128
        Ly = 128
        
        rho0 = 5
        phi = 5

        Dt = 0.07
        Dr = 0.05
        v0 = 1.76
        eta = -2
        
        dt = 5000
        t_beg = 0
        seed = 1000
        if (platform.system() == "Linux"):
            folder = f"{root_sohrab}/latQS/{Lx:d}_{Ly:d}_q6"
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

        while f.tell() < filesize:
            i_frame = f.tell() // frame_size
            if savefig:
                fout = f"{outfolder}/{i_frame:05d}.jpg"
                if os.path.exists(fout):
                    f.seek(frame_size, 1)
                    continue
            buf = f.read(frame_size)
            data = np.array(struct.unpack("%dB" % (Lx * Ly * 6), buf)).reshape(Ly, Lx, 6)
            rho = np.sum(data, axis=2)
            mx = data[:, : , 0] + 0.5 * data[:, :, 1] - 0.5 * data[:, :, 2] - data[:, :, 3] - 0.5 * data[:, :, 4] + 0.5 * data[:, :, 5]
            my = np.sqrt(3) / 2 * (data[:, :, 1] + data[:, :, 2] - data[:, :, 4] - data[:, :, 5])
            print("frame", f.tell()//frame_size, "max particle number:", data.max())
 
            fig, axes = plt.subplots(1, 3, figsize=(12, 4.5), constrained_layout=True, sharex=True, sharey=True)
            im1 = axes[0].tripcolor(tri, rho.flat, vmin=0, vmax=3 * rho0)

            px = np.zeros((Ly, Lx))
            py = np.zeros((Ly, Lx))
            mask = rho != 0
            px[mask] = mx[mask] / rho[mask]
            py[mask] = my[mask] / rho[mask]

            im2 = axes[1].tripcolor(tri, px[:, :].flat, vmin=-1, vmax=1, cmap="bwr")
            im3 = axes[2].tripcolor(tri, py[:, :].flat, vmin=-1, vmax=1, cmap="bwr")

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


def update_figs():
    folder = f"{root_sohrab}/latQS/256_256_q6"
    files = glob.glob(f"{folder}/*")
    for fin in files:
        show_q6(fin, only_rho=False, savefig=True)


if __name__ == "__main__":
    show_q6(savefig=False)
    # update_figs()