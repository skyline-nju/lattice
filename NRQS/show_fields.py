import numpy as np
import matplotlib.pyplot as plt
import struct
import platform

root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def show_q2():
    Lx = 128
    Ly = 16
    
    rho0 = 10
    phiA = rho0
    phiB = rho0

    Dt = 0.1
    Dr = 0.01
    v0 = 1

    etaAA = 0
    etaBB = 0
    etaAB = 1
    etaBA = -etaAB
    
    dt = 5000
    t_beg = 0
    seed = 1001

    if (platform.system() == "Linux"):
        folder = f"{root_sohrab}/lat_NRQS/{Lx:d}_{Ly:d}"
    elif (platform.system() == "Windows"):
        folder = "data"
    fin = "%s/L%d_%d_Dr%g_Dt%g_e%g_%g_J%g_%g_v%g_r%g_s%d_dt%d_t%d.bin" % (folder, Lx, Ly, Dr, Dt, etaAA, etaBB, etaAB, etaBA, v0, rho0, seed, dt, t_beg)
    frame_size = Lx * Ly * 2 * 2 * 2
    with open(fin, "rb") as f:
        f.seek(0, 2)
        filesize = f.tell()

        n_frames = filesize//frame_size
        print("find", n_frames, "frames")
        if n_frames > 10:
            f.seek(frame_size * (n_frames - 10))
        else:
            f.seek(0)

        while f.tell() < filesize:
            buf = f.read(frame_size)
            data = np.array(struct.unpack("%dH%dh" % (Lx * Ly * 2, Lx * Ly * 2), buf)).reshape(2, Ly, Lx, 2)
            rho, m = data

            fig, axes = plt.subplots(2, 2, figsize=(9, 3), constrained_layout=True, sharex=True, sharey=True)
            
            
            pA = m[:, :, 0]
            pB = m[:, :, 1]

            pA, pB = np.zeros((2, Ly, Lx))
            mask_A = rho[:, :, 0] > 0
            pA[mask_A] = m[:, :, 0][mask_A] / rho[:, :, 0][mask_A]
            mask_B = rho[:, :, 1] > 0
            pB[mask_B] = m[:, :, 1][mask_B] / rho[:, :, 1][mask_B]

            axes[0, 0].imshow(rho[:, :, 0], origin="lower", vmin=0, vmax=rho0*3, aspect="auto")
            axes[1, 0].imshow(rho[:, :, 1], origin="lower", vmin=0, vmax=rho0*3, aspect="auto")
            axes[0, 1].imshow(pA, origin="lower", vmin=-1, vmax=1, cmap="bwr", aspect="auto")
            axes[1, 1].imshow(pB, origin="lower", vmin=-1, vmax=1, cmap="bwr", aspect="auto")

            plt.show()
            plt.close()


def show_q4(only_rho=False):
    Lx = 512
    Ly = 512
    
    rho0 = 10
    phiA = rho0
    phiB = rho0

    Dt = 0.1
    Dr = 0.1
    v0 = 1

    etaAA = -1
    etaBB = -1
    etaAB = 1
    etaBA = -etaAB
    
    dt = 10000
    t_beg = 0
    seed = 1000
    if (platform.system() == "Linux"):
        folder = f"{root_sohrab}/lat_NRQS/{Lx:d}_{Ly:d}_q4"
    elif (platform.system() == "Windows"):
        folder = "data"
    fin = "%s/L%d_%d_Dr%g_Dt%g_e%g_%g_J%g_%g_v%g_r%g_s%d_dt%d_t%d.bin" % (folder, Lx, Ly, Dr, Dt, etaAA, etaBB, etaAB, etaBA, v0, rho0, seed, dt, t_beg)

    print(fin)
    frame_size = Lx * Ly * 8 * 2
    with open(fin, "rb") as f:
        f.seek(0, 2)
        filesize = f.tell()

        n_frames = filesize//frame_size
        print("find", n_frames, "frames")
        if n_frames > 10:
            f.seek(frame_size * (n_frames - 10))
        else:
            f.seek(0)
        # f.seek(0)

        while f.tell() < filesize:
            buf = f.read(frame_size)
            data = np.array(struct.unpack("%dH" % (Lx * Ly * 8), buf)).reshape(Ly, Lx, 2, 4)
            rho = np.sum(data, axis=3)
            mx = data[:, :, :, 0] - data[:, :, :, 2]
            my = data[:, :, :, 1] - data[:, :, :, 3]
            print("frame", f.tell()//frame_size, "max particle number:", data.max())
 
            if not only_rho:
                fig, axes = plt.subplots(2, 3, figsize=(9, 3), constrained_layout=True, sharex=True, sharey=True)
                axes[0, 0].imshow(rho[:, :, 0], origin="lower", vmin=0, vmax=rho0*3)
                axes[1, 0].imshow(rho[:, :, 1], origin="lower", vmin=0, vmax=rho0*3)
                axes[0, 1].imshow(mx[:, :, 0], origin="lower", vmin=-1, vmax=1, cmap="bwr")
                axes[1, 1].imshow(mx[:, :, 1], origin="lower", vmin=-1, vmax=1, cmap="bwr")
                axes[0, 2].imshow(my[:, :, 0], origin="lower", vmin=-1, vmax=1, cmap="bwr")
                axes[1, 2].imshow(my[:, :, 1], origin="lower", vmin=-1, vmax=1, cmap="bwr")
                plt.show()
                plt.close()
            else:
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5), constrained_layout=True, sharex=True, sharey=True)
                im1 = ax1.imshow(rho[:, :, 0], origin="lower", vmin=0, vmax=rho0*3)
                im2 = ax2.imshow(rho[:, :, 1], origin="lower", vmin=0, vmax=rho0*3)
                fig.colorbar(im1, ax=ax1, orientation="horizontal", extend="max")
                fig.colorbar(im2, ax=ax2, orientation="horizontal", extend="max")
                ax1.set_title(r"$\rho_A(\mathbf{r})$", fontsize="x-large")
                ax2.set_title(r"$\rho_B(\mathbf{r})$", fontsize="x-large")

                plt.show()
                plt.close()

            

if __name__ == "__main__":
    print(platform.system())
#    show_q4(only_rho=True)
