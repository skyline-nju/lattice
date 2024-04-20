import numpy as np
import matplotlib.pyplot as plt
import struct
import platform

root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"

def show_q6(only_rho=False):
    Lx = 256
    Ly = 256
    
    rho0 = 10
    phiA = rho0
    phiB = rho0

    Dt = 0.1
    Dr = 0.1
    v0 = 1

    etaAA = -1
    etaBB = -2
    etaAB = 1
    etaBA = -etaAB
    
    dt = 5000
    t_beg = 0
    seed = 1000
    if (platform.system() == "Linux"):
        folder = f"{root_sohrab}/lat_NRQS/{Lx:d}_{Ly:d}_q4"
    elif (platform.system() == "Windows"):
        folder = "data"
    fin = "%s/L%d_%d_Dr%g_Dt%g_e%g_%g_J%g_%g_v%g_r%g_s%d_dt%d_t%d.bin" % (folder, Lx, Ly, Dr, Dt, etaAA, etaBB, etaAB, etaBA, v0, rho0, seed, dt, t_beg)

    print(fin)
    frame_size = Lx * Ly * 12 * 2
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
            data = np.array(struct.unpack("%dH" % (Lx * Ly * 12), buf)).reshape(Ly, Lx, 2, 6)
            rho = np.sum(data, axis=3)
            # mx = data[:, :, :, 0] - data[:, :, :, 2]
            # my = data[:, :, :, 1] - data[:, :, :, 3]
            print("frame", f.tell()//frame_size, "max particle number:", data.max())
 
            if not only_rho:
                fig, axes = plt.subplots(2, 3, figsize=(9, 3), constrained_layout=True, sharex=True, sharey=True)
                axes[0, 0].imshow(rho[:, :, 0], origin="lower", vmin=0, vmax=rho0*3)
                axes[1, 0].imshow(rho[:, :, 1], origin="lower", vmin=0, vmax=rho0*3)
                # axes[0, 1].imshow(mx[:, :, 0], origin="lower", vmin=-1, vmax=1, cmap="bwr")
                # axes[1, 1].imshow(mx[:, :, 1], origin="lower", vmin=-1, vmax=1, cmap="bwr")
                # axes[0, 2].imshow(my[:, :, 0], origin="lower", vmin=-1, vmax=1, cmap="bwr")
                # axes[1, 2].imshow(my[:, :, 1], origin="lower", vmin=-1, vmax=1, cmap="bwr")
                plt.show()
                plt.close()
            else:
                fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 5), constrained_layout=True, sharex=True, sharey=True)

                x = np.arange(Lx) + 0.5
                y = (np.arange(Ly) + 0.5) * np.sqrt(3) / 2

                xx, yy = np.meshgrid(x, y)

                for i in range(y.size):
                    xx[i] -= i * 0.5
                    mask = xx[i] < 0
                    xx[i][mask] += Lx

                ax1.scatter(xx, yy, s=5, c=rho[:, :, 0], marker="o")
                ax2.scatter(xx, yy, s=5, c=rho[:, :, 1], marker="o")
                # ax2.contourf(xx, yy, rho[:, :, 0], vmin=0, vmax=30)
                ax1.axis("scaled")
                ax2.axis("scaled")
                ax1.set_xlim(0, Lx)
                ax1.set_ylim(0, Ly * np.sqrt(3)/2)


                # im1 = ax1.imshow(rho[:, :, 0], origin="lower", vmin=0, vmax=rho0*3)
                # im2 = ax2.imshow(rho[:, :, 1], origin="lower", vmin=0, vmax=rho0*3)
                # fig.colorbar(im1, ax=ax1, orientation="horizontal", extend="max")
                # fig.colorbar(im2, ax=ax2, orientation="horizontal", extend="max")
                # ax1.set_title(r"$\rho_A(\mathbf{r})$", fontsize="x-large")
                # ax2.set_title(r"$\rho_B(\mathbf{r})$", fontsize="x-large")

                plt.show()
                plt.close()


if __name__ == "__main__":
    show_q6(True)