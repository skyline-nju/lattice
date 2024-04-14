import numpy as np
import matplotlib.pyplot as plt
import struct


if __name__ == "__main__":
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
    fin = "data/L%d_%d_Dr%g_Dt%g_e%g_%g_J%g_%g_v%g_r%g_s%d_dt%d_t%d.bin" % (Lx, Ly, Dr, Dt, etaAA, etaBB, etaAB, etaBA, v0, rho0, seed, dt, t_beg)

    print(fin)
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
            
