import numpy as np
import matplotlib.pyplot as plt
import struct


if __name__ == "__main__":
    Lx = 128
    Ly = 16
    
    rho0 = 80
    phi = rho0

    Dt = 0.3
    Dr = 0.01
    v0 = 1

    eta = 2
    alpha = 1
    
    dt = 1000
    t_beg = 0
    seed = 1001
    fin = "data/L%d_%d_Dr%g_Dt%g_e%g_a%g_v%g_r%g_s%d_dt%d_t%d.bin" % (Lx, Ly, Dr, Dt, eta, alpha, v0, rho0, seed, dt, t_beg)

    print(fin)
    frame_size = Lx * Ly * 2 * 2
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
            i_frame = f.tell() // frame_size
            print("frame", i_frame)
            buf = f.read(frame_size)
            data = np.array(struct.unpack("%dH%dh" % (Lx * Ly, Lx * Ly), buf)).reshape(2, Ly, Lx)
            rho, m = data
            fig, axes = plt.subplots(2, 2, figsize=(9, 3), constrained_layout=True, sharex=True)
            
            p = np.zeros((Ly, Lx))
            mask = rho > 0
            p[mask] = m[mask] / rho[mask]
            extent = [0, Lx, 0, Ly]
            im1 = axes[0, 0].imshow(rho, origin="lower", vmin=0, vmax=rho0*2, extent=extent)
            im2 = axes[0, 1].imshow(p, origin="lower", vmin=-1, vmax=1, cmap="bwr", extent=extent)

            fig.colorbar(im1, ax=axes[0, 0], orientation="horizontal", extend="max")
            fig.colorbar(im2, ax=axes[0, 1], orientation="horizontal")

            x = np.arange(Lx) + 0.5
            rho_x = np.mean(rho, axis=0)
            m_x = np.mean(m, axis=0)
            axes[1, 0].plot(x, rho_x)
            axes[1, 1].plot(x, m_x/rho_x)
            axes[0, 0].set_title(r"density")
            axes[0, 1].set_title(r"polarity")
            axes[1, 0].set_xlabel(r"$x$")
            axes[1, 1].set_xlabel(r"$x$")
            axes[1, 1].set_ylim(-1, 1)

            title = r"$\eta=%g,D_r=%g, D_t=%g, L_x=%g, L_y=%g$" % (eta, Dr, Dt, Lx, Ly)
            fig.suptitle(title, fontsize="x-large")
            plt.show()
            plt.close()
            
