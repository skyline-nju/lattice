# plot phase diagram in the composition plane
import os
import numpy as np
import matplotlib.pyplot as plt
import struct

root_sohrab = "/run/user/1148/gvfs/sftp:host=sohrab003,user=yduan/scratch03.local/yduan"
root_tahmineh = "/run/user/1148/gvfs/sftp:host=tahmineh002,user=yduan/scratch03.local/yduan"


def plot_one_panel(fname, Lx, Ly, i_frame=-1, rho_max=None, ax=None, ms=0.05):
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=(6, 6), constrained_layout=True)
        flag_show = True
    else:
        flag_show = False
    frame_size = Lx * Ly * 8 * 2
    with open(fname, "rb") as f:
        f.seek(0, 2)
        filesize = f.tell()
        n_frames = filesize//frame_size
        print("find", n_frames, "frames")
        if i_frame == -1:
            i_frame += n_frames
        f.seek(frame_size * i_frame)

        buf = f.read(frame_size)
        data = np.array(struct.unpack("%dH" % (Lx * Ly * 8), buf)).reshape(Ly, Lx, 2, 4)
        rho = np.sum(data, axis=3)
        mx = data[:, :, :, 0] - data[:, :, :, 2]
        my = data[:, :, :, 1] - data[:, :, :, 3]
        ax.imshow(rho[:, :, 1], origin="lower", vmin=0, vmax=rho_max, extent=[0, Lx, 0, Ly])

        if flag_show:
            plt.show()
            plt.close()


def phiA_phiB_plane(phiA_arr, phiB_arr, Lx, Ly, Dr, Dt, etaAA, etaBB, etaAB, etaBA, v0, rho0, seed, dn_out=50000, t0=0):
    ncols = phiA_arr.size
    nrows = phiB_arr.size
    folder = f"{root_sohrab}/lat_NRQS/{Lx:d}_{Ly:d}_q4"
    figsize = (ncols * 1 + 0.25, nrows * 1 + 0.25)
    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=True, sharey=True)
    for j, pB in enumerate(phiB_arr[::-1]):
        for i, pA in enumerate(phiA_arr):
            fin = f"{folder}/L{Lx:d}_{Ly:d}_Dr{Dr:g}_Dt{Dt:g}_e{etaAA:g}_{etaBB:g}_J{etaAB}_{etaBA}_v{v0:g}_r{rho0:g}_p{pA:g}_{pB:g}_s{seed:d}_dt{dn_out:d}_t{t0:d}.bin"
            if os.path.exists(fin):
                print(os.path.basename(fin))
                plot_one_panel(fin, Lx, Ly, ax=axes[j, i], rho_max=rho0*3)
            axes[j, i].axis("off")
    rect = [-0.01, -0.01, 1.01, 1.01]
    plt.tight_layout(h_pad=0.1, w_pad=0.1, pad=1.1, rect=rect)
    plt.show()
    plt.close()



if __name__ == "__main__":
    Lx = 128
    Ly = 128
    Dr = 0.1
    Dt = 0.1
    etaAA = etaBB = -1
    etaAB = 1
    etaBA = -etaAB
    v0 = 1
    rho0 = 10
    seed = 1000
    dt = 50000
    t0 = 0

    phiA_arr = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])
    phiB_arr = np.array([0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20])

    phiA_phiB_plane(phiA_arr, phiB_arr, Lx, Ly, Dr, Dt, etaAA, etaBB, etaAB, etaBA, v0, rho0, seed, dt, t0)
    # folder = f"{root_sohrab}/lat_NRQS/{Lx:d}_{Ly:d}_q4"
    # fname = f"{folder}/L128_128_Dr0.1_Dt0.1_e-1_-1_J1_-1_v1_r10_p6_10_s1000_dt50000_t0.bin"
    # plot_one_panel(fname, Lx, Ly)