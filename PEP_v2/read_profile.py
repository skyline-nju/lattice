import numpy as np
import matplotlib.pyplot as plt
import struct
import os
import glob


def read_profile_2d(fin):
    basename = os.path.basename(fin)
    str_list = basename.rstrip(".bin").split("_")
    Lx = int(str_list[1].lstrip("Lx"))
    with open(fin, "rb") as f:
        f.seek(0, 2)
        filesize = f.tell()
        f.seek(0, 0)
        while f.tell() < filesize:
            buf = f.read(4)
            t, = struct.unpack("i", buf)
            buf = f.read(Lx * 2)
            data = struct.unpack("%dB" % (Lx * 2), buf)
            profile1, profile2 = np.array(data).reshape(2, Lx)
            yield t, profile1, profile2


def count_empty_sites(profile):
    count = 0
    for i in profile:
        if i == 0:
            count += 1
    return count


def cal_order_para(Lx, Ly, phi=None, alpha=None, N=2, d=2):
    os.chdir(r"D:/data/PEP/d=%d" % d)
    if phi is not None:
        pat = "NM%d_Lx%d_Ly%d_a*_p%g.bin" % (N, Lx, Ly, phi)
        const_para = "phi"
    elif alpha is not None:
        pat = "NM%d_Lx%d_Ly%d_a%g_p*.bin" % (N, Lx, Ly, alpha)
        const_para = "alpha"
    else:
        const_para = ""
    files = glob.glob(pat)
    for f in files:
        fout = open("order_para" + os.path.sep + f.replace("bin", "dat"), "w")
        n_sum = 0
        count = 0
        frames = read_profile_2d(f)
        for frame in frames:
            t, profile1, profile2 = frame
            n1 = count_empty_sites(profile1)
            n2 = count_empty_sites(profile2)
            fout.write("%d\t%d\t%d\n" % (t, n1, n2))
            if t >= 400000:
                n_sum += (n1 + n2)
                count += 1
        fout.close()
        order_para = n_sum / count / 2 / Lx
        if const_para == "phi":
            alpha = float(f.split("_")[3].lstrip("a"))
            print("%g\t%g" % (alpha, order_para))
        elif const_para == "alpha":
            phi = float(f.split("_")[4].lstrip("p").rstrip(".bin"))
            print("%g\t%g" % (phi, order_para))


if __name__ == "__main__":
    # f1 = "NM4_Lx200_Ly200_a0.001_p0.01.bin"
    # frames = read_profile_2d(f1)
    # x = np.arange(0, 1000)

    # for frame in frames:
    #     t, profile1, profile2 = frame
    #     n1 = count_empty_sites(profile1)
    #     n2 = count_empty_sites(profile2)
    #     print(t, n1, n2)
    #     # plt.plot(x, profile1, x, profile2)
    #     # plt.show()
    #     # plt.close()
    cal_order_para(1000, 12000, alpha=0.0024)
