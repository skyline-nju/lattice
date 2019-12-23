import numpy as np
import matplotlib.pyplot as plt
import struct
import os


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


if __name__ == "__main__":
    f1 = "NM4_Lx200_Ly200_a0.001_p0.01.bin"
    frames = read_profile_2d(f1)
    x = np.arange(0, 1000)

    for frame in frames:
        t, profile1, profile2 = frame
        n1 = count_empty_sites(profile1)
        n2 = count_empty_sites(profile2)
        print(t, n1, n2)
        # plt.plot(x, profile1, x, profile2)
        # plt.show()
        # plt.close()
