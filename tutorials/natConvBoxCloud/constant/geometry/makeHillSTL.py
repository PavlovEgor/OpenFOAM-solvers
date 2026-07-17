#!/usr/bin/env python3
"""Generate hill.stl: a symmetric, axisymmetric raised-cosine hill heightfield
used by system/blockMeshDict to project the bottom patch into a genuine 3D
dome (instead of a 2D ridge extruded along Z).

h(r) = H/2 * (1 + cos(pi*r/(W/2)))  for r <= W/2, else 0

where r is the radial distance from the hill centre (xc, zc) in the X-Z
plane. L, Lz, H and W are not duplicated here: they are read straight out of
system/blockMeshDict (the single source of truth) via foamDictionary, so
this script always matches whatever blockMeshDict currently defines.
"""
import math
import os
import struct
import subprocess

CASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BLOCK_MESH_DICT = os.path.join(CASE_DIR, "system", "blockMeshDict")


def read_dict_entry(entry):
    value = subprocess.run(
        ["foamDictionary", "-entry", entry, "-value", BLOCK_MESH_DICT],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    return float(value)


L = read_dict_entry("L")
Lz = read_dict_entry("Lz")
H = read_dict_entry("H")
W = read_dict_entry("W")
xc = L / 2.0
zc = Lz / 2.0

# STL grid resolution (kept finer than the block mesh's N x Nz for a smooth
# projection surface).
nx = 200
nz = 200


def height(x, z):
    r = math.hypot(x - xc, z - zc)
    if r <= W / 2.0:
        return 0.5 * H * (1.0 + math.cos(math.pi * r / (W / 2.0)))
    return 0.0


def normal(a, b, c):
    ux, uy, uz = b[0] - a[0], b[1] - a[1], b[2] - a[2]
    vx, vy, vz = c[0] - a[0], c[1] - a[1], c[2] - a[2]
    nx_ = uy * vz - uz * vy
    ny_ = uz * vx - ux * vz
    nz_ = ux * vy - uy * vx
    mag = math.sqrt(nx_ * nx_ + ny_ * ny_ + nz_ * nz_) or 1.0
    return (nx_ / mag, ny_ / mag, nz_ / mag)


def main():
    pts = [
        [(x, height(x, z), z) for z in (Lz * j / nz for j in range(nz + 1))]
        for x in (L * i / nx for i in range(nx + 1))
    ]

    triangles = []
    for i in range(nx):
        for j in range(nz):
            p00 = pts[i][j]
            p10 = pts[i + 1][j]
            p11 = pts[i + 1][j + 1]
            p01 = pts[i][j + 1]
            triangles.append((p00, p10, p11))
            triangles.append((p00, p11, p01))

    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hill.stl")
    with open(out_path, "wb") as f:
        header = b"hill heightfield STL"
        f.write(header + b" " * (80 - len(header)))
        f.write(struct.pack("<I", len(triangles)))
        for a, b, c in triangles:
            n = normal(a, b, c)
            f.write(struct.pack("<3f", *n))
            f.write(struct.pack("<3f", *a))
            f.write(struct.pack("<3f", *b))
            f.write(struct.pack("<3f", *c))
            f.write(struct.pack("<H", 0))

    print(f"wrote {len(triangles)} triangles to {out_path}")


if __name__ == "__main__":
    main()
