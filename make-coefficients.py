import sys, struct, re
src, out = sys.argv[1], sys.argv[2]
sizes = [("IKIM",60),("DUD",300),("FAM",168),("SYS",864),("FAKP",2784),("FAKABP",12),("XFM3CF",882),
         ("ANEW",10),("FAKMAP",464),("ABMAP",6),("F2D",576),("PERR",216),("XESMCF",854),("XPMAP",928),
         ("XESLCF",550),("XESUCF",550),("XERCOF",396)]
def fixed(line, w):
    line = line.rstrip("\n")
    return [line[i:i+w] for i in range(0, len(line), w) if line[i:i+w].strip()]
fof2 = open(src + "/fof2CCIR.asc").read().split("\n")[1:]
f2 = {}
i = 0
for m in range(1, 13):
    assert int(fof2[i][6:8]) == m; i += 1
    vals = []
    while len(vals) < 1976:
        vals += [float(x) for x in fixed(fof2[i], 15)]; i += 1
    assert len(vals) == 1976
    f2[m] = vals
blob = b""
for m in range(1, 13):
    lines = open(src + "/coeff%02d.asc" % m).read().split("\n")
    assert int(lines[0][6:8]) == m
    i = 1
    for name, n in sizes:
        assert lines[i].startswith(name), (lines[i], name); i += 1
        vals = []
        w = 5 if name == "IKIM" else 15
        while len(vals) < n:
            vals += fixed(lines[i], w); i += 1
        assert len(vals) == n, (name, len(vals))
        if name == "IKIM": blob += struct.pack("<%di" % n, *[int(v) for v in vals])
        else: blob += struct.pack("<%df" % n, *[float(v) for v in vals])
    blob += struct.pack("<1976f", *f2[m])
open(out, "wb").write(blob)
print(len(blob), len(blob) / 12 / 4)
