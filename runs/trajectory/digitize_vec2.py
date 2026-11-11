"""Extract the solid (undashed) curve(s) of a vector PDF plot panel.
Axis calibration from the centres of the tick-label words (labels centred on their ticks).
usage: digitize_vec2.py pdf page "frame x0,y0,x1,y1" "xlab box" "ylab box" out.csv [yscale]"""
import sys, pymupdf, numpy as np
src, pno = sys.argv[1], int(sys.argv[2])
fr = pymupdf.Rect(*map(float, sys.argv[3].split(",")))
xb = pymupdf.Rect(*map(float, sys.argv[4].split(",")))
yb = pymupdf.Rect(*map(float, sys.argv[5].split(",")))
out = sys.argv[6]
p = pymupdf.open(src)[pno - 1]
xs, ys = [], []
for w in p.get_text("words"):
    r = pymupdf.Rect(w[:4]); c = (r.tl + r.br) / 2
    try:
        v = float(w[4])
    except ValueError:
        continue
    if xb.contains(c):
        xs.append((c.x, v))
    elif yb.contains(c):
        ys.append((c.y, v))
xs, ys = np.array(sorted(xs)), np.array(sorted(ys))
fx = np.polyfit(xs[:, 0], xs[:, 1], 1); fy = np.polyfit(ys[:, 0], ys[:, 1], 1)
rx = np.abs(np.polyval(fx, xs[:, 0]) - xs[:, 1]).max(); ry = np.abs(np.polyval(fy, ys[:, 0]) - ys[:, 1]).max()
print("x labels", xs[:, 1].tolist(), "max resid", rx); print("y labels", ys[:, 1].tolist(), "max resid", ry)
pts = []
for g in p.get_drawings():
    if str(g.get("dashes")) != "[] 0" or not fr.intersects(g["rect"]):
        continue
    for it in g["items"]:
        if it[0] != "l":
            continue
        a, b = it[1], it[2]
        if abs(a.x - b.x) < 1e-3 and (abs(a.x - fr.x0) < 1 or abs(a.x - fr.x1) < 1):
            continue  # frame / ticks on vertical frame lines
        if abs(a.y - b.y) < 1e-3 and (abs(a.y - fr.y0) < 1 or abs(a.y - fr.y1) < 1):
            continue
        if (abs(a.x - b.x) < 1e-3 or abs(a.y - b.y) < 1e-3) and abs(a - b) < 4 and (
                min(abs(a.x - fr.x0), abs(a.x - fr.x1)) < 4.5 or min(abs(a.y - fr.y0), abs(a.y - fr.y1)) < 4.5):
            continue  # tick marks
        for q in (a, b):
            if fr.x0 + 0.3 < q.x < fr.x1 - 0.3 and fr.y0 + 0.3 < q.y < fr.y1 - 0.3:
                pts.append((q.x, q.y))
P = np.unique(np.round(np.array(pts), 3), axis=0)
t, v = np.polyval(fx, P[:, 0]), np.polyval(fy, P[:, 1])
o = np.argsort(t)
np.savetxt(out, np.c_[t[o], v[o]], delimiter=",", fmt="%.5g",
           header=f"x,y: solid vector paths of {src.split('/')[-1]} page {pno}; frame {list(fr)}; x=fx(px) {fx.tolist()}, y=fy(py) {fy.tolist()}; label-fit residual x {rx:.3g} y {ry:.3g}")
print(len(t), "points; x", t.min(), t.max(), "y", v.min(), v.max())
