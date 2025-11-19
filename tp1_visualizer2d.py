#!/usr/bin/env python3
"""
tp1_visualizador2d.py
Visualizador 2D para TP1 (leitura VTK legacy ASCII simples, transformacoes,
projecao ortografica, crescimento incremental e recorte Cohen-Sutherland).
Uso:
    python tp1_visualizador2d.py [arquivo.vtk]
Se arquivo não for informado, gera um exemplo sintético.
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
import math
import re

# ---------------------------
# Util: Leitor simplificado VTK legacy ASCII (PolyData style)
# ---------------------------
def read_vtk_legacy_2d(path):
    """
    Lê um VTK legacy ASCII básico contendo:
      - POINTS n float
      - LINES m indices...
      - POSSÍVEL: CELL_DATA or POINT_DATA with a scalar array 'radius' (tenta encontrar)
    Retorna: points (Nx3), lines (list of (i0,i1)), radius_per_line (list) or None
    Observação: Variações do formato existem; este leitor busca padrões comuns.
    """
    with open(path, 'r') as f:
        text = f.read()
    tokens = text.split()
    # naive parse: find "POINTS"
    pts = []
    lines = []
    radii_line = None

    m = re.search(r'POINTS\s+(\d+)\s+(\w+)', text)
    if m:
        npts = int(m.group(1))
        # find where POINTS line starts, then read following numbers
        idx = text.find(m.group(0))
        after = text[idx + len(m.group(0)):]
        nums = re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', after)
        # take first 3*npts numbers
        nums = nums[:3*npts]
        pts = np.array(nums, dtype=float).reshape((-1,3))
    else:
        raise ValueError("Formato VTK sem 'POINTS' detectado.")

    # find LINES or CELLS (lines could be in various forms)
    m2 = re.search(r'\nLINES\s+(\d+)\s+(\d+)\s', text)
    if m2:
        nlines = int(m2.group(1))
        # locate start
        idx = text.find(m2.group(0))
        after = text[idx + len(m2.group(0)):]
        nums = re.findall(r'\d+', after)
        # each line: k i1 i2 ... where k typically 2 for segments
        nums = nums[:sum([int(nums[i]) + 1 for i in range(0, len(nums), int(nums[0])+1)])] if False else nums
        # simpler: parse sequentially nlines groups
        p = 0
        for _ in range(nlines):
            k = int(nums[p]); p+=1
            if k < 2:
                # ignore degenerate
                p += k
                continue
            indices = [int(nums[p+i]) for i in range(k)]
            p += k
            # if k>2, split into segments between consecutive indices
            if k == 2:
                lines.append((indices[0], indices[1]))
            else:
                for i in range(k-1):
                    lines.append((indices[i], indices[i+1]))
    else:
        # try CELLS with VTK polygon/lines inside
        m3 = re.search(r'\nCELL_DATA\s+(\d+)', text)
        if m3:
            # fallback: look for lines as pairs of consecutive points in points list
            # This is a heuristic
            for i in range(len(pts)-1):
                lines.append((i, i+1))
        else:
            # fallback: no explicit lines -> connect points in order
            for i in range(len(pts)-1):
                lines.append((i, i+1))

    # attempt to find radius: search for "radius" or "radii" in the file and read next numbers array
    radius_per_point = None
    # search for scalar array name radius in POINT_DATA or CELL_DATA
    mrad = re.search(r'(?:POINT_DATA|CELL_DATA)[\s\S]*?(?:SCALARS|LOOKUP_TABLE)\s+([A-Za-z_0-9]+)', text, re.IGNORECASE)
    if mrad:
        name = mrad.group(1)
        # try to find array with that name further and get numbers
        # simpler: find all floating numbers after this position limited to number of points
        pos = text.find(mrad.group(0))
        after = text[pos+len(mrad.group(0)):]
        nums = re.findall(r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', after)
        if len(nums) >= len(pts):
            radius_per_point = np.array(nums[:len(pts)], dtype=float)
    # map radius per line (mean of endpoints) if available
    radius_per_line = None
    if radius_per_point is not None:
        radius_per_line = [0.5*(radius_per_point[i]+radius_per_point[j]) for (i,j) in lines]

    return pts, lines, radius_per_line

# ---------------------------
# Geometria 2D: matrizes homogêneas 3x3
# ---------------------------
def translation_matrix(tx, ty):
    return np.array([[1,0,tx],[0,1,ty],[0,0,1]], dtype=float)

def rotation_matrix(theta_rad):
    c = math.cos(theta_rad); s = math.sin(theta_rad)
    return np.array([[c,-s,0],[s,c,0],[0,0,1]], dtype=float)

def scale_matrix(sx, sy):
    return np.array([[sx,0,0],[0,sy,0],[0,0,1]], dtype=float)

def apply_transform(points2d, M):
    """
    points2d: (N,2) array
    M: 3x3 homogeneous transform
    """
    n = points2d.shape[0]
    homog = np.hstack([points2d, np.ones((n,1))])
    transformed = (M @ homog.T).T
    return transformed[:,:2]

# ---------------------------
# Cohen-Sutherland clipping for segments
# ---------------------------
INSIDE = 0  # 0000
LEFT = 1    # 0001
RIGHT = 2   # 0010
BOTTOM = 4  # 0100
TOP = 8     # 1000

def compute_out_code(x, y, xmin, xmax, ymin, ymax):
    code = INSIDE
    if x < xmin: code |= LEFT
    elif x > xmax: code |= RIGHT
    if y < ymin: code |= BOTTOM
    elif y > ymax: code |= TOP
    return code

def cohen_sutherland_clip(x0,y0,x1,y1, xmin,xmax,ymin,ymax):
    out0 = compute_out_code(x0,y0,xmin,xmax,ymin,ymax)
    out1 = compute_out_code(x1,y1,xmin,xmax,ymin,ymax)
    accept = False
    while True:
        if not (out0 | out1):
            accept = True
            break
        elif out0 & out1:
            break
        else:
            # choose one outside
            outcode_out = out0 if out0 else out1
            if outcode_out & TOP:
                x = x0 + (x1-x0)*(ymax-y0)/(y1-y0 if y1!=y0 else 1e-9)
                y = ymax
            elif outcode_out & BOTTOM:
                x = x0 + (x1-x0)*(ymin-y0)/(y1-y0 if y1!=y0 else 1e-9)
                y = ymin
            elif outcode_out & RIGHT:
                y = y0 + (y1-y0)*(xmax-x0)/(x1-x0 if x1!=x0 else 1e-9)
                x = xmax
            elif outcode_out & LEFT:
                y = y0 + (y1-y0)*(xmin-x0)/(x1-x0 if x1!=x0 else 1e-9)
                x = xmin
            # replace
            if outcode_out == out0:
                x0,y0 = x,y
                out0 = compute_out_code(x0,y0,xmin,xmax,ymin,ymax)
            else:
                x1,y1 = x,y
                out1 = compute_out_code(x1,y1,xmin,xmax,ymin,ymax)
    if accept:
        return True, (x0,y0,x1,y1)
    else:
        return False, None

# ---------------------------
# Visualizador Matplotlib interativo
# ---------------------------
class Visualizador2D:
    def __init__(self, points, lines, radii=None):
        # store original geometry (2D)
        self.orig_pts = points[:,:2].copy()
        self.lines = lines
        self.radii = radii
        # transform parameters
        self.tx = 0.0; self.ty = 0.0
        self.theta = 0.0   # radians
        self.s = 1.0
        # clip window (in world coords)
        self.clip_enabled = False
        self.clip_rect = None  # (xmin,xmax,ymin,ymax) computed from view
        # growth control
        self.n_segments_draw = len(lines)
        self.max_segments = len(lines)

        # figure
        self.fig, self.ax = plt.subplots(figsize=(14, 7))
        plt.subplots_adjust(bottom=0.25)
        self.ax.set_aspect('equal', 'box')
        self.ax.set_title("Visualizador TP1 — teclas: arrows(transl), r/R(rot), +/- (scale), c(toggle clip)")

        # initial draw
        self.line_collection = []
        self._draw()

        # slider para crescimento
        axcolor = 'lightgoldenrodyellow'
        axseg = plt.axes([0.15, 0.1, 0.7, 0.03], facecolor=axcolor)
        self.slider = Slider(axseg, 'Crescimento', 1, self.max_segments, valinit=self.max_segments, valstep=1)
        self.slider.on_changed(self.on_slider)

        # button reset
        axreset = plt.axes([0.8, 0.025, 0.1, 0.04])
        self.btn_reset = Button(axreset, 'Reset')
        self.btn_reset.on_clicked(self.reset)

        # key events
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)

    def _current_transform_matrix(self):
        T = translation_matrix(self.tx, self.ty)
        R = rotation_matrix(self.theta)
        S = scale_matrix(self.s, self.s)
        # order: scale, rotate, translate (world <- T * R * S * p)
        return T @ R @ S

    def _draw(self):
        self.ax.clear()
        self.ax.set_aspect('equal', 'box')
        M = self._current_transform_matrix()
        pts_t = apply_transform(self.orig_pts, M)

        # update clip rect to view limits unless previously set
        xlim = self.ax.get_xlim()
        ylim = self.ax.get_ylim()
        # if axes are empty (first draw) set reasonable limits
        if np.allclose(xlim, (0.0, 1.0)) and np.allclose(ylim, (0.0, 1.0)):
            # compute bounds from pts_t
            minx = pts_t[:,0].min(); maxx = pts_t[:,0].max()
            miny = pts_t[:,1].min(); maxy = pts_t[:,1].max()
            dx = maxx-minx; dy = maxy-miny
            padx = max(0.1*dx, 0.5)
            pady = max(0.1*dy, 0.5)
            self.ax.set_xlim(minx-padx, maxx+padx)
            self.ax.set_ylim(miny-pady, maxy+pady)
            xlim = self.ax.get_xlim(); ylim = self.ax.get_ylim()

        xmin, xmax = self.ax.get_xlim()
        ymin, ymax = self.ax.get_ylim()
        self.clip_rect = (xmin, xmax, ymin, ymax)

        # draw segments up to n_segments_draw
        count = 0
        for idx, (i,j) in enumerate(self.lines):
            if count >= self.n_segments_draw:
                break
            p0 = pts_t[i]; p1 = pts_t[j]
            x0,y0 = p0; x1,y1 = p1
            if self.clip_enabled:
                ok, seg = cohen_sutherland_clip(x0,y0,x1,y1, xmin,xmax,ymin,ymax)
                if ok:
                    x0,y0,x1,y1 = seg
                    self.ax.plot([x0,x1],[y0,y1], linewidth= max(0.5, (self.radii[idx] if self.radii is not None and idx < len(self.radii) else 0.8)), alpha=0.9)
            else:
                self.ax.plot([x0,x1],[y0,y1], linewidth= max(0.5, (self.radii[idx] if self.radii is not None and idx < len(self.radii) else 0.8)), alpha=0.9)
            count += 1

        # draw clip rectangle if enabled
        if self.clip_enabled:
            xmin,xmax,ymin,ymax = self.clip_rect
            self.ax.add_patch(plt.Rectangle((xmin,ymin), xmax-xmin, ymax-ymin, fill=False, linestyle='--', linewidth=1.0))

        # status text
        stat = f"tx={self.tx:.2f}, ty={self.ty:.2f}, rot={math.degrees(self.theta):.1f}°, scale={self.s:.2f}, segs={self.n_segments_draw}/{self.max_segments}, clip={'ON' if self.clip_enabled else 'OFF'}"
        self.ax.text(0.01, 0.01, stat, transform=self.ax.transAxes, fontsize=8, verticalalignment='bottom', bbox=dict(facecolor='white', alpha=0.6))

        self.fig.canvas.draw_idle()

    def on_slider(self, val):
        self.n_segments_draw = int(val)
        self._draw()

    def reset(self, event):
        self.tx = 0.0; self.ty = 0.0; self.theta = 0.0; self.s = 1.0
        self.slider.set_val(self.max_segments)
        self._draw()

    def on_key(self, event):
        step = 0.2
        rot_step = math.radians(5)
        if event.key == 'left':
            self.tx -= step
        elif event.key == 'right':
            self.tx += step
        elif event.key == 'up':
            self.ty += step
        elif event.key == 'down':
            self.ty -= step
        elif event.key == 'r':
            self.theta += rot_step
        elif event.key == 'R':
            self.theta -= rot_step
        elif event.key == '+':
            self.s *= 1.1
        elif event.key == '-':
            self.s /= 1.1
        elif event.key == 'c':
            self.clip_enabled = not self.clip_enabled
        elif event.key == 'pageup':
            self.n_segments_draw = min(self.max_segments, self.n_segments_draw + 1)
            self.slider.set_val(self.n_segments_draw)
        elif event.key == 'pagedown':
            self.n_segments_draw = max(1, self.n_segments_draw - 1)
            self.slider.set_val(self.n_segments_draw)
        # redraw
        self._draw()

def generate_example_tree(num_segments=64):
    # gera uma "árvore" 2D sintética: zigzag + ramificações simples
    pts = []
    lines = []
    angle = 0.0
    x,y = 0.0,0.0
    pts.append([x,y,0.0])
    idx = 0
    rng = np.random.RandomState(123)
    for s in range(1, num_segments+1):
        # each new segment: step with small random angle and occasional branch
        ang = (rng.rand()-0.5)*math.pi/4 + (0.1*math.sin(s*0.12))
        L = 0.8 + 0.2*rng.rand()
        nx = x + L*math.cos(ang)
        ny = y + L*math.sin(ang)
        pts.append([nx,ny,0.0])
        lines.append((idx, idx+1))
        idx += 1
        x,y = nx,ny
        # occasional branch: create short side edge
        if s%8 == 0:
            bx = x + 0.5*math.cos(ang+math.pi/4)
            by = y + 0.5*math.sin(ang+math.pi/4)
            pts.append([bx,by,0.0])
            lines.append((idx, idx+1))
            idx += 1
    pts = np.array(pts)
    # radii synthetic
    radii = [max(0.5, 2.0 - 0.02*i) for i in range(len(lines))]
    return pts, lines, radii

def main():
    if len(sys.argv) >= 2:
        path = sys.argv[1]
        try:
            pts, lines, radii = read_vtk_legacy_2d(path)
        except Exception as e:
            print("Erro ao ler VTK:", e)
            print("Gerando exemplo sintético.")
            pts, lines, radii = generate_example_tree(64)
    else:
        pts, lines, radii = generate_example_tree(128)

    vis = Visualizador2D(pts, lines, radii)
    plt.show()

if __name__ == "__main__":
    main()
