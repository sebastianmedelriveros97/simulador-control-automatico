"""Renderiza pantallas de la app a PNG (aprox. de la API grafica del Nspire)."""
import os, sys, json
from lupa import lua51
from PIL import Image, ImageDraw, ImageFont
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
S = 3
W, H = 318, 212
FONTS = "/usr/share/fonts/truetype/dejavu/"
def font(sz, bold):
    f = FONTS + ("DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf")
    try: return ImageFont.truetype(f, int(sz * 1.25 * S))
    except Exception: return ImageFont.load_default()

class GC:
    def __init__(self):
        self.img = Image.new("RGB", (W*S, H*S), "white"); self.d = ImageDraw.Draw(self.img)
        self.c = (0,0,0); self.f = font(10, False); self.sz = 10
    def setFont(self, _, st, sz): self.f = font(sz, st == "b"); self.sz = sz
    def setColorRGB(self, r, g, b): self.c = (int(r), int(g), int(b))
    def fillRect(self, x, y, w, h): self.d.rectangle([x*S, y*S, (x+w)*S, (y+h)*S], fill=self.c)
    def drawRect(self, x, y, w, h): self.d.rectangle([x*S, y*S, (x+w)*S, (y+h)*S], outline=self.c, width=S)
    def drawLine(self, a, b, c, d): self.d.line([a*S, b*S, c*S, d*S], fill=self.c, width=S)
    def _arc(self, x, y, w, h, st, ang, fill):
        a0, a1 = -(st + ang), -st
        if ang < 0: a0, a1 = -st, -(st + ang)
        box = [x*S, y*S, (x+w)*S, (y+h)*S]
        if fill: self.d.pieslice(box, a0, a1, fill=self.c)
        else: self.d.arc(box, a0, a1, fill=self.c, width=S)
    def drawArc(self, *a): self._arc(*a, False)
    def fillArc(self, *a): self._arc(*a, True)
    def getStringWidth(self, s): return self.d.textlength(s, font=self.f) / S
    def getStringHeight(self, s): return {7: 11, 9: 14, 10: 17}.get(self.sz, self.sz + 5)
    def drawString(self, s, x, y, *_): self.d.text((x*S, y*S), s, fill=self.c, font=self.f)

L = lua51.LuaRuntime(unpack_returned_tuples=True)
L.execute("platform={window={width=function() return %d end,height=function() return %d end,invalidate=function() end}} on={}" % (W, H))
L.execute(open(os.path.join(ROOT, "dist", "ControlAuto.lua"), encoding="utf-8").read())
g = L.globals()

def shot(name):
    gc = GC(); g.on.paint(gc); gc.img.save(name); print("->", name)

if __name__ == "__main__":
    out = sys.argv[1]; mod = sys.argv[2]; vals = json.loads(sys.argv[3])
    g.on.charIn(mod)
    # escribir campos en orden
    fields = g.MODULES[int(mod)].fields
    for i in range(1, len(fields) + 1):
        key = fields[i][1]
        g.on.clearKey()
        txt = vals.get(key, fields[i][4] if len(fields[i]) >= 4 and fields[i][4] else "")
        for ch in txt: g.on.charIn(ch)
        g.on.arrowDown()
    g.on.enterKey()
    shot(out + "_res.png")
    g.on.tabKey()
    shot(out + "_graf.png")
