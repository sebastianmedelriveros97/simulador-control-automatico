"""Prueba la interfaz completa con una API Nspire simulada (Lua 5.1)."""
import os
from lupa import lua51
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
L = lua51.LuaRuntime(unpack_returned_tuples=True)
L.execute(r'''
DRAWN = {}
platform = {window = {width=function() return 318 end, height=function() return 212 end,
            invalidate=function() end}}
on = {}
toolpalette = {register=function(t) MENU = t end}
GC = {}
function GC:setFont() end
function GC:setColorRGB() end
function GC:fillRect() end
function GC:drawRect() end
function GC:drawLine() end
function GC:drawArc() end
function GC:fillArc() end
function GC:getStringWidth(s) return #s * 6 end
function GC:getStringHeight(s) return 13 end
function GC:drawString(s, x, y) assert(type(s)=="string", "drawString no string"); DRAWN[#DRAWN+1] = s end
''')
L.execute(open(os.path.join(ROOT, "dist", "ControlAuto.lua"), encoding="utf-8").read())
g = L.globals()
def paint():
    L.execute("DRAWN = {}"); g.on.paint(g.GC); return list(g.DRAWN.values())
def typ(text):
    for ch in text: g.on.charIn(ch)
print(paint()[:4])
g.on.charIn("7")                       # modulo 7 (localizacion de polos)
typ("4/((s+1)(s+2))"); g.on.arrowDown()
typ("2"); g.on.arrowDown(); g.on.arrowDown()
typ("16.3")
print(paint()[:12])
g.on.enterKey()
out = paint()
print("\n".join(out))
for i in range(30): g.on.arrowDown()
print("--- scrolled ---"); print("\n".join(paint()[-6:]))
g.on.escapeKey(); g.on.escapeKey()
for i in range(1, 10):
    g.on.charIn(str(i)); g.on.enterKey(); paint(); g.on.escapeKey(); g.on.escapeKey()
# teclas especiales del Nspire
g.on.charIn("6"); typ("1/(s(s+1))"); g.on.backspaceKey(); typ(")"); g.on.arrowLeft(); g.on.deleteKey(); typ(")")
g.on.clearKey(); typ("−")
print("valor tras clear + menos unicode:", repr(L.eval("VALS ~= nil")))
st = g.on.save(); g.on.restore(st)
g.on.mouseDown(10, 40)
# grafico del plano s (modulo 7)
g.on.escapeKey(); g.on.escapeKey(); g.on.charIn("7")
g.on.clearKey(); typ("4/((s+1)(s+2))")
g.on.arrowDown(); g.on.clearKey(); typ("2"); g.on.arrowDown(); g.on.arrowDown(); g.on.clearKey(); typ("16.3"); g.on.enterKey()
g.on.tabKey(); gr = paint(); assert any("plano" in x.lower() for x in gr), gr[:3]
g.on.arrowDown(); paint(); g.on.tabKey(); paint()
print("menu:", len(list(g.MENU.values())))
print("OK UI")
