"""Ejecuta los modulos de calculo con casos de las guias (Lua 5.1 via lupa)."""
import sys, os
from lupa import lua51
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
L = lua51.LuaRuntime(unpack_returned_tuples=True)
for f in ("core.lua", "modules.lua"):
    L.execute(open(os.path.join(ROOT, "src", f), encoding="utf-8").read())
MOD = {"spec": 1, "polos": 2, "error": 3, "routh": 4, "freq": 5, "place": 6, "emp": 7}

def run(mod, **vals):
    m = L.globals().MODULES[MOD[mod]]
    t = L.table_from(vals)
    for f in m.fields.values():
        if len(f) >= 4 and f[1] not in vals:
            t[f[1]] = f[4]
    ok, res = L.globals().pcall(m.run, t)
    if not ok:
        return ["ERROR: " + str(res)]
    return list(res.values())

if __name__ == "__main__":
    import json
    cases = json.loads(sys.argv[1])
    for c in cases:
        mod = c.pop("mod")
        print("#####", mod, c)
        for line in run(mod, **c):
            print("  ", line)
