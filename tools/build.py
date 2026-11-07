#!/usr/bin/env python3
"""Une src/*.lua en dist/ControlAuto.lua y genera dist/ControlAuto.tns"""
import os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import make_tns

parts = ["core.lua", "modules.lua", "ui.lua"]
out = ["platform.apilevel = '2.0'\n",
       "-- Control Automatico - asistente de examen (TI-Nspire CX II)\n",
       "-- Generado por tools/build.py; editar los archivos en src/\n"]
for p in parts:
    out.append("\n-- ===== %s =====\n" % p)
    out.append(open(os.path.join(ROOT, "src", p), encoding="utf-8").read())
src = "".join(out)
os.makedirs(os.path.join(ROOT, "dist"), exist_ok=True)
lua_path = os.path.join(ROOT, "dist", "ControlAuto.lua")
open(lua_path, "w", encoding="utf-8").write(src)
tns = make_tns.build_tns([
    ("Document.xml", make_tns.DEFAULT_DOCUMENT_XML, 0x0D),
    ("Problem1.xml", make_tns.lua_to_problem(src.encode("utf-8")), 0x0D),
])
tns_path = os.path.join(ROOT, "dist", "ControlAuto.tns")
open(tns_path, "wb").write(tns)
print("OK:", lua_path, "|", tns_path, len(tns), "bytes")
