--------------------------------------------------------------------------
-- modules.lua : modulos de calculo (cada uno devuelve lista de lineas)
--------------------------------------------------------------------------
local abs, sqrt, atan2, log, exp = math.abs, math.sqrt, math.atan2, math.log, math.exp
local pi = math.pi
local INF = 1/0
local D2R, R2D = pi/180, 180/pi

local function blank(s) return s == nil or s:match("^%s*$") ~= nil end

-- numero desde un campo (acepta expresiones: pi/3, sqrt(2), 60*pi/180 ...)
local function num(v, key, name)
  local s = v[key]
  if blank(s) then return nil end
  local r = parse(s, nil)
  if not risconst(r) then error((name or key) .. ": debe ser un numero") end
  return rval(r)
end

local function req(x, name)
  if x == nil then error("Falta dato: " .. name) end
  return x
end

-- funcion de transferencia desde campo (con k opcional)
local function tf(v, key, kval, default)
  local s = v[key]
  if blank(s) then
    if default then s = default else error("Falta dato: " .. key) end
  end
  local ok, r = pcall(parse, s, kval)
  if not ok then
    if tostring(r):find("NEEDK") then error(key .. " usa k: ingrese un valor en el campo k") end
    error(key .. ": " .. tostring(r):gsub("^.-:%d+: ", ""))
  end
  return r
end

local function add(L, s) L[#L + 1] = s end
local function sep(L, t) L[#L + 1] = "== " .. t .. " ==" end

local Z, W, WN, PHI, TAU, SIG, DEG, PM = U.zeta, U.w, U.wn, U.phi, U.tau, U.sigma, U.deg, U.pm

--------------------------------------------------------------------------
-- reporte de especificaciones desde el par/polo dominante
--------------------------------------------------------------------------
local function specs_pair(L, sig, wd)
  local wn = sqrt(sig * sig + wd * wd)
  local z = sig / wn
  add(L, SIG .. " = " .. fmt(sig) .. "   " .. W .. "d = " .. fmt(wd))
  add(L, WN .. " = " .. U.sqrt .. "(" .. SIG .. U.sq .. "+" .. W .. "d" .. U.sq .. ") = " .. fmt(wn) .. "   " .. Z .. " = " .. SIG .. "/" .. WN .. " = " .. fmt(z))
  add(L, "ts(2%) = 4/" .. SIG .. " = " .. fmt(4 / sig) .. " s   ts(5%) = 3/" .. SIG .. " = " .. fmt(3 / sig) .. " s")
  add(L, "Mp = 100e^(-" .. Z .. U.pi .. "/" .. U.sqrt .. "(1-" .. Z .. U.sq .. ")) = " .. fmt(Mp_from_zeta(z)) .. " %")
  add(L, "tp = " .. U.pi .. "/" .. W .. "d = " .. fmt(pi / wd) .. " s   tr = (" .. U.pi .. "-acos" .. Z .. ")/" .. W .. "d = " .. fmt((pi - math.acos(z)) / wd) .. " s")
end

local function specs_real(L, p)
  local tau = 1 / abs(p)
  add(L, "Polo real " .. fmt(p) .. ":  " .. TAU .. " = 1/|p| = " .. fmt(tau) .. " s")
  add(L, "ts(2%) = 4" .. TAU .. " = " .. fmt(4 * tau) .. " s   ts(5%) = 3" .. TAU .. " = " .. fmt(3 * tau) .. " s   Mp = 0 %")
end

-- analisis dinamico (polos dominantes / reduccion de orden)
local function dyn_report(L, T, poles, zeros)
  local _, stab = stabilitytext(poles)
  if not stab then
    add(L, "(no se calculan ts/Mp: el sistema no es estable)")
    return
  end
  local np = #poles
  if np == 0 then add(L, "Ganancia pura (sin polos)"); return end
  local K0 = T.n[1] / T.d[1]
  local dom = poles[1]
  local nd = (dom[2] ~= 0) and 2 or 1
  if np == nd then
    if nd == 1 then specs_real(L, dom[1]) else specs_pair(L, -dom[1], abs(dom[2])) end
  else
    local nxt = poles[nd + 1]
    local ratio = abs(nxt[1]) / abs(dom[1])
    add(L, "Polo dominante: " .. fmtc(dom[1], abs(dom[2])) .. (nd == 2 and (" (par " .. PM .. ")") or ""))
    add(L, "Separacion: |Re sig.|/|Re dom.| = " .. fmt(abs(nxt[1])) .. "/" .. fmt(abs(dom[1])) .. " = " .. fmt(ratio, 4))
    if np == 2 and nd == 1 and ratio < 5 then
      -- segundo orden sobreamortiguado
      local p1, p2 = poles[1][1], poles[2][1]
      local wn = sqrt(p1 * p2)
      local z = -(p1 + p2) / (2 * wn)
      add(L, "2o orden sobreamortiguado: " .. WN .. " = " .. fmt(wn) .. "  " .. Z .. " = " .. fmt(z))
      add(L, "ts(2%) " .. U.approx .. " 4/|p dom| = " .. fmt(4 / abs(p1)) .. " s  (5%: " .. fmt(3 / abs(p1)) .. " s)  Mp = 0 %")
      return
    end
    if ratio >= 5 then
      add(L, "Separacion " .. U.ge .. " 5 " .. U.arrow .. " se puede reducir a orden " .. nd)
      local dp = {dom}
      if nd == 2 then dp[2] = {dom[1], -dom[2]} end
      local den = pfromroots(dp, 1)
      local Kr = K0 * den[1]
      add(L, "Gr(s) = " .. fmt(Kr) .. " / (" .. pstr(den) .. ")   [Gr(0)=G(0)=" .. fmt(K0) .. "]")
    else
      add(L, "Separacion < 5 " .. U.arrow .. " NO conviene reducir (valores solo referenciales):")
    end
    if nd == 1 then specs_real(L, dom[1]) else specs_pair(L, -dom[1], abs(dom[2])) end
  end
  if #zeros > 0 then add(L, "Ojo: hay ceros (" .. rootsstr(zeros, 4) .. "); pueden alterar Mp y tr.") end
end

--------------------------------------------------------------------------
-- 1) Especificaciones de segundo orden
--------------------------------------------------------------------------
local function m_specs(v)
  local L = {}
  local Mp, ts, crit = num(v, "Mp"), num(v, "ts"), num(v, "crit") or 2
  local z, wn, tp, sig, wd = num(v, "zeta"), num(v, "wn"), num(v, "tp"), num(v, "sigma"), num(v, "wd")
  local c = (crit == 5) and 3 or 4
  sep(L, "Datos " .. U.arrow .. " parametros")
  if sig and wd then
    wn = wn or sqrt(sig * sig + wd * wd)
    z = z or sig / wn
    add(L, "Desde polos: " .. WN .. " = " .. fmt(wn) .. ", " .. Z .. " = " .. fmt(z))
  end
  if not z and Mp then
    local l = log(Mp / 100)
    z = -l / sqrt(pi * pi + l * l)
    add(L, Z .. " = -ln(Mp/100)/" .. U.sqrt .. "(" .. U.pi .. U.sq .. "+ln" .. U.sq .. "(Mp/100))")
    add(L, "  = -ln(" .. fmt(Mp / 100) .. ")/" .. U.sqrt .. "(" .. fmt(pi * pi) .. "+" .. fmt(l * l) .. ") = " .. fmt(z))
  end
  if not sig and ts then
    sig = c / ts
    add(L, SIG .. " = " .. Z .. WN .. " = " .. c .. "/ts = " .. c .. "/" .. fmt(ts) .. " = " .. fmt(sig) .. "  (crit " .. (c == 4 and "2" or "5") .. "%)")
  end
  if not wd and tp then
    wd = pi / tp
    add(L, W .. "d = " .. U.pi .. "/tp = " .. fmt(wd))
  end
  if not wn then
    if sig and z then wn = sig / z; add(L, WN .. " = " .. SIG .. "/" .. Z .. " = " .. fmt(wn))
    elseif wd and z and z < 1 then wn = wd / sqrt(1 - z * z); add(L, WN .. " = " .. W .. "d/" .. U.sqrt .. "(1-" .. Z .. U.sq .. ") = " .. fmt(wn))
    elseif sig and wd then wn = sqrt(sig * sig + wd * wd) end
  end
  if not z and sig and wn then z = sig / wn; add(L, Z .. " = " .. SIG .. "/" .. WN .. " = " .. fmt(z)) end
  if not z and wd and wn and wd < wn then z = sqrt(1 - (wd / wn) ^ 2); add(L, Z .. " = " .. U.sqrt .. "(1-(" .. W .. "d/" .. WN .. ")" .. U.sq .. ") = " .. fmt(z)) end
  if not z then error("Faltan datos: ingrese Mp o " .. Z) end
  sep(L, "Resultados")
  add(L, Z .. " = " .. fmt(z) .. "   Mp = " .. fmt(Mp_from_zeta(z)) .. " %")
  if z > 0 and z < 1 then
    add(L, U.phi .. " = acos(" .. Z .. ") = " .. fmt(math.acos(z) * R2D) .. DEG .. "  (angulo del polo)")
  end
  if wn then
    sig = sig or z * wn
    add(L, WN .. " = " .. fmt(wn) .. " rad/s   " .. SIG .. " = " .. Z .. WN .. " = " .. fmt(sig))
    if z < 1 then
      wd = wd or wn * sqrt(1 - z * z)
      add(L, W .. "d = " .. WN .. U.sqrt .. "(1-" .. Z .. U.sq .. ") = " .. fmt(wd) .. " rad/s")
      add(L, "Polos deseados: s = -" .. fmt(sig) .. " " .. PM .. " j" .. fmt(wd))
      add(L, "ts(2%) = " .. fmt(4 / sig) .. " s   ts(5%) = " .. fmt(3 / sig) .. " s")
      add(L, "tp = " .. fmt(pi / wd) .. " s   tr = " .. fmt((pi - math.acos(z)) / wd) .. " s")
    else
      local r = sqrt(z * z - 1)
      add(L, "Polos reales: " .. fmt(-z * wn + wn * r) .. ", " .. fmt(-z * wn - wn * r))
    end
    add(L, "Pol. deseado: s" .. U.sq .. " + " .. fmt(2 * sig) .. "s + " .. fmt(wn * wn))
    add(L, "G std = " .. fmt(wn * wn) .. "K/(s" .. U.sq .. "+" .. fmt(2 * sig) .. "s+" .. fmt(wn * wn) .. ")")
  else
    add(L, "(ingrese ts, " .. WN .. " o tp para obtener polos y tiempos)")
  end
  return L
end

--------------------------------------------------------------------------
-- 2) Polos, lazo cerrado y respuesta temporal
--------------------------------------------------------------------------
local function linstr(c0, c1)
  -- texto de c0 + c1 k
  local a, b = abs(c0) > 1e-12, abs(c1) > 1e-12
  if not a and not b then return "0" end
  local kt = (abs(c1 - 1) < 1e-12) and "k" or (fmt(c1) .. "k")
  if not a then return kt end
  if not b then return fmt(c0) end
  return fmt(c0) .. (c1 < 0 and " - " or " + ") .. ((abs(abs(c1) - 1) < 1e-12) and "k" or (fmt(abs(c1)) .. "k"))
end

local function polyk_str(P0, P1)
  -- polinomio con coeficientes lineales en k
  local n = math.max(#P0, #P1)
  local parts = {}
  for i = n, 1, -1 do
    local c0, c1 = P0[i] or 0, P1[i] or 0
    if abs(c0) > 1e-12 or abs(c1) > 1e-12 then
      local e = i - 1
      local cs = linstr(c0, c1)
      local v = (e == 0) and "" or ((e == 1) and "s" or ("s^" .. e))
      if e > 0 and cs == "1" then cs = "" elseif e > 0 and (abs(c0) > 1e-12 and abs(c1) > 1e-12) then cs = "(" .. cs .. ")" end
      parts[#parts + 1] = cs .. v
    end
  end
  if #parts == 0 then return "0" end
  return table.concat(parts, " + ")
end

local function m_polos_k(v, lazo)
  local L = {}
  local function G0at(kk) return rmul(tf(v, "Gc", kk, "1"), tf(v, "Gp", kk)) end
  local g0, g1, g2 = G0at(0), G0at(1), G0at(2)
  local function parts(f)
    local p0, p1, p2 = f(g0), f(g1), f(g2)
    local d1 = psub(p1, p0)
    local chk = psub(psub(p2, p0), pscale(d1, 2))
    for i = 1, #chk do if abs(chk[i]) > 1e-9 * math.max(1, abs(p2[i] or 0)) then error("k aparece en forma no lineal: ingrese un valor de k") end end
    return p0, d1
  end
  -- normalizar: mismo denominador en los tres (no simplificar)
  local N0, N1 = parts(function(g) return g.n end)
  local D0, D1 = parts(function(g) return g.d end)
  sep(L, "En funcion de k")
  add(L, "G0(s) = [" .. polyk_str(N0, N1) .. "] / [" .. polyk_str(D0, D1) .. "]")
  local C0, C1 = padd(D0, N0), padd(D1, N1)
  if lazo == 1 then
    add(L, "GLC(s) = G0/(1+G0) = [" .. polyk_str(N0, N1) .. "] / [" .. polyk_str(C0, C1) .. "]")
    add(L, "Ec. caract.: " .. polyk_str(C0, C1) .. " = 0")
    local n = math.max(#C0, #C1)
    local lead0, lead1 = C0[n] or 0, C1[n] or 0
    if n == 3 and abs(lead1) < 1e-12 and abs(lead0) > 1e-12 then
      local l = lead0
      add(L, "Comparando con s" .. U.sq .. " + 2" .. Z .. WN .. "s + " .. WN .. U.sq .. ":")
      add(L, "  " .. WN .. U.sq .. " = " .. linstr((C0[1] or 0) / l, (C1[1] or 0) / l) .. "   " .. U.arrow .. " " .. WN .. " = " .. U.sqrt .. "(" .. linstr((C0[1] or 0) / l, (C1[1] or 0) / l) .. ")")
      add(L, "  2" .. Z .. WN .. " = " .. linstr((C0[2] or 0) / l, (C1[2] or 0) / l) .. "   " .. U.arrow .. " " .. Z .. " = (" .. linstr((C0[2] or 0) / l, (C1[2] or 0) / l) .. ")/(2" .. WN .. ")")
      if abs(C1[2] or 0) < 1e-12 then
        add(L, "  " .. Z .. WN .. " = " .. fmt((C0[2] or 0) / l / 2) .. " fijo " .. U.arrow .. " ts no depende de k;")
        add(L, "  al subir k sube " .. WN .. ", baja " .. Z .. " y sube Mp.")
      end
    end
  end
  add(L, "Ingrese k para evaluar polos, ts y Mp.")
  return L
end

local function m_polos(v)
  local lazo0 = num(v, "lazo") or 1
  if num(v, "k") == nil and (usesk(v.Gc) or usesk(v.Gp)) then return m_polos_k(v, lazo0) end
  local L = {}
  local k = num(v, "k")
  local Gc = tf(v, "Gc", k, "1")
  local Gp = tf(v, "Gp", k)
  local lazo = num(v, "lazo") or 1
  local A = num(v, "A") or 1
  local G0 = rsimplify(rmul(Gc, Gp))
  sep(L, "Lazo abierto")
  add(L, "G0 = Gc" .. "\194\183" .. "Gp = " .. rstr(G0))
  local T
  if lazo == 1 then
    T = closedloop(G0)
    sep(L, "Lazo cerrado (retro unitaria)")
    add(L, "GLC = G0/(1+G0) = " .. rstr(T))
    add(L, "Ec. caract.: " .. pstr(T.d) .. " = 0")
  else
    T = G0
    sep(L, "Funcion analizada: G = Gc" .. "\194\183" .. "Gp")
  end
  local poles, zeros = proots(T.d), proots(T.n)
  add(L, "Polos: " .. rootsstr(poles))
  add(L, "Ceros: " .. rootsstr(zeros))
  local st, ok = stabilitytext(poles)
  add(L, st)
  local d0 = T.d[1]
  if abs(d0) > 1e-12 then
    local K0 = T.n[1] / d0
    add(L, "Ganancia estatica G(0) = " .. fmt(K0))
    if ok then
      add(L, "y(" .. U.inf .. ") = A" .. "\194\183" .. "G(0) = " .. fmt(A * K0))
      if lazo == 1 then add(L, "e(" .. U.inf .. ") escalon = A - y(" .. U.inf .. ") = " .. fmt(A - A * K0)) end
    end
  else
    add(L, "G tiene integrador: y(" .. U.inf .. ") no acotada ante escalon")
  end
  sep(L, "Respuesta temporal")
  add(L, "Regla: si |Re| del siguiente polo " .. U.ge .. " 5" .. "\194\183" .. "|Re| del dominante, se desprecia (reduccion).")
  dyn_report(L, T, poles, zeros)
  return L
end

--------------------------------------------------------------------------
-- 2) Respuesta temporal y(t) por fracciones parciales
--------------------------------------------------------------------------
local function expo(pr)
  if abs(pr) < 1e-11 then return "" end
  if abs(pr + 1) < 1e-12 then return "e^(-t)" end
  if abs(pr - 1) < 1e-12 then return "e^(t)" end
  return "e^(" .. fmt(pr) .. "t)"
end

local function m_resp(v)
  local L = {}
  local k = num(v, "k")
  local G = tf(v, "G", k)
  local lazo = num(v, "lazo") or 0
  local A = num(v, "A") or 1
  local ent = num(v, "entrada") or 1
  G = rsimplify(G)
  if lazo == 1 then G = closedloop(G) end
  local R
  if ent == 0 then R = rconst(A) elseif ent == 2 then R = rat({A}, {0, 0, 1}) else R = rat({A}, {0, 1}) end
  local names = {[0] = "impulso R(s) = A", "escalon R(s) = A/s", "rampa R(s) = A/s" .. U.sq}
  sep(L, "Paso 1: Y(s) = G(s)" .. "\194\183" .. "R(s)")
  add(L, (lazo == 1 and "GLC" or "G") .. "(s) = " .. rstr(G))
  add(L, "Entrada: " .. (names[ent] or names[1]) .. ",  A = " .. fmt(A))
  local Y = rat(pmul(G.n, R.n), pmul(G.d, R.d))
  add(L, "Y(s) = " .. rstr(Y))
  local gp = proots(G.d)
  add(L, "Polos de G: " .. rootsstr(gp))
  local groups = partialfractions(Y.n, Y.d)
  sep(L, "Paso 2: Fracciones parciales")
  add(L, "Y(s) = " .. "\206\163" .. " Ai/(s - pi),  Ai = lim(s" .. U.arrow .. "pi) (s - pi)Y(s)")
  local terms = {}
  local function addterm(c, body)
    if abs(c) < 1e-11 then return end
    local sgn = (c < 0) and "-" or "+"
    local mag = abs(c)
    local ms = fmt(mag)
    if body ~= "" and abs(mag - 1) < 1e-12 then ms = "" end
    terms[#terms + 1] = {sgn, ms .. body}
  end
  local function fact(n) local f = 1 for i = 2, n do f = f * i end return f end
  for _, g in ipairs(groups) do
    local pr, pim = g.p[1], g.p[2]
    if pim >= 0 then
      for j = 1, g.m do
        local c = g.c[j]
        local den = "(s" .. ((abs(pr) < 1e-11 and pim == 0) and "" or (" - (" .. fmtc(pr, pim, 4) .. ")")) .. ")" .. (j > 1 and ("^" .. j) or "")
        if abs(pr) < 1e-11 and pim == 0 then den = "s" .. (j > 1 and ("^" .. j) or "") end
        local cs = fmtc(c[1], c[2])
        if c[2] ~= 0 and abs(c[1]) > 1e-11 then cs = "(" .. cs .. ")" end
        add(L, cs .. " / " .. den .. ((pim > 0) and "  (+ conjugado)" or ""))
        local tpow = j - 1
        local tt = (tpow == 0) and "" or ((tpow == 1) and "t" or ("t^" .. tpow))
        local f = fact(tpow)
        if pim == 0 then
          local ex = expo(pr)
          local body = tt .. ((tt ~= "" and ex ~= "") and "\194\183" or "") .. ex
          addterm(c[1] / f, body)
        else
          local a, b = c[1] / f, c[2] / f
          local M, ph = 2 * sqrt(a * a + b * b), atan2(b, a) * R2D
          local ex = expo(pr)
          local body = tt .. ex .. "cos(" .. fmt(pim) .. "t" .. (ph < 0 and " - " or " + ") .. fmt(abs(ph), 4) .. DEG .. ")"
          addterm(M, body)
        end
      end
    end
  end
  sep(L, "Paso 3: Antitransformada y(t)  (t " .. U.ge .. " 0)")
  add(L, "A/(s-p) " .. U.arrow .. " A e^(pt);  par complejo " .. U.arrow .. " 2|A|e^(" .. U.sigma .. "t)cos(" .. W .. "t+" .. U.ang .. "A)")
  local txt = "y(t) ="
  for i, t in ipairs(terms) do
    if i == 1 then txt = txt .. " " .. (t[1] == "-" and "-" or "") .. t[2]
    else txt = txt .. " " .. t[1] .. " " .. t[2] end
  end
  if #terms == 0 then txt = txt .. " 0" end
  add(L, txt)
  -- modos
  sep(L, "Paso 4: Modos de respuesta (polos de G)")
  add(L, "Cada polo de G da un modo e^(pt). Re(p)<0: decae; Re(p)>0: diverge.")
  add(L, "Mas veloz = el que decae primero = mayor |Re(p)|.")
  local fast, slow, diverge = nil, nil, false
  local seen = {}
  for _, r in ipairs(gp) do
    local key = fmt(r[1], 6) .. "," .. fmt(abs(r[2]), 6)
    if r[2] >= 0 and not seen[key] then
      seen[key] = true
      local nm = (r[2] > 0) and (expo(r[1]) .. "cos(" .. fmt(r[2]) .. "t)") or expo(r[1])
      if nm == "" then nm = "constante (polo en 0)" end
      if r[1] > 1e-9 then
        add(L, nm .. ": DIVERGE (polo en SPD) " .. U.arrow .. " y(t)" .. U.arrow .. U.inf)
        diverge = true
      elseif r[1] >= -1e-9 then
        add(L, nm .. ": no decae (polo en eje j" .. W .. ")")
        diverge = true
      else
        local tau = 1 / abs(r[1])
        add(L, nm .. ": converge, " .. TAU .. " = 1/" .. fmt(abs(r[1])) .. " = " .. fmt(tau) .. " s, se extingue en " .. U.approx .. "4" .. TAU .. " = " .. fmt(4 * tau) .. " s" .. ((r[2] > 0) and " (oscilatorio)" or ""))
        if not fast or abs(r[1]) > abs(fast[1]) then fast = r end
        if not slow or abs(r[1]) < abs(slow[1]) then slow = r end
      end
    end
  end
  if fast then add(L, "Modo mas veloz (decae primero): polo " .. fmtc(fast[1], fast[2], 4) .. " (mayor |Re|)") end
  if slow then add(L, "Modo mas lento (dominante, fija ts): polo " .. fmtc(slow[1], slow[2], 4)) end
  if diverge then
    add(L, "SISTEMA INESTABLE: hay modos que no decaen;")
    add(L, "y(t) no alcanza estado estable.")
  else
    add(L, "Sistema estable: ts(2%) " .. U.approx .. " 4/|Re dom| = " .. fmt(4 / abs(slow[1])) .. " s")
    if ent == 1 and abs(G.d[1]) > 1e-12 then add(L, "y(" .. U.inf .. ") = A" .. "\194\183" .. "G(0) = " .. fmt(A * G.n[1] / G.d[1])) end
  end
  return L
end

--------------------------------------------------------------------------
-- 3) Error en estado estacionario
--------------------------------------------------------------------------
local function errlines(L, e, A)
  add(L, "Tipo del sistema N = " .. e.N .. "  (n" .. "\194\186" .. " de polos en s=0 de G0)")
  add(L, "(tipo 0: ep finito; tipo 1: ep=0, ev finito; tipo 2: ev=0, ea finito)")
  add(L, "Kp = lim G0 = " .. fmt(e.Kp) .. "   " .. U.arrow .. " ep = A/(1+Kp) = " .. fmt(e.ep))
  add(L, "Kv = lim sG0 = " .. fmt(e.Kv) .. "   " .. U.arrow .. " ev = A/Kv = " .. fmt(e.ev))
  add(L, "Ka = lim s" .. U.sq .. "G0 = " .. fmt(e.Ka) .. "   " .. U.arrow .. " ea = A/Ka = " .. fmt(e.ea))
  if e.ep ~= INF and e.ep > 0 then add(L, "ep = " .. fmt(100 * e.ep / A) .. " % de la referencia") end
end

local function m_error(v)
  local L = {}
  local A = num(v, "A") or 1
  local k = num(v, "k")
  local eobj = num(v, "eobj")
  local sig = num(v, "senal") or 1
  local needk = (usesk(v.Gc) or usesk(v.Gp)) and k == nil
  if needk then
    if not eobj then error("Hay k sin valor: ingrese k, o 'e obj' para despejar k") end
    local function E(kk)
      local G0 = rmul(tf(v, "Gc", kk, "1"), tf(v, "Gp", kk))
      local e = sserrors(G0, A)
      if sig == 1 then return e.ep elseif sig == 2 then return e.ev else return e.ea end
    end
    local names = {"escalon (ep)", "rampa (ev)", "parabola (ea)"}
    sep(L, "Despeje de k para " .. names[sig] .. " = " .. fmt(eobj))
    local e1, e2 = E(1.37), E(2.91)
    if abs(e1 - e2) < 1e-12 * math.max(1, abs(e1)) then
      add(L, "El error NO depende de k (e = " .. fmt(e1) .. "): no se puede cumplir")
      add(L, "con este parametro.")
      return L
    end
    local ks = {}
    for e = -40, 50 do ks[#ks + 1] = 10 ^ (e / 10) end
    local sols = {}
    for sgn = -1, 1, 2 do
      local prev, pk
      for i = 1, #ks do
        local kk = sgn * ks[i]
        local ok, val = pcall(E, kk)
        if ok and val ~= INF and val == val then
          local f = val - eobj
          if prev and (prev > 0) ~= (f > 0) then
            local a, b, fa = pk, kk, prev
            for it = 1, 80 do
              local c = (a + b) / 2
              local okc, vc = pcall(E, c)
              local fc = (okc and vc ~= INF) and (vc - eobj) or fa
              if (fc > 0) == (fa > 0) then a, fa = c, fc else b = c end
            end
            sols[#sols + 1] = (a + b) / 2
          end
          prev, pk = f, kk
        else
          prev = nil
        end
      end
    end
    if #sols == 0 then add(L, "No se encontro k en [-1e5, 1e5]") end
    for _, kk in ipairs(sols) do
      local G0 = rmul(tf(v, "Gc", kk, "1"), tf(v, "Gp", kk))
      local T = closedloop(G0)
      local st = stabilitytext(proots(T.d))
      add(L, "k = " .. fmt(kk, 6) .. "   " .. U.arrow .. " " .. st)
    end
    add(L, "Verifique luego con k en el campo k.")
    return L
  end
  local Gc = tf(v, "Gc", k, "1")
  local Gp = tf(v, "Gp", k)
  local G0 = rsimplify(rmul(Gc, Gp))
  sep(L, "Lazo abierto")
  add(L, "G0 = " .. rstr(G0))
  local T = closedloop(G0)
  local st, ok = stabilitytext(proots(T.d))
  add(L, "Lazo cerrado: " .. st)
  if not ok then add(L, "OJO: teorema valor final NO aplica (inestable)") end
  sep(L, "Constantes y errores (A = " .. fmt(A) .. ")")
  local e = sserrors(G0, A)
  errlines(L, e, A)
  add(L, "Parabola: r(t)=A t" .. U.sq .. "/2 (R=A/s" .. "\194\179" .. ")")
  return L
end

--------------------------------------------------------------------------
-- 4) Routh-Hurwitz
--------------------------------------------------------------------------
local function charpoly(v, kk)
  if not blank(v.P) then
    local r = tf(v, "P", kk)
    return r.n
  end
  local G0 = rmul(tf(v, "Gc", kk, "1"), tf(v, "Gp", kk))
  return padd(G0.d, G0.n)
end

local function m_routh(v)
  local L = {}
  local k = num(v, "k")
  local hask = usesk(v.P) or (blank(v.P) and (usesk(v.Gc) or usesk(v.Gp)))
  if hask and k == nil then
    local P0, P1, P2 = charpoly(v, 0), charpoly(v, 1), charpoly(v, 2)
    local d1 = psub(P1, P0)
    local lin = true
    local chk = psub(psub(P2, P0), pscale(d1, 2))
    for i = 1, #chk do if abs(chk[i]) > 1e-9 * math.max(1, abs(P2[i] or 0)) then lin = false end end
    sep(L, "Polinomio caracteristico en funcion de k")
    if lin then
      add(L, "P(s) = " .. pstr(P0) .. "  +  k(" .. pstr(d1) .. ")")
    else
      add(L, "P(s) depende de k en forma no lineal")
      add(L, "P(s)|k=1 = " .. pstr(P1))
    end
    sep(L, "Rango de k estable (Routh)")
    add(L, "Se exige 1a columna de Routh > 0; los limites son los k donde")
    add(L, "un elemento se anula (polos sobre el eje j" .. W .. ").")
    local ivs, bounds = krange(function(kk) return charpoly(v, kk) end)
    if #ivs == 0 then add(L, "Ningun k en [-1e5,1e5] estabiliza el lazo") end
    for _, iv in ipairs(ivs) do
      local a = (iv[1] == -INF) and ("-" .. U.inf) or fmt(iv[1], 6)
      local b = (iv[2] == INF) and U.inf or fmt(iv[2], 6)
      add(L, a .. " < k < " .. b)
      if iv[1] < 0 and iv[2] > 0 then
        add(L, "  (fisico k>0: 0 < k < " .. b .. ")")
      end
    end
    sep(L, "Limites (estabilidad marginal)")
    for _, kb in ipairs(bounds) do
      local r = proots(charpoly(v, kb))
      local cross = {}
      for _, q in ipairs(r) do
        if abs(q[1]) < 1e-4 * math.max(1, abs(q[2])) and q[2] >= 0 then cross[#cross + 1] = q end
      end
      local txt = "k = " .. fmt(kb, 6)
      if #cross > 0 then
        local w = cross[1][2]
        if w > 1e-6 then
          txt = txt .. ":  s = " .. PM .. "j" .. fmt(w) .. "  (" .. W .. " = " .. fmt(w) .. " rad/s)"
        else
          txt = txt .. ":  raiz en s = 0"
        end
      end
      add(L, txt)
    end
    add(L, "(ingrese k para ver la tabla de Routh)")
    return L
  end
  local P = charpoly(v, k)
  sep(L, "Polinomio caracteristico")
  add(L, "P(s) = " .. pstr(P))
  local rows, notes, ch, _, _, epscells = routh(P)
  sep(L, "Tabla de Routh")
  local n = #P - 1
  for i, r in ipairs(rows) do
    local cells = {}
    for j = 1, #r do
      if epscells and epscells[i] and j == 1 then cells[#cells + 1] = U.eps
      elseif j == 1 or abs(r[j]) > 0 or j < #r then cells[#cells + 1] = fmt(r[j], 5) end
    end
    while #cells > 1 and cells[#cells] == "0" do cells[#cells] = nil end
    add(L, "s^" .. (n - i + 1) .. " | " .. table.concat(cells, "   "))
  end
  for _, s in ipairs(notes) do add(L, s) end
  local fc = {}
  for i = 1, #rows do
    if epscells and epscells[i] then fc[#fc + 1] = U.eps .. "(+)"
    else fc[#fc + 1] = fmt(rows[i][1], 4) end
  end
  add(L, "1a columna: " .. table.concat(fc, ", "))
  add(L, "Cambios de signo: " .. ch .. "  " .. U.arrow .. " " .. ch .. " raiz(es) en SPD")
  add(L, "(estable si toda la 1a columna tiene el mismo signo)")
  local roots = proots(P)
  add(L, "Raices: " .. rootsstr(roots))
  add(L, (stabilitytext(roots)))
  return L
end

--------------------------------------------------------------------------
-- 5) Respuesta en frecuencia y margenes
--------------------------------------------------------------------------
-- expresiones de |G(jw)| y angulo G(jw) a partir de polos y ceros
local function freqformula(o)
  local num, den, phz, php = {}, {}, {}, {}
  local n0z, n0p = 0, 0
  local function fac(r, isz)
    local m, a = isz and num or den, isz and phz or php
    if abs(r[1]) < 1e-11 and r[2] == 0 then
      if isz then n0z = n0z + 1 else n0p = n0p + 1 end
    elseif r[2] == 0 then
      local c = fmt(abs(r[1]), 4)
      m[#m + 1] = U.sqrt .. "(" .. W .. U.sq .. "+" .. c .. U.sq .. ")"
      if r[1] < 0 then a[#a + 1] = "atan(" .. W .. "/" .. c .. ")"
      else a[#a + 1] = "(180" .. DEG .. "-atan(" .. W .. "/" .. c .. "))" end
    elseif r[2] > 0 then
      local b, c = fmt(r[1] * r[1] + r[2] * r[2], 4), fmt(-2 * r[1], 4)
      m[#m + 1] = U.sqrt .. "((" .. b .. "-" .. W .. U.sq .. ")" .. U.sq .. "+(" .. c .. W .. ")" .. U.sq .. ")"
      a[#a + 1] = "atan2(" .. c .. W .. ", " .. b .. "-" .. W .. U.sq .. ")"
    end
  end
  for _, r in ipairs(o.z) do fac(r, true) end
  for _, r in ipairs(o.p) do fac(r, false) end
  local K = fmt(abs(o.lead))
  local ns = K .. ((#num > 0) and ("\194\183" .. table.concat(num, "\194\183")) or "")
  local dpart = {}
  if n0p - n0z > 0 then dpart[#dpart + 1] = W .. ((n0p - n0z > 1) and ("^" .. (n0p - n0z)) or "") end
  for _, x in ipairs(den) do dpart[#dpart + 1] = x end
  local ms = "|G0(j" .. W .. ")| = " .. ns .. ((#dpart > 0) and (" / [" .. table.concat(dpart, "\194\183") .. "]") or "")
  local ph = {}
  local nint = n0p - n0z
  if nint ~= 0 then ph[#ph + 1] = fmt(-90 * nint) .. DEG end
  if o.lead < 0 then ph[#ph + 1] = "-180" .. DEG end
  for _, x in ipairs(phz) do ph[#ph + 1] = "+ " .. x end
  for _, x in ipairs(php) do ph[#ph + 1] = "- " .. x end
  if o.L ~= 0 then ph[#ph + 1] = "- " .. fmt(o.L) .. W .. "(180/" .. U.pi .. ")" end
  local ps = U.ang .. "G0(j" .. W .. ") = " .. ((#ph > 0) and table.concat(ph, " "):gsub("^%+ ", "") or "0")
  return ms, ps
end

-- margenes explicados paso a paso
local function marginlines(L, o, label, maxpc, quiet)
  local m = margins(o)
  if label then add(L, label) end
  if not quiet then
    local ms, ps = freqformula(o)
    add(L, "  " .. ms)
    add(L, "  " .. ps)
  end
  -- paso 1: margen de fase
  if #m.gc == 0 then
    add(L, "  a) |G0|=1 no tiene solucion " .. U.arrow .. " MF = " .. U.inf)
  end
  for i, c in ipairs(m.gc) do
    add(L, "  a) MF: resolver |G0(j" .. W .. ")| = 1 " .. U.arrow .. " " .. W .. "cg = " .. fmt(c.w) .. " rad/s")
    add(L, "     " .. U.ang .. "G0(j" .. W .. "cg) = " .. fmt(c.ph) .. DEG .. "  " .. U.arrow .. "  MF = 180" .. DEG .. " + (" .. fmt(c.ph) .. DEG .. ") = " .. fmt(c.MF) .. DEG)
  end
  -- paso 2: margen de ganancia
  if #m.pc == 0 then
    add(L, "  b) " .. U.ang .. "G0 nunca llega a -180" .. DEG .. " " .. U.arrow .. " MG = " .. U.inf)
  end
  for i, c in ipairs(m.pc) do
    if i <= (maxpc or 1) then
      local tgt = fmt(c.ph, 4)
      add(L, "  b) MG: resolver " .. U.ang .. "G0(j" .. W .. ") = " .. tgt .. DEG .. " " .. U.arrow .. " " .. W .. "cf = " .. fmt(c.w) .. " rad/s")
      add(L, "     |G0(j" .. W .. "cf)| = " .. fmt(c.mag) .. "  " .. U.arrow .. "  MG = 1/" .. fmt(c.mag) .. " = " .. fmt(c.MG))
      add(L, "     en dB: 20" .. "\194\183" .. "log10(" .. fmt(c.MG) .. ") = " .. fmt(db(c.MG), 4) .. " dB")
    end
  end
  local mf = m.gc[1] and m.gc[1].MF or INF
  local mg = m.pc[1] and m.pc[1].MG or INF
  if mf > 0 and mg > 1 then
    add(L, "  c) MF > 0" .. DEG .. " y MG > 1 (" .. U.ge .. "0 dB) " .. U.arrow .. " lazo cerrado ESTABLE")
  else
    add(L, "  c) MF " .. U.le .. " 0" .. DEG .. " o MG " .. U.le .. " 1 " .. U.arrow .. " lazo cerrado INESTABLE")
  end
  return m
end

local function m_freq(v)
  local L = {}
  local k = num(v, "k")
  local Ld = num(v, "L") or 0
  local Gc = tf(v, "Gc", k, "1")
  local Gp = tf(v, "Gp", k)
  local G0 = rsimplify(rmul(Gc, Gp))
  sep(L, "Lazo abierto")
  add(L, "G0 = " .. rstr(G0) .. (Ld ~= 0 and ("  " .. "\194\183" .. " e^(-" .. fmt(Ld) .. "s)") or ""))
  local o = freqobj(G0, Ld)
  sep(L, "Margenes de estabilidad")
  add(L, "MF: fase extra que tolera el lazo antes de inestabilizarse.")
  add(L, "MG: cuantas veces puede crecer la ganancia antes de inestabilizarse.")
  local m = marginlines(L, o, nil, 2)
  local mf = m.gc[1] and m.gc[1].MF or INF
  local mg = m.pc[1] and m.pc[1].MG or INF
  if m.pc[1] then
    add(L, "Ganancia critica: Kcrit = MG" .. "\194\183" .. "(ganancia actual) " .. U.arrow .. " factor " .. fmt(mg))
    if blank(v.Gc) or v.Gc:match("^%s*1%s*$") then add(L, "ZN (Gc=1): Ku = " .. fmt(mg) .. ",  Pu = 2" .. U.pi .. "/" .. W .. "cf = " .. fmt(2 * pi / m.pc[1].w) .. " s") end
  end
  sep(L, "Estabilidad lazo cerrado")
  if Ld == 0 then
    local T = closedloop(G0)
    local r = proots(T.d)
    add(L, "Polos LC: " .. rootsstr(r, 4))
    add(L, (stabilitytext(r)))
  else
    if mf > 0 and mg > 1 then add(L, "MF > 0 y MG > 1 " .. U.arrow .. " ESTABLE (criterio de margenes)")
    else add(L, "MF " .. U.le .. " 0 o MG " .. U.le .. " 1 " .. U.arrow .. " INESTABLE (criterio de margenes)") end
  end
  local w = num(v, "w")
  if w then
    local mag, ph = freqat(o, w)
    sep(L, "G0(j" .. W .. ") en " .. W .. " = " .. fmt(w))
    add(L, "|G0| = " .. fmt(mag) .. "  (" .. fmt(db(mag), 4) .. " dB)   " .. U.ang .. "G0 = " .. fmt(ph) .. DEG)
    local A = num(v, "A")
    if A then
      add(L, "u = " .. fmt(A) .. " sin(" .. fmt(w) .. "t)  " .. U.arrow)
      add(L, "yss = " .. fmt(A * mag) .. " sin(" .. fmt(w) .. "t " .. (ph < 0 and "- " or "+ ") .. fmt(abs(ph)) .. DEG .. ")")
    end
  end
  return L
end

--------------------------------------------------------------------------
-- 6) Diseno por localizacion de polos
--------------------------------------------------------------------------
local function m_place(v)
  local L = {}
  local Gp = rsimplify(tf(v, "Gp", nil))
  local sig, wd = num(v, "sigma"), num(v, "wd")
  local crit = num(v, "crit") or 2
  local ts, Mp = num(v, "ts"), num(v, "Mp")
  local enulo = num(v, "enulo") or 0
  local estr = num(v, "estr") or 0
  local p3 = num(v, "p3")
  sep(L, "Paso 1: Polos deseados pd")
  add(L, "(polos de LC que producen el ts y Mp pedidos)")
  local z, wn
  if sig and wd then
    wn = sqrt(sig * sig + wd * wd); z = sig / wn
    add(L, "Dados: " .. SIG .. " = " .. fmt(sig) .. ", " .. W .. "d = " .. fmt(wd) .. "  (" .. Z .. " = " .. fmt(z) .. ", " .. WN .. " = " .. fmt(wn) .. ")")
  else
    req(ts, "ts"); req(Mp, "Mp")
    local l = log(Mp / 100)
    z = -l / sqrt(pi * pi + l * l)
    local c = (crit == 5) and 3 or 4
    sig = c / ts
    wn = sig / z
    wd = wn * sqrt(1 - z * z)
    add(L, Z .. " = -ln(Mp/100)/" .. U.sqrt .. "(" .. U.pi .. U.sq .. "+ln" .. U.sq .. ") = " .. fmt(z))
    add(L, SIG .. " = " .. Z .. WN .. " = " .. c .. "/ts = " .. fmt(sig) .. "   " .. WN .. " = " .. SIG .. "/" .. Z .. " = " .. fmt(wn))
    add(L, W .. "d = " .. WN .. U.sqrt .. "(1-" .. Z .. U.sq .. ") = " .. fmt(wd))
  end
  add(L, "pd = -" .. fmt(sig) .. " " .. PM .. " j" .. fmt(wd))
  add(L, "Pol. deseado: s" .. U.sq .. " + " .. fmt(2 * sig) .. "s + " .. fmt(sig * sig + wd * wd))

  local zp, pp, lead = proots(Gp.n), proots(Gp.d), plead(Gp.n) / plead(Gp.d)
  local Np = pzeroroots(Gp.d) - pzeroroots(Gp.n)
  sep(L, "Paso 2: Planta")
  add(L, "Gp = " .. rstr(Gp))
  add(L, "Polos: " .. rootsstr(pp, 4) .. "   Ceros: " .. rootsstr(zp, 4))
  add(L, "Tipo de la planta: " .. Np)

  -- PID algebraico con tercer polo fijo
  if p3 then
    sep(L, "PID algebraico (tercer polo en -" .. fmt(p3) .. ")")
    if #Gp.n ~= 1 or #Gp.d ~= 3 then error("PID algebraico: requiere Gp = b/(s^2+a1 s+a0)") end
    local b = Gp.n[1]; local a1, a0 = Gp.d[2], Gp.d[1]
    local des = pmul({sig * sig + wd * wd, 2 * sig, 1}, {p3, 1})
    add(L, "s(s" .. U.sq .. "+a1s+a0)+b(Kd s" .. U.sq .. "+Kp s+Ki) = " .. pstr(des))
    local Kd, Kp, Ki = (des[3] - a1) / b, (des[2] - a0) / b, des[1] / b
    add(L, "Kd = (" .. fmt(des[3]) .. "-" .. fmt(a1) .. ")/" .. fmt(b) .. " = " .. fmt(Kd))
    add(L, "Kp = (" .. fmt(des[2]) .. "-" .. fmt(a0) .. ")/" .. fmt(b) .. " = " .. fmt(Kp))
    add(L, "Ki = " .. fmt(des[1]) .. "/" .. fmt(b) .. " = " .. fmt(Ki))
    add(L, "Gc = " .. fmt(Kd) .. "s + " .. fmt(Kp) .. " + " .. fmt(Ki) .. "/s")
    add(L, "Ideal: Kp = " .. fmt(Kp) .. "  Ti = Kp/Ki = " .. fmt(Kp / Ki) .. "  Td = Kd/Kp = " .. fmt(Kd / Kp))
    local Gc = rat({Ki, Kp, Kd}, {0, 1})
    local G0 = rmul(Gc, Gp)
    local r = proots(padd(G0.d, G0.n))
    add(L, "Polos LC: " .. rootsstr(r, 4))
    return L
  end

  -- eleccion de integrador
  sep(L, "Paso 3: Integrador? (error estacionario)")
  local useI
  if estr == 3 or estr == 4 then useI = true
  elseif estr == 1 or estr == 2 then useI = false
  else
    if enulo == 1 then
      if Np >= 1 then
        useI = false
        add(L, "Se exige error nulo, pero la planta ya es tipo " .. Np .. " " .. U.arrow .. " no se agrega integrador")
      else
        useI = true
        add(L, "Se exige error nulo y planta tipo 0 " .. U.arrow .. " se agrega integrador 1/s")
      end
    else
      useI = false
      add(L, "No se exige error nulo " .. U.arrow .. " sin integrador (se prueba P)")
    end
  end
  if useI then add(L, "(el polo en s=0 del PI/PID hace el sistema tipo 1: ep = 0)") end
  -- angulos de cada polo y cero hacia pd
  sep(L, "Paso 4: Condicion de angulo en pd")
  add(L, "Para que pd sea polo de LC: " .. U.ang .. "Gc(pd)Gp(pd) = -180" .. DEG)
  add(L, U.ang .. "G(pd) = " .. "\206\163" .. U.theta .. "ceros - " .. "\206\163" .. U.theta .. "polos,  " .. U.theta .. " = angulo del vector (pd - raiz)")
  local mag, ph = gainphase_at(zp, pp, lead, -sig, wd)
  local sump, sumz = 0, 0
  local dp, dz_ = {}, {}
  local function th(r) return atan2(wd - r[2], -sig - r[1]) * R2D end
  local function dist(r) return sqrt((-sig - r[1]) ^ 2 + (wd - r[2]) ^ 2) end
  local np = 0
  for _, r in ipairs(pp) do
    np = np + 1
    local t = th(r)
    sump = sump + t
    dp[#dp + 1] = dist(r)
    add(L, "  " .. U.theta .. "p" .. np .. " (polo " .. fmtc(r[1], r[2], 4) .. ") = atan2(" .. fmt(wd - r[2], 4) .. ", " .. fmt(-sig - r[1], 4) .. ") = " .. fmt(t) .. DEG)
  end
  if useI then
    np = np + 1
    local t = atan2(wd, -sig) * R2D
    sump = sump + t
    dp[#dp + 1] = sqrt(sig * sig + wd * wd)
    add(L, "  " .. U.theta .. "p" .. np .. " (integrador, polo 0) = " .. fmt(t) .. DEG)
    ph = ph - t
    mag = mag / sqrt(sig * sig + wd * wd)
  end
  local nzr = 0
  for _, r in ipairs(zp) do
    nzr = nzr + 1
    local t = th(r)
    sumz = sumz + t
    dz_[#dz_ + 1] = dist(r)
    add(L, "  " .. U.theta .. "z" .. nzr .. " (cero " .. fmtc(r[1], r[2], 4) .. ") = " .. fmt(t) .. DEG)
  end
  add(L, U.ang .. "G(pd) = " .. fmt(sumz) .. DEG .. " - " .. fmt(sump) .. DEG .. ((lead < 0) and " - 180" .. DEG .. " (K<0)" or "") .. " = " .. fmt(ph) .. DEG)
  local phi = -180 - ph
  while phi > 180 do phi = phi - 360 end
  while phi <= -180 do phi = phi + 360 end
  add(L, "Deficiencia: " .. PHI .. "d = -180" .. DEG .. " - (" .. fmt(ph) .. DEG .. ") = " .. fmt(phi) .. DEG)
  add(L, "(" .. PHI .. "d = fase que deben aportar los ceros del controlador)")
  sep(L, "Paso 5: Estructura minima")
  local nz
  if estr == 1 then nz = 0
  elseif estr == 2 or estr == 3 then nz = 1
  elseif estr == 4 then nz = 2
  else
    if abs(phi) < 0.5 then nz = 0; add(L, "|" .. PHI .. "d| " .. U.approx .. " 0 " .. U.arrow .. " la condicion de angulo ya se cumple: sin ceros")
    elseif phi > 0 and phi < 90 then nz = 1; add(L, "0" .. DEG .. " < " .. PHI .. "d < 90" .. DEG .. " " .. U.arrow .. " basta UN cero (un cero aporta de 0" .. DEG .. " a 90" .. DEG .. ")")
    elseif phi >= 90 and phi < 180 then nz = 2; add(L, PHI .. "d " .. U.ge .. " 90" .. DEG .. " " .. U.arrow .. " un cero no basta: DOS ceros iguales (" .. PHI .. "d/2 c/u)")
    else
      add(L, PHI .. "d < 0: sobra fase; un cero no sirve.")
      add(L, "Pruebe otros polos o estructura forzada.")
      return L
    end
  end
  local names
  if useI then names = {[0] = "I  (Kc/s)", "PI", "PID"} else names = {[0] = "P", "PD", "PD doble (impropio)"} end
  add(L, "Estructura minima: " .. names[nz])
  if nz == 0 and abs(phi) >= 0.5 then
    add(L, "OJO: con estructura forzada sin ceros NO se cumple el angulo;")
    add(L, "los polos LC no pasaran por pd.")
  end

  local zc = nil
  local cmag = 1
  if nz > 0 then
    sep(L, "Paso 6: Ubicacion del cero en s = -a")
    local phz = phi / nz
    zc = sig + wd / math.tan(phz * D2R)
    add(L, "El cero debe 've' a pd con angulo " .. PHI .. " = " .. fmt(phz) .. DEG .. ":")
    add(L, "tan(" .. PHI .. ") = Im pd/(a - |Re pd|)  " .. U.arrow .. "  a = |Re pd| + Im pd/tan(" .. PHI .. ")")
    add(L, "a = " .. fmt(sig) .. " + " .. fmt(wd) .. "/tan(" .. fmt(phz) .. DEG .. ")")
    add(L, "   a = " .. fmt(zc) .. "   (cero en s = " .. fmt(-zc) .. ")")
    if zc < 0 then add(L, "OJO: cero en SPD") end
    local dz = sqrt((zc - sig) ^ 2 + wd * wd)
    cmag = dz ^ nz
    add(L, "|pd + a| = " .. fmt(dz) .. (nz == 2 and ("  " .. U.arrow .. " |pd+a|" .. U.sq .. " = " .. fmt(cmag)) or ""))
  end
  local Kc = 1 / (mag * cmag)
  sep(L, "Paso 7: Ganancia Kc (condicion de magnitud)")
  add(L, "|Gc(pd)Gp(pd)| = 1,  |s - r| = distancia de pd a cada raiz")
  local dps = {}
  for _, d in ipairs(dp) do dps[#dps + 1] = fmt(d, 4) end
  local nums = {}
  for _, d in ipairs(dz_) do nums[#nums + 1] = fmt(d, 4) end
  if nz > 0 then nums[#nums + 1] = "|pd+a|" .. (nz == 2 and U.sq or "") .. "=" .. fmt(cmag, 4) end
  add(L, "Kc = (dist. a polos) / (K" .. "\194\183" .. "dist. a ceros)")
  local nn = {fmt(abs(lead), 4)}
  for _, d in ipairs(dz_) do nn[#nn + 1] = fmt(d, 4) end
  if nz > 0 then nn[#nn + 1] = fmt(cmag, 4) .. ((nz == 2) and " [=|pd+a|" .. U.sq .. "]" or " [=|pd+a|]") end
  add(L, "Kc = (" .. table.concat(dps, "\194\183") .. ") / (" .. table.concat(nn, "\194\183") .. ") = " .. fmt(Kc))
  -- controlador
  local cn = {Kc}
  if nz >= 1 then cn = pmul(cn, {zc, 1}) end
  if nz >= 2 then cn = pmul(cn, {zc, 1}) end
  local Gc = rat(cn, useI and {0, 1} or {1})
  sep(L, "Paso 8: Controlador")
  local fact = fmt(Kc)
  if nz == 1 then fact = fact .. "(s + " .. fmt(zc) .. ")" elseif nz == 2 then fact = fact .. "(s + " .. fmt(zc) .. ")" .. U.sq end
  if useI then fact = fact .. "/s" end
  add(L, "Gc(s) = " .. fact)
  if useI and nz == 0 then
    add(L, "Ki = " .. fmt(Kc))
  elseif useI and nz == 1 then
    add(L, "Ideal: Kp = Kc = " .. fmt(Kc) .. "   Ti = 1/a = " .. fmt(1 / zc) .. " s")
    add(L, "Paralelo: Kp = " .. fmt(Kc) .. "   Ki = Kc" .. "\194\183" .. "a = " .. fmt(Kc * zc))
  elseif useI and nz == 2 then
    local Kp, Ki, Kd = Kc * 2 * zc, Kc * zc * zc, Kc
    add(L, "Ideal: Kp = 2Kc" .. "\194\183" .. "a = " .. fmt(Kp) .. "  Ti = 2/a = " .. fmt(2 / zc) .. "  Td = 1/(2a) = " .. fmt(1 / (2 * zc)))
    add(L, "Paralelo: Kp = " .. fmt(Kp) .. "  Ki = " .. fmt(Ki) .. "  Kd = " .. fmt(Kd))
  elseif not useI and nz == 1 then
    add(L, "Ideal Kp(1+Td s): Kp = Kc" .. "\194\183" .. "a = " .. fmt(Kc * zc) .. "   Td = 1/a = " .. fmt(1 / zc) .. " s")
    add(L, "Paralelo: Kp = " .. fmt(Kc * zc) .. "   Kd = " .. fmt(Kc))
  else
    add(L, "Kp = " .. fmt(Kc))
  end
  -- verificacion
  sep(L, "Paso 9: Verificacion")
  add(L, "Polos LC = raices de 1 + Gc" .. "\194\183" .. "Gp = 0; deben incluir pd.")
  local G0 = rmul(Gc, Gp)
  local cp = padd(G0.d, G0.n)
  add(L, "1+GcGp=0: " .. pstr(pscale(cp, 1 / plead(cp))) .. " = 0")
  local r = proots(cp)
  add(L, "Polos LC: " .. rootsstr(r, 5))
  local best = INF
  for _, q in ipairs(r) do
    local d = abs(q[1] + sig) + abs(abs(q[2]) - wd)
    if d < best then best = d end
  end
  if best < 1e-3 * math.max(1, wd) then add(L, "pd pertenece a los polos LC  " .. U.arrow .. " OK")
  else add(L, "pd NO es polo LC (revise estructura)") end
  local others = {}
  for _, q in ipairs(r) do
    if not (abs(q[1] + sig) < 1e-3 * math.max(1, sig) and abs(abs(q[2]) - wd) < 1e-3 * math.max(1, wd)) then others[#others + 1] = q end
  end
  for _, q in ipairs(others) do
    add(L, "Polo adicional " .. fmtc(q[1], q[2], 4) .. ": |Re|/" .. SIG .. " = " .. fmt(abs(q[1]) / sig, 3) .. ((abs(q[1]) / sig >= 5) and " (" .. U.ge .. "5, despreciable)" or " (<5, influye)"))
  end
  add(L, (stabilitytext(r)))
  sep(L, "Paso 10: Error estacionario (A=1)")
  local e = sserrors(rsimplify(G0), 1)
  errlines(L, e, 1)
  if e.ep > 0 and e.ep < INF then
    add(L, "y(" .. U.inf .. ") = 1 - ep = " .. fmt(1 - e.ep) .. " (escalon unitario)")
  end
  return L
end

--------------------------------------------------------------------------
-- 7) Sintonia empirica (ZN, rele, Ho-Hang-Cao)
--------------------------------------------------------------------------
local function pidtf(Kp, Ti, Td)
  -- Kp(1 + 1/(Ti s) + Td s)
  if not Ti then
    if Td then return rat({Kp, Kp * Td}, {1}) end
    return rconst(Kp)
  end
  if Td then return rat({Kp, Kp * Ti, Kp * Ti * Td}, {0, Ti}) end
  return rat({Kp, Kp * Ti}, {0, Ti})
end

local function m_empir(v)
  local L = {}
  local K, tau, Ld, tau2 = num(v, "K"), num(v, "tau"), num(v, "L"), num(v, "tau2")
  local MG, MF = num(v, "MG"), num(v, "MF")
  local Ku, Pu, d, a = num(v, "Ku"), num(v, "Pu"), num(v, "d"), num(v, "a")
  local fK = num(v, "fKp")
  local haveModel = K and tau and Ld
  if haveModel then
    if K == 0 then error("K no puede ser 0") end
    if tau <= 0 then error(TAU .. " debe ser > 0") end
    if Ld <= 0 then error("L (retardo) debe ser > 0 para ZN/HHC") end
    if tau2 and tau2 < 0 then error(TAU .. "2 debe ser >= 0") end
  end
  if MG and MG <= 1 then error("MG debe ser > 1 (escala lineal)") end
  local Gp
  if haveModel then
    local den = {1, tau}
    if tau2 then den = pmul(den, {1, tau2}) end
    Gp = rat({K}, den)
    sep(L, "Modelo")
    add(L, "Gp = " .. fmt(K) .. "e^(-" .. fmt(Ld) .. "s)/" .. (tau2 and "[" or "") .. "(" .. pstr({1, tau}) .. ")" .. (tau2 and ("(" .. pstr({1, tau2}) .. ")]") or ""))
    add(L, "Retardo normalizado L/" .. TAU .. " = " .. fmt(Ld / tau, 4))
  end
  local verify = {}
  -- ZN lazo abierto
  if haveModel then
    sep(L, "ZN curva de reaccion (lazo abierto)")
    local P = tau / (K * Ld)
    add(L, "P  : Kp = " .. TAU .. "/(KL) = " .. fmt(P))
    add(L, "PI : Kp = 0.9" .. TAU .. "/(KL) = " .. fmt(0.9 * P) .. "   Ti = 3.33L = " .. fmt(3.33 * Ld))
    add(L, "PID: Kp = 1.2" .. TAU .. "/(KL) = " .. fmt(1.2 * P) .. "  Ti = 2L = " .. fmt(2 * Ld) .. "  Td = 0.5L = " .. fmt(0.5 * Ld))
    verify[#verify + 1] = {"ZN-LA PI", 0.9 * P, 3.33 * Ld, nil}
    verify[#verify + 1] = {"ZN-LA PID", 1.2 * P, 2 * Ld, 0.5 * Ld}
  end
  -- ganancia ultima
  local src
  if Ku then src = "dato"
  elseif d and a then
    Ku = 4 * d / (pi * a); src = "rele"
  elseif haveModel then
    local o = freqobj(Gp, Ld)
    local m = margins(o)
    if m.pc[1] then
      Ku = m.pc[1].MG; Pu = Pu or 2 * pi / m.pc[1].w; src = "modelo"
    end
  end
  if Ku then
    sep(L, "ZN oscilacion sostenida (lazo cerrado)")
    if src == "rele" then add(L, "Rele: Ku " .. U.approx .. " 4d/(" .. U.pi .. "a) = 4" .. "\194\183" .. fmt(d) .. "/(" .. U.pi .. "" .. "\194\183" .. fmt(a) .. ") = " .. fmt(Ku)) end
    if src == "modelo" then add(L, "Del modelo (cruce de fase): Ku = " .. fmt(Ku) .. ", Pu = " .. fmt(Pu) .. " s") end
    add(L, "P  : Kp = 0.5Ku = " .. fmt(0.5 * Ku))
    if Pu then
      add(L, "PI : Kp = 0.45Ku = " .. fmt(0.45 * Ku) .. "   Ti = Pu/1.2 = " .. fmt(Pu / 1.2))
      add(L, "PID: Kp = 0.6Ku = " .. fmt(0.6 * Ku) .. "  Ti = Pu/2 = " .. fmt(Pu / 2) .. "  Td = Pu/8 = " .. fmt(Pu / 8))
      verify[#verify + 1] = {"ZN-LC PI", 0.45 * Ku, Pu / 1.2, nil}
      verify[#verify + 1] = {"ZN-LC PID", 0.6 * Ku, Pu / 2, Pu / 8}
    else
      add(L, "PI : Kp = 0.45Ku = " .. fmt(0.45 * Ku) .. "  (falta Pu para Ti)")
    end
  end
  -- Ho-Hang-Cao
  local hhc
  if haveModel and MG and MF then
    sep(L, "Ho-Hang-Cao")
    local mf = MF * D2R
    add(L, "MF = " .. fmt(MF) .. DEG .. " = " .. fmt(mf) .. " rad,  MG = " .. fmt(MG) .. " (lineal)")
    local wp = (MG * (mf - pi / 2) + (pi / 2) * MG * MG) / ((MG * MG - 1) * Ld)
    add(L, W .. "p = [MG(MF-" .. U.pi .. "/2)+(" .. U.pi .. "/2)MG" .. U.sq .. "]/[(MG" .. U.sq .. "-1)td]")
    add(L, "   = [" .. fmt(MG * (mf - pi / 2)) .. " + " .. fmt((pi / 2) * MG * MG) .. "]/" .. fmt((MG * MG - 1) * Ld) .. " = " .. fmt(wp) .. " rad/s")
    local Kp = wp * tau / (MG * K)
    add(L, "Kp = " .. W .. "p" .. TAU .. "/(MG K) = " .. fmt(wp * tau) .. "/" .. fmt(MG * K) .. " = " .. fmt(Kp))
    local t1, t2, t3 = 2 * wp, 4 * wp * wp * Ld / pi, 1 / tau
    local Ti = 1 / (t1 - t2 + t3)
    add(L, "Ti = 1/(2" .. W .. "p - 4" .. W .. "p" .. U.sq .. "td/" .. U.pi .. " + 1/" .. TAU .. ")")
    add(L, "   = 1/(" .. fmt(t1) .. " - " .. fmt(t2) .. " + " .. fmt(t3) .. ") = " .. fmt(Ti) .. " s")
    if Ti <= 0 then add(L, "OJO: Ti <= 0 (especificaciones no factibles)") end
    add(L, W .. "g = Kp K/" .. TAU .. " = " .. fmt(Kp * K / tau) .. " rad/s")
    add(L, "Ki = Kp/Ti = " .. fmt(Kp / Ti))
    if tau2 then add(L, "PID: Td = " .. TAU .. "2 = " .. fmt(tau2) .. " s") end
    add(L, "Gc = " .. fmt(Kp) .. "(1 + 1/(" .. fmt(Ti) .. "s)" .. (tau2 and (" + " .. fmt(tau2) .. "s") or "") .. ")")
    hhc = {"HHC " .. (tau2 and "PID" or "PI"), Kp, Ti, tau2}
    table.insert(verify, 1, hhc)
  end
  -- verificacion con retardo exacto
  if haveModel then
    sep(L, "Verificacion (retardo exacto)")
    add(L, "Con G0 = Gc" .. "\194\183" .. "Gp y el retardo e^(-L" .. W .. "): a) MF en " .. W .. "cg, b) MG en " .. W .. "cf, c) estabilidad.")
    marginlines(L, freqobj(Gp, Ld), "Planta sola (Gc = 1):", 1)
    local base
    for i, c in ipairs(verify) do
      local G0 = rmul(pidtf(c[2], c[3], c[4]), Gp)
      local m = marginlines(L, freqobj(G0, Ld), c[1] .. " (Kp=" .. fmt(c[2], 4) .. ", Ti=" .. fmt(c[3], 4) .. (c[4] and (", Td=" .. fmt(c[4], 4)) or "") .. "):", 1)
      if i == 1 and c == hhc then base = m end
    end
    if fK and hhc then
      local Kp2 = fK * hhc[2]
      sep(L, "Efecto de Kp x " .. fmt(fK) .. " (Ti fijo)")
      local G0 = rmul(pidtf(Kp2, hhc[3], hhc[4]), Gp)
      local m = marginlines(L, freqobj(G0, Ld), "Kp = " .. fmt(Kp2) .. ", Ti = " .. fmt(hhc[3]) .. ":", 1)
      local mf0 = base and base.gc[1] and base.gc[1].MF
      local mg0 = base and base.pc[1] and base.pc[1].MG
      local mf1 = m.gc[1] and m.gc[1].MF
      local mg1 = m.pc[1] and m.pc[1].MG
      if mf0 and mf1 then add(L, "MF: " .. fmt(mf0, 4) .. DEG .. " " .. U.arrow .. " " .. fmt(mf1, 4) .. DEG) end
      if mg0 and mg1 then add(L, "MG: " .. fmt(mg0, 4) .. " " .. U.arrow .. " " .. fmt(mg1, 4) .. "  (MG se divide por " .. fmt(fK) .. ")") end
      if mf0 and mf1 then
        if mf1 < mf0 then
          add(L, "Aumenta " .. W .. "cg (respuesta mas rapida) pero bajan MF y MG:")
          add(L, "menor estabilidad relativa, mas oscilacion/sobrepaso.")
        else
          add(L, "Los margenes aumentan: mayor estabilidad relativa.")
        end
        if (mf1 <= 0) or (mg1 and mg1 <= 1) then add(L, "MF<=0 o MG<=1 " .. U.arrow .. " lazo INESTABLE.") end
      end
    end
  end
  if #L == 0 then error("Ingrese K, " .. TAU .. ", L (modelo) y/o Ku, Pu, d, a") end
  return L
end

--------------------------------------------------------------------------
-- ayuda
--------------------------------------------------------------------------
local HELP = {
  "== Que modulo usar segun el enunciado ==",
  "'Mp y ts dados, hallar " .. U.zeta .. ", " .. U.wn .. ", polos' " .. U.arrow .. " 1",
  "'Determine y(t) ante escalon de magnitud A / modos' " .. U.arrow .. " 2 (G, A=magnitud, entrada=1, lazo=0)",
  "'ts y Mp de Y/R con K dado / polos dominantes / modelo reducido' " .. U.arrow .. " 3 (Gc=K, Gp)",
  "'Y/R en funcion de K' " .. U.arrow .. " 3 con Gc = k y campo k vacio",
  "'error estacionario / K para un error' " .. U.arrow .. " 4",
  "'rango de K estable / Routh' " .. U.arrow .. " 5 (Gc = k, k vacio)",
  "'MF, MG, frecuencias de cruce, |G(j" .. U.w .. ")|' " .. U.arrow .. " 6",
  "'Disene P/PD/PI/PID para ts y Mp' " .. U.arrow .. " 7 (error nulo=1 si pide respuesta exacta)",
  "'Ziegler-Nichols / rele / Ho-Hang-Cao' " .. U.arrow .. " 8 (K, " .. U.tau .. ", L del modelo)",
  "== Glosario ==",
  U.w .. "cg: cruce de GANANCIA (|G|=1, 0 dB). 'Frecuencia de corte'. Ahi se mide MF.",
  "MF = 180" .. U.deg .. " + " .. U.ang .. "G(j" .. U.w .. "cg). Estable si MF > 0.",
  U.w .. "cf: cruce de FASE (" .. U.ang .. "G = -180" .. U.deg .. "). Ahi se mide MG.",
  "MG = 1/|G(j" .. U.w .. "cf)| (lineal); en dB: 20log(MG). Estable si MG > 1.",
  "Kcrit = Ku: ganancia que lleva a oscilacion (= MG con Gc=1). Pu = Tcrit = 2" .. U.pi .. "/" .. U.w .. "cf.",
  "HHC: " .. U.w .. "p = " .. U.w .. "cf esperada,  " .. U.w .. "g = " .. U.w .. "cg esperada (diseno).",
  U.sigma .. " = |Re| del polo = " .. U.zeta .. U.wn .. ";  " .. U.w .. "d = Im del polo = " .. U.wn .. U.sqrt .. "(1-" .. U.zeta .. U.sq .. ").",
  U.zeta .. ": amortiguamiento;  " .. U.wn .. ": frecuencia natural;  pd: polo deseado.",
  U.phi .. "d: deficiencia angular (fase que debe aportar el controlador).",
  "SPI/SPD: semiplano izquierdo (estable) / derecho (inestable).",
  "Tipo N: n" .. "\194\186" .. " de integradores (polos en s=0) de G0 = Gc" .. "\194\183" .. "Gp.",
  "Kp, Kv, Ka: constantes de error (posicion, velocidad, aceleracion).",
  "ep, ev, ea: error estacionario ante escalon, rampa, parabola.",
  U.tau .. ": constante de tiempo; L = td: retardo; K: ganancia estatica.",
  "== Teclas ==",
  "Flechas arriba/abajo: cambiar campo",
  "Flechas izq/der: mover cursor   del: borrar",
  "clear: vaciar campo   enter: CALCULAR",
  "esc: volver   menu: elegir modulo",
  "En resultados: flechas para desplazar",
  "== Como escribir G(s) ==",
  "Use s como variable. Multiplicacion implicita:",
  "4/((s+1)(s+2))     2/(s^2+2s+17)",
  "(s+9)/s     k/(s(s+1)(s+10))     10s/(s^2+4s+8)",
  "kp(1+1/(ti s)) -> escriba numeros: 1.3(1+1/(5s))",
  "Funciones: sqrt ln log exp sin cos tan, pi, e",
  "El retardo e^(-Ls) va en el campo L (no en G).",
  "k = parametro libre: dejelo vacio para",
  "obtener el rango (Routh) o despejarlo (error).",
  "== Convenios ==",
  "ts(2%)=4/(" .. U.zeta .. U.wn .. "), ts(5%)=3/(" .. U.zeta .. U.wn .. "). MF en grados, MG lineal y dB.",
  "Retro unitaria negativa. Campos en blanco = desconocidos.",
  "Estructura (mod. 7): 0 auto, 1 P, 2 PD, 3 PI, 4 PID.",
  "error nulo: 1 = se exige seguimiento exacto (escalon).",
  "Senal (mod. 4): 1 escalon, 2 rampa, 3 parabola.",
}

--------------------------------------------------------------------------
-- definicion de modulos para la interfaz
--------------------------------------------------------------------------
MODULES = {
  {title = "Especificaciones 2o orden", short = "Mp, ts, " .. U.zeta .. ", " .. U.wn .. ", polos", run = m_specs, fields = {
    {"Mp", "Mp [%]", "Sobrepaso en %. Ej: 16.3"},
    {"ts", "ts [s]", "Tiempo de asentamiento"},
    {"crit", "criterio ts", "2 = 2% (4/" .. U.sigma .. ")   5 = 5% (3/" .. U.sigma .. ")", "2"},
    {"zeta", U.zeta, "Amortiguamiento (si lo conoce)"},
    {"wn", U.wn .. " [rad/s]", "Frecuencia natural (si la conoce)"},
    {"tp", "tp [s]", "Tiempo de peak (opcional)"},
    {"sigma", U.sigma .. " (=|Re p|)", "Parte real del polo (opcional)"},
    {"wd", U.w .. "d (=Im p)", "Parte imaginaria del polo (opcional)"},
  }},
  {title = "Respuesta y(t)", short = "fracciones parciales, modos", run = m_resp, fields = {
    {"G", "G(s)", "La G(s) del enunciado. Ej: 2/(s^3+4s^2+s+6)"},
    {"A", "magnitud A", "'escalon de magnitud 5' " .. U.arrow .. " escriba 5 (unitario = 1)", "1"},
    {"entrada", "tipo entrada", "1 = escalon (lo normal), 2 = rampa, 0 = impulso", "1"},
    {"lazo", "lazo", "0 = usar G tal cual (lo normal); 1 = cerrar lazo G/(1+G)", "0"},
    {"k", "k", "Solo si G contiene k; si no, dejar vacio"},
  }},
  {title = "Polos y respuesta (LA/LC)", short = "GLC, polos, ts, Mp, reduccion", run = m_polos, fields = {
    {"Gc", "Gc(s)", "Controlador. Ej: 9   ,  3(s+9.389)   (vacio = 1)", "1"},
    {"Gp", "Gp(s)", "Planta. Ej: 7/(s(s+10)(s+1))"},
    {"lazo", "lazo", "1 = cerrado retro unitaria,  0 = solo Gc" .. "\194\183" .. "Gp", "1"},
    {"A", "A escalon", "Amplitud del escalon", "1"},
    {"k", "k", "Valor del parametro k (si se usa)"},
  }},
  {title = "Error estado estacionario", short = "tipo, Kp Kv Ka, despeje", run = m_error, fields = {
    {"Gc", "Gc(s)", "Controlador (vacio = 1)", "1"},
    {"Gp", "Gp(s)", "Planta. Ej: 63/(s(s+10)(s+1)) ; K/(s+1)"},
    {"A", "A", "Amplitud de la referencia", "1"},
    {"k", "k", "Valor de k (vacio + e obj = despejar k)"},
    {"eobj", "e obj", "Error deseado (ej 0.1) para despejar k"},
    {"senal", "senal", "1 escalon  2 rampa  3 parabola (para despeje)", "1"},
  }},
  {title = "Routh-Hurwitz", short = "tabla, rango de k", run = m_routh, fields = {
    {"P", "P(s)", "Pol. caracteristico (opcional). Ej: s^3+4s^2+5s+k"},
    {"Gc", "Gc(s)", "Si P vacio: P = den+num de Gc" .. "\194\183" .. "Gp. Ej: k", "k"},
    {"Gp", "Gp(s)", "Planta. Ej: 1/(s(s^2+5s+1))"},
    {"k", "k", "Vacio = rango estable;  numero = tabla"},
  }},
  {title = "Frecuencia y margenes", short = "MF, MG, " .. U.w .. "cg, " .. U.w .. "cf, |G(j" .. U.w .. ")|", run = m_freq, fields = {
    {"Gc", "Gc(s)", "Controlador (vacio = 1). Ej: (s+9)/s", "1"},
    {"Gp", "Gp(s)", "Planta sin retardo. Ej: 20/((s+1)^2(s+6))"},
    {"L", "retardo L", "Retardo e^(-Ls) en s (0 si no hay)", "0"},
    {"w", U.w .. " eval.", "Frecuencia para evaluar G(j" .. U.w .. ") (opcional)"},
    {"A", "A seno", "Amplitud de u=A sin(" .. U.w .. "t) (opcional)"},
    {"k", "k", "Valor de k (si se usa)"},
  }},
  {title = "Localizacion de polos", short = "P, PD, PI, PID por angulo", run = m_place, fields = {
    {"Gp", "Gp(s)", "Planta. Ej: 4/((s+1)(s+2))"},
    {"ts", "ts [s]", "Tiempo de asentamiento"},
    {"crit", "criterio ts", "2 = 2%   5 = 5%", "2"},
    {"Mp", "Mp [%]", "Sobrepaso maximo en %"},
    {"enulo", "error nulo", "1 = se exige error nulo (seguim. exacto), 0 = no", "0"},
    {"estr", "estructura", "0 auto  1 P  2 PD  3 PI  4 PID", "0"},
    {"p3", "3er polo", "Solo PID algebraico Kd s+Kp+Ki/s (ej 5)"},
    {"sigma", U.sigma .. " (opc)", "Polo deseado directo: |Re pd| (opcional)"},
    {"wd", U.w .. "d (opc)", "Polo deseado directo: Im pd (opcional)"},
  }},
  {title = "Sintonia empirica PI/PID", short = "ZN, rele, Ho-Hang-Cao", run = m_empir, fields = {
    {"K", "K", "Ganancia del modelo Ke^(-Ls)/(" .. U.tau .. "s+1)"},
    {"tau", U.tau, "Constante de tiempo (" .. U.tau .. "1)"},
    {"L", "L = td", "Retardo"},
    {"tau2", U.tau .. "2", "Solo 2o polo (" .. U.tau .. "2 s+1): HHC PID Td=" .. U.tau .. "2"},
    {"MG", "MG", "Margen de ganancia LINEAL deseado (ej 3)"},
    {"MF", "MF [" .. U.deg .. "]", "Margen de fase deseado en grados (ej 60)"},
    {"Ku", "Ku", "Ganancia ultima (si es dato)"},
    {"Pu", "Pu", "Periodo ultimo (si es dato)"},
    {"d", "rele d", "Amplitud del rele " .. U.pm .. "d"},
    {"a", "rele a", "Amplitud de oscilacion de la salida"},
    {"fKp", "factor Kp", "Ej 2 = duplicar Kp (Ti fijo) y comparar"},
  }},
  {title = "Ayuda", short = "teclas y sintaxis", help = HELP, fields = {}},
}
