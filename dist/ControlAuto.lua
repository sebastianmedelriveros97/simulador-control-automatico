platform.apilevel = '2.0'
-- Control Automatico - asistente de examen (TI-Nspire CX II)
-- Generado por tools/build.py; editar los archivos en src/

-- ===== core.lua =====
--------------------------------------------------------------------------
-- core.lua : nucleo matematico (Lua 5.1, sin dependencias del Nspire)
-- polinomios, numeros complejos, parser de G(s), raices, Routh, frecuencia
--------------------------------------------------------------------------
local pi = math.pi
local abs, sqrt, atan2, log, exp, floor = math.abs, math.sqrt, math.atan2, math.log, math.exp, math.floor
local INF = 1/0
local D2R, R2D = pi/180, 180/pi

-- simbolos UTF-8 (escritos con escapes para mantener el fuente en ASCII)
U = {
  zeta="\206\182", wn="\207\137n", w="\207\137", phi="\207\134", tau="\207\132",
  sigma="\207\131", pi="\207\128", deg="\194\176", pm="\194\177", inf="\226\136\158",
  approx="\226\137\136", ang="\226\136\160", sq="\194\178", sqrt="\226\136\154",
  arrow="\226\134\146", theta="\206\184", le="\226\137\164", ge="\226\137\165", eps="\206\181",
}

--------------------------------------------------------------------------
-- formato numerico
--------------------------------------------------------------------------
function fmt(x, d)
  d = d or 5
  if x == nil then return "?" end
  if x ~= x then return "NaN" end
  if x == INF then return U.inf end
  if x == -INF then return "-" .. U.inf end
  if abs(x) < 1e-11 then return "0" end
  local s = string.format("%." .. d .. "g", x)
  -- quitar exponentes feos tipo e+00
  s = s:gsub("e%+0*(%d)", "e%1"):gsub("e%-0*(%d)", "e-%1")
  return s
end

function fmtc(re, im, d)
  if abs(im) < 1e-11 then return fmt(re, d) end
  if abs(re) < 1e-11 then
    return (im < 0 and "-" or "") .. "j" .. fmt(abs(im), d)
  end
  return fmt(re, d) .. (im < 0 and " - j" or " + j") .. fmt(abs(im), d)
end

--------------------------------------------------------------------------
-- polinomios: arreglo ascendente p[1]=a0, p[2]=a1 ...
--------------------------------------------------------------------------
function ptrim(p)
  local m = 0
  for i = 1, #p do if abs(p[i]) > m then m = abs(p[i]) end end
  local tol = m * 1e-13
  local n = #p
  while n > 1 and abs(p[n]) <= tol do p[n] = nil; n = n - 1 end
  if n == 0 then p[1] = 0 end
  return p
end

function pdeg(p) return #p - 1 end

function pcopy(p) local q = {} for i = 1, #p do q[i] = p[i] end return q end

function padd(a, b)
  local r = {}
  for i = 1, math.max(#a, #b) do r[i] = (a[i] or 0) + (b[i] or 0) end
  return ptrim(r)
end

function pscale(a, k)
  local r = {}
  for i = 1, #a do r[i] = a[i] * k end
  return ptrim(r)
end

function psub(a, b) return padd(a, pscale(b, -1)) end

function pmul(a, b)
  local r = {}
  for i = 1, #a + #b - 1 do r[i] = 0 end
  for i = 1, #a do
    local ai = a[i]
    if ai ~= 0 then
      for j = 1, #b do r[i + j - 1] = r[i + j - 1] + ai * b[j] end
    end
  end
  return ptrim(r)
end

function ppow(a, n)
  local r = {1}
  for i = 1, n do r = pmul(r, a) end
  return r
end

function piszero(p) return #p == 1 and p[1] == 0 end

function plead(p) return p[#p] end

-- evaluacion real
function peval(p, x)
  local r = 0
  for i = #p, 1, -1 do r = r * x + p[i] end
  return r
end

-- evaluacion compleja (Horner)
function pevalc(p, xr, xi)
  local rr, ri = 0, 0
  for i = #p, 1, -1 do
    rr, ri = rr * xr - ri * xi + p[i], rr * xi + ri * xr
  end
  return rr, ri
end

function pderiv(p)
  local r = {}
  for i = 2, #p do r[i - 1] = p[i] * (i - 1) end
  if #r == 0 then r[1] = 0 end
  return r
end

-- cuenta raices en el origen (coeficientes bajos nulos)
function pzeroroots(p)
  local m = 0
  for i = 1, #p do if abs(p[i]) > m then m = abs(p[i]) end end
  local n = 0
  while n < #p - 1 and abs(p[n + 1]) <= m * 1e-13 do n = n + 1 end
  return n
end

function pshiftdown(p, n)
  local r = {}
  for i = n + 1, #p do r[#r + 1] = p[i] end
  if #r == 0 then r[1] = 0 end
  return r
end

-- texto descendente:  s^3 + 11s^2 + 10s + 63
function pstr(p, v)
  v = v or "s"
  local parts = {}
  for i = #p, 1, -1 do
    local c = p[i]
    if abs(c) > 1e-12 or #p == 1 then
      local e = i - 1
      local mag = abs(c)
      local cs
      if e > 0 and abs(mag - 1) < 1e-12 then cs = "" else cs = fmt(mag) end
      local term
      if e == 0 then term = fmt(mag)
      elseif e == 1 then term = cs .. v
      else term = cs .. v .. "^" .. e end
      if #parts == 0 then
        parts[1] = (c < 0 and "-" or "") .. term
      else
        parts[#parts + 1] = (c < 0 and " - " or " + ") .. term
      end
    end
  end
  if #parts == 0 then return "0" end
  return table.concat(parts)
end

--------------------------------------------------------------------------
-- raices (Durand-Kerner + pulido de Newton + agrupacion de multiples)
--------------------------------------------------------------------------
local function polish(p, dp, zr, zi)
  for it = 1, 6 do
    local fr, fi = pevalc(p, zr, zi)
    local dr, di = pevalc(dp, zr, zi)
    local den = dr * dr + di * di
    if den < 1e-300 then break end
    local sr = (fr * dr + fi * di) / den
    local si = (fi * dr - fr * di) / den
    local nr, ni = zr - sr, zi - si
    local f2r, f2i = pevalc(p, nr, ni)
    if f2r * f2r + f2i * f2i < fr * fr + fi * fi then zr, zi = nr, ni else break end
  end
  return zr, zi
end

-- devuelve lista de {re, im}
function proots(pin, init)
  local p = ptrim(pcopy(pin))
  local roots = {}
  local nz = pzeroroots(p)
  for i = 1, nz do roots[#roots + 1] = {0, 0} end
  p = pshiftdown(p, nz)
  local n = #p - 1
  if n <= 0 then return roots end
  if n == 1 then
    roots[#roots + 1] = {-p[1] / p[2], 0}
  elseif n == 2 then
    local a, b, c = p[3], p[2], p[1]
    local disc = b * b - 4 * a * c
    local scale = b * b + abs(4 * a * c)
    if abs(disc) <= 1e-12 * scale then disc = 0 end
    if disc >= 0 then
      local q = -0.5 * (b + (b >= 0 and 1 or -1) * sqrt(disc))
      local r1 = q / a
      local r2 = (q ~= 0) and c / q or -b / (2 * a)
      if disc == 0 then r1 = -b / (2 * a); r2 = r1 end
      roots[#roots + 1] = {r1, 0}
      roots[#roots + 1] = {r2, 0}
    else
      local re, im = -b / (2 * a), sqrt(-disc) / (2 * abs(a))
      roots[#roots + 1] = {re, im}
      roots[#roots + 1] = {re, -im}
    end
  else
    -- monico
    local a = {}
    for i = 1, n + 1 do a[i] = p[i] / p[n + 1] end
    local R = 0
    for i = 1, n do
      local v = abs(a[i]) ^ (1 / (n - i + 1))
      if v > R then R = v end
    end
    R = math.max(R, 1e-3)
    local zr, zi = {}, {}
    local useinit = init and #init == n and nz == 0
    for k = 1, n do
      if useinit then
        -- continuacion: partir de las raices anteriores (perturbadas)
        zr[k] = init[k][1] + 1e-7 * k
        zi[k] = init[k][2] + 1.3e-7 * k
      else
        local ang = 2 * pi * (k - 1) / n + 0.4
        zr[k] = R * math.cos(ang)
        zi[k] = R * math.sin(ang)
      end
    end
    for it = 1, 800 do
      local maxstep = 0
      for k = 1, n do
        local fr, fi = pevalc(a, zr[k], zi[k])
        local qr, qi = 1, 0
        for j = 1, n do
          if j ~= k then
            local dr, di = zr[k] - zr[j], zi[k] - zi[j]
            qr, qi = qr * dr - qi * di, qr * di + qi * dr
          end
        end
        local den = qr * qr + qi * qi
        if den < 1e-300 then den = 1e-300 end
        local sr = (fr * qr + fi * qi) / den
        local si = (fi * qr - fr * qi) / den
        zr[k] = zr[k] - sr
        zi[k] = zi[k] - si
        local st = abs(sr) + abs(si)
        local sc = math.max(1, abs(zr[k]) + abs(zi[k]))
        if st / sc > maxstep then maxstep = st / sc end
      end
      if maxstep < 1e-13 then break end
    end
    local dp = pderiv(a)
    for k = 1, n do zr[k], zi[k] = polish(a, dp, zr[k], zi[k]) end
    -- agrupar raices multiples (promedio del cumulo)
    local used = {}
    for k = 1, n do
      if not used[k] then
        local cl = {k}
        local sc = math.max(1, abs(zr[k]) + abs(zi[k]))
        for j = k + 1, n do
          if not used[j] and abs(zr[j] - zr[k]) + abs(zi[j] - zi[k]) < 2e-3 * sc then
            cl[#cl + 1] = j
          end
        end
        if #cl > 1 then
          local mr, mi = 0, 0
          for _, j in ipairs(cl) do mr = mr + zr[j]; mi = mi + zi[j] end
          mr, mi = mr / #cl, mi / #cl
          for _, j in ipairs(cl) do zr[j] = mr; zi[j] = mi; used[j] = true end
        end
      end
    end
    for k = 1, n do roots[#roots + 1] = {zr[k], zi[k]} end
  end
  -- limpiar parte imaginaria y forzar pares conjugados
  for _, r in ipairs(roots) do
    local sc = math.max(1, abs(r[1]))
    if abs(r[2]) < 1e-7 * sc then r[2] = 0 end
  end
  for i, r in ipairs(roots) do
    if r[2] > 0 and not r.paired then
      local best, bd = nil, INF
      for j, q in ipairs(roots) do
        if j ~= i and q[2] < 0 and not q.paired then
          local d = abs(q[1] - r[1]) + abs(q[2] + r[2])
          if d < bd then bd, best = d, q end
        end
      end
      if best then
        local mr, mi = (r[1] + best[1]) / 2, (r[2] - best[2]) / 2
        r[1], r[2] = mr, mi
        best[1], best[2] = mr, -mi
        r.paired = true; best.paired = true
      end
    end
  end
  for _, r in ipairs(roots) do r.paired = nil end
  table.sort(roots, function(x, y)
    if abs(x[1] - y[1]) > 1e-9 * math.max(1, abs(x[1])) then return x[1] > y[1] end
    return x[2] > y[2]
  end)
  return roots
end

-- polinomio desde raices (coef. reales)
function pfromroots(roots, lead)
  local pr, pim = {1}, {0}
  for _, r in ipairs(roots) do
    local nr, ni = {}, {}
    for i = 1, #pr + 1 do nr[i] = 0; ni[i] = 0 end
    for i = 1, #pr do
      -- (s - r) * p
      nr[i + 1] = nr[i + 1] + pr[i]
      ni[i + 1] = ni[i + 1] + pim[i]
      nr[i] = nr[i] - (r[1] * pr[i] - r[2] * pim[i])
      ni[i] = ni[i] - (r[1] * pim[i] + r[2] * pr[i])
    end
    pr, pim = nr, ni
  end
  for i = 1, #pr do pr[i] = pr[i] * (lead or 1) end
  return ptrim(pr)
end

-- lista de raices en texto, agrupando pares conjugados
function rootsstr(roots, d)
  local out = {}
  local skip = {}
  for i, r in ipairs(roots) do
    if not skip[i] then
      if r[2] > 0 then
        if abs(r[1]) < 1e-11 then out[#out + 1] = U.pm .. "j" .. fmt(r[2], d)
        else out[#out + 1] = fmt(r[1], d) .. " " .. U.pm .. " j" .. fmt(r[2], d) end
        for j = i + 1, #roots do
          if not skip[j] and abs(roots[j][2] + r[2]) < 1e-9 * math.max(1, r[2]) and abs(roots[j][1] - r[1]) < 1e-9 * math.max(1, abs(r[1])) then
            skip[j] = true; break
          end
        end
      else
        out[#out + 1] = fmtc(r[1], r[2], d)
      end
    end
  end
  if #out == 0 then return "(ninguno)" end
  return table.concat(out, ";  ")
end

--------------------------------------------------------------------------
-- funciones racionales  {n=poly, d=poly}
--------------------------------------------------------------------------
function rat(n, d) return {n = ptrim(n), d = ptrim(d or {1})} end
function rconst(c) return rat({c}, {1}) end
function risconst(r) return #r.n == 1 and #r.d == 1 end
function rval(r) return r.n[1] / r.d[1] end

local function peq(a, b)
  if #a ~= #b then return false end
  local m = 0
  for i = 1, #a do m = math.max(m, abs(a[i]), abs(b[i])) end
  for i = 1, #a do if abs(a[i] - b[i]) > 1e-12 * m then return false end end
  return true
end

function radd(a, b)
  if peq(a.d, b.d) then return rat(padd(a.n, b.n), pcopy(a.d)) end
  return rat(padd(pmul(a.n, b.d), pmul(b.n, a.d)), pmul(a.d, b.d))
end
function rneg(a) return rat(pscale(a.n, -1), pcopy(a.d)) end
function rsub(a, b) return radd(a, rneg(b)) end
function rmul(a, b) return rat(pmul(a.n, b.n), pmul(a.d, b.d)) end
function rdiv(a, b)
  if piszero(b.n) then error("division por cero") end
  return rat(pmul(a.n, b.d), pmul(a.d, b.n))
end
function rpow(a, k)
  if k >= 0 then
    return rat(ppow(a.n, k), ppow(a.d, k))
  else
    if piszero(a.n) then error("division por cero") end
    return rat(ppow(a.d, -k), ppow(a.n, -k))
  end
end

-- normaliza: denominador monico y cancela factores comunes
function rsimplify(r)
  local n, d = ptrim(pcopy(r.n)), ptrim(pcopy(r.d))
  if piszero(n) then return rat({0}, {1}) end
  -- cancelar s^k
  local zn, zd = pzeroroots(n), pzeroroots(d)
  local m = math.min(zn, zd)
  if m > 0 then n = pshiftdown(n, m); d = pshiftdown(d, m) end
  -- cancelar raices comunes (numericamente)
  if #n > 1 and #d > 1 then
    local rn, rd = proots(n), proots(d)
    local usedn, usedd, any = {}, {}, false
    for i, a in ipairs(rn) do
      for j, b in ipairs(rd) do
        if not usedd[j] then
          local sc = math.max(1, abs(b[1]) + abs(b[2]))
          if abs(a[1] - b[1]) + abs(a[2] - b[2]) < 1e-6 * sc then
            usedn[i] = true; usedd[j] = true; any = true; break
          end
        end
      end
    end
    if any then
      local kn, kd = {}, {}
      for i, a in ipairs(rn) do if not usedn[i] then kn[#kn + 1] = a end end
      for j, b in ipairs(rd) do if not usedd[j] then kd[#kd + 1] = b end end
      n = pfromroots(kn, plead(n))
      d = pfromroots(kd, plead(d))
    end
  end
  local l = plead(d)
  return rat(pscale(n, 1 / l), pscale(d, 1 / l))
end

function rstr(r)
  local ns, ds = pstr(r.n), pstr(r.d)
  if #r.d == 1 and abs(r.d[1] - 1) < 1e-12 then return ns end
  return "(" .. ns .. ") / (" .. ds .. ")"
end

--------------------------------------------------------------------------
-- parser de expresiones en s (y parametro k)
--------------------------------------------------------------------------
local REPL = {
  {"\226\136\146", "-"}, {"\207\128", "pi"}, {"\195\151", "*"}, {"\194\183", "*"},
  {"\226\136\154", "sqrt"}, {"\225\180\135", "E"}, {"\194\178", "^2"}, {"\194\179", "^3"},
  {"\226\128\162", "*"}, {"\226\136\153", "*"}, {"\226\136\151", "*"},
  {"\195\183", "/"}, {"\226\136\149", "/"}, {"\226\128\147", "-"}, {"\226\128\148", "-"},
  {"\239\188\136", "("}, {"\239\188\137", ")"}, {"\194\160", " "},
}

function normalize(str)
  for _, rp in ipairs(REPL) do
    local a, b = rp[1], rp[2]
    local i = 1
    while true do
      local s, e = string.find(str, a, i, true)
      if not s then break end
      str = str:sub(1, s - 1) .. b .. str:sub(e + 1)
      i = s + #b
    end
  end
  return str
end

local FUNCS = {
  sqrt = math.sqrt, ln = math.log, log = math.log10, exp = math.exp,
  sin = function(x) return math.sin(x) end, cos = math.cos, tan = math.tan,
  abs = math.abs, atan = math.atan, asin = math.asin, acos = math.acos,
  sind = function(x) return math.sin(x * D2R) end,
  cosd = function(x) return math.cos(x * D2R) end,
  tand = function(x) return math.tan(x * D2R) end,
  rad = function(x) return x * D2R end, deg = function(x) return x * R2D end,
}

local function tokenize(str)
  local toks = {}
  local i, n = 1, #str
  while i <= n do
    local c = str:sub(i, i)
    if c:match("%s") then
      i = i + 1
    elseif c:match("[%d%.]") then
      local j = i
      while j <= n and str:sub(j, j):match("[%d%.]") do j = j + 1 end
      if j <= n and str:sub(j, j) == "E" then
        local k = j + 1
        if str:sub(k, k):match("[%+%-]") then k = k + 1 end
        if str:sub(k, k):match("%d") then
          while k <= n and str:sub(k, k):match("%d") do k = k + 1 end
          j = k
        end
      end
      local v = tonumber(str:sub(i, j - 1))
      if not v then error("numero invalido: " .. str:sub(i, j - 1)) end
      toks[#toks + 1] = {t = "num", v = v}
      i = j
    elseif c:match("%a") then
      local j = i
      while j <= n and str:sub(j, j):match("[%a_]") do j = j + 1 end
      local w = str:sub(i, j - 1):lower()
      if FUNCS[w] or w == "pi" then
        toks[#toks + 1] = {t = "id", v = w}
      else
        -- separar letras sueltas: "ks" -> k * s
        local k = 1
        while k <= #w do
          local rest = w:sub(k)
          local matched = false
          for _, nm in ipairs({"sqrt", "ln", "log", "exp", "pi"}) do
            if rest:sub(1, #nm) == nm then
              toks[#toks + 1] = {t = "id", v = nm}; k = k + #nm; matched = true; break
            end
          end
          if not matched then
            local ch = w:sub(k, k)
            if ch == "s" or ch == "k" or ch == "e" then
              toks[#toks + 1] = {t = "id", v = ch}
            else
              error("simbolo desconocido '" .. ch .. "' (use s, k, pi, e)")
            end
            k = k + 1
          end
        end
      end
      i = j
    elseif c:match("[%+%-%*/%^%(%)]") then
      toks[#toks + 1] = {t = "op", v = c}
      i = i + 1
    elseif c == "[" or c == "{" then
      toks[#toks + 1] = {t = "op", v = "("}; i = i + 1
    elseif c == "]" or c == "}" then
      toks[#toks + 1] = {t = "op", v = ")"}; i = i + 1
    else
      error("caracter no valido: '" .. c .. "'")
    end
  end
  return toks
end

-- parse(str, kval) -> rational ; usa_k (bool)
function parse(str, kval)
  if str == nil or str:match("^%s*$") then error("campo vacio") end
  local toks = tokenize(normalize(str))
  local pos = 1
  local usesk = false
  local function peek() return toks[pos] end
  local function nextt() pos = pos + 1; return toks[pos - 1] end
  local function isop(t, v) return t and t.t == "op" and t.v == v end
  local expr

  local function primary()
    local t = nextt()
    if not t then error("expresion incompleta") end
    if t.t == "num" then return rconst(t.v) end
    if t.t == "id" then
      if t.v == "s" then return rat({0, 1}, {1}) end
      if t.v == "k" then
        usesk = true
        if kval == nil then error("NEEDK") end
        return rconst(kval)
      end
      if t.v == "pi" then return rconst(pi) end
      if t.v == "e" then return rconst(math.exp(1)) end
      local f = FUNCS[t.v]
      if f then
        local arg
        if isop(peek(), "(") then
          nextt(); arg = expr()
          if not isop(nextt(), ")") then error("falta ')'") end
        else
          arg = primary()
        end
        if not risconst(arg) then
          error(t.v .. "() solo con constantes")
        end
        return rconst(f(rval(arg)))
      end
      error("simbolo " .. t.v)
    end
    if isop(t, "(") then
      local v = expr()
      if not isop(nextt(), ")") then error("falta ')'") end
      return v
    end
    error("simbolo inesperado '" .. tostring(t.v) .. "'")
  end

  local unary
  local function power()
    local b = primary()
    if isop(peek(), "^") then
      nextt()
      local e = unary()
      if not risconst(e) then error("exponente debe ser constante") end
      local ev = rval(e)
      if risconst(b) then return rconst(rval(b) ^ ev) end
      if ev ~= floor(ev) then error("exponente de s debe ser entero") end
      if abs(ev) > 30 then error("exponente muy grande") end
      return rpow(b, ev)
    end
    return b
  end
  unary = function()
    if isop(peek(), "-") then nextt(); return rneg(unary()) end
    if isop(peek(), "+") then nextt(); return unary() end
    return power()
  end
  local function term()
    local v = unary()
    while true do
      local t = peek()
      if isop(t, "*") then nextt(); v = rmul(v, unary())
      elseif isop(t, "/") then nextt(); v = rdiv(v, unary())
      elseif t and (t.t == "num" or t.t == "id" or isop(t, "(")) then v = rmul(v, power())
      else break end
    end
    return v
  end
  expr = function()
    local v = term()
    while true do
      local t = peek()
      if isop(t, "+") then nextt(); v = radd(v, term())
      elseif isop(t, "-") then nextt(); v = rsub(v, term())
      else break end
    end
    return v
  end
  local v = expr()
  if pos <= #toks then error("sobra texto despues de la expresion") end
  if piszero(v.d) then error("denominador nulo") end
  return v, usesk
end

-- detecta si la expresion contiene k (sin evaluarla)
function usesk(str)
  if not str or str:match("^%s*$") then return false end
  local ok, err = pcall(parse, str, nil)
  if not ok and tostring(err):find("NEEDK") then return true end
  return false
end

--------------------------------------------------------------------------
-- utilidades de sistemas
--------------------------------------------------------------------------
function isstable(roots, tol)
  tol = tol or 1e-9
  for _, r in ipairs(roots) do if r[1] >= -tol then return false end end
  return true
end

function stabilitytext(roots)
  local rhp, axis = 0, 0
  for _, r in ipairs(roots) do
    if r[1] > 1e-7 then rhp = rhp + 1
    elseif r[1] >= -1e-7 then axis = axis + 1 end
  end
  if rhp > 0 then return "INESTABLE (" .. rhp .. " polo(s) en SPD)", false end
  if axis > 0 then return "CRITICAMENTE ESTABLE (polos en eje j" .. U.w .. ")", false end
  return "ESTABLE (todos los polos en SPI)", true
end

-- lazo cerrado con retro unitaria: G/(1+G)
function closedloop(G)
  local n, d = G.n, G.d
  return rat(pcopy(n), padd(d, n))
end

-- angulo (grados) y magnitud de G en s = (sr + j si) usando raices
function gainphase_at(zeros, poles, lead, sr, si)
  local mag, ph = abs(lead), (lead < 0) and -180 or 0
  for _, z in ipairs(zeros) do
    local dr, di = sr - z[1], si - z[2]
    mag = mag * sqrt(dr * dr + di * di)
    ph = ph + atan2(di, dr) * R2D
  end
  for _, p in ipairs(poles) do
    local dr, di = sr - p[1], si - p[2]
    mag = mag / sqrt(dr * dr + di * di)
    ph = ph - atan2(di, dr) * R2D
  end
  return mag, ph
end

-- objeto de respuesta en frecuencia (con retardo L)
function freqobj(G, L)
  local o = {G = G, L = L or 0}
  o.z = proots(G.n)
  o.p = proots(G.d)
  o.lead = plead(G.n) / plead(G.d)
  return o
end

-- magnitud lineal y fase (grados, continua) en w
function freqat(o, w)
  local mag, ph = gainphase_at(o.z, o.p, o.lead, 0, w)
  ph = ph - o.L * w * R2D
  return mag, ph
end

-- busca cruces de ganancia y de fase
function margins(o)
  local res = {gc = {}, pc = {}}
  local N = 1600
  local lw0, lw1 = -4, 4
  local function f(lw)
    local m, ph = freqat(o, 10 ^ lw)
    return m, ph
  end
  local prevm, prevp, prevlw
  for i = 0, N do
    local lw = lw0 + (lw1 - lw0) * i / N
    local m, ph = f(lw)
    local gm = (m > 0) and log(m) or -1e9
    local pv = (ph + 180) / 360
    if prevlw then
      -- cruce de ganancia |G|=1
      if (prevm > 0) ~= (gm > 0) and prevm ~= 0 then
        local a, b = prevlw, lw
        for it = 1, 60 do
          local c = (a + b) / 2
          local mc = f(c)
          local gmc = log(mc)
          if (gmc > 0) == (prevm > 0) then a = c else b = c end
        end
        local w = 10 ^ ((a + b) / 2)
        local mm, phh = freqat(o, w)
        res.gc[#res.gc + 1] = {w = w, ph = phh, MF = 180 + phh}
      end
      -- cruce de fase (fase = -180 + 360m)
      local k0, k1 = floor(prevp), floor(pv)
      if k0 ~= k1 then
        local target = (k0 > k1) and k0 or k1   -- multiplo cruzado
        local a, b = prevlw, lw
        local sa = (prevp - target) > 0
        for it = 1, 60 do
          local c = (a + b) / 2
          local _, phc = f(c)
          if (((phc + 180) / 360 - target) > 0) == sa then a = c else b = c end
        end
        local w = 10 ^ ((a + b) / 2)
        local mm, phh = freqat(o, w)
        res.pc[#res.pc + 1] = {w = w, mag = mm, ph = phh, MG = 1 / mm}
      end
    end
    prevm, prevp, prevlw = gm, pv, lw
  end
  return res
end

function db(x) return 20 * math.log10(x) end

--------------------------------------------------------------------------
-- Routh-Hurwitz numerico
--------------------------------------------------------------------------
-- devuelve tabla (lista de filas), notas, cambios de signo
function routh(pin)
  local p = ptrim(pcopy(pin))
  local n = #p - 1
  local notes = {}
  if n < 1 then return {}, {"grado < 1"}, 0 end
  -- descendente
  local a = {}
  for i = n + 1, 1, -1 do a[#a + 1] = p[i] end
  local w = floor(n / 2) + 1
  local rows = {}
  local r1, r2 = {}, {}
  for j = 1, w do
    r1[j] = a[2 * j - 1] or 0
    r2[j] = a[2 * j] or 0
  end
  rows[1] = r1; rows[2] = r2
  local scale = 0
  for i = 1, #a do scale = math.max(scale, abs(a[i])) end
  local tol = 1e-10 * math.max(scale, 1e-300)
  local epsused, auxused = false, false
  local epscells = {}
  for i = 3, n + 1 do
    local A, B = rows[i - 2], rows[i - 1]
    -- fila nula?
    local allz = true
    for j = 1, w do if abs(B[j]) > tol then allz = false end end
    if allz then
      local pw = n - (i - 2) + 1 -- potencia de la fila A
      local pa = {}
      for j = 1, w do pa[j] = A[j] end
      local newB = {}
      for j = 1, w do newB[j] = A[j] * (pw - 2 * (j - 1)) end
      rows[i - 1] = newB
      B = newB
      auxused = true
      local auxp = {}
      for j = 1, w do
        local e = pw - 2 * (j - 1)
        if e >= 0 then auxp[e + 1] = A[j] end
      end
      for e = 1, pw + 1 do auxp[e] = auxp[e] or 0 end
      notes[#notes + 1] = "Fila s^" .. (pw - 1) .. " nula " .. U.arrow .. " pol. auxiliar A(s) = " .. pstr(ptrim(auxp))
      notes[#notes + 1] = "  raices de A(s): " .. rootsstr(proots(auxp))
      notes[#notes + 1] = "  se reemplaza la fila por dA/ds"
    end
    if abs(B[1]) <= tol then
      B[1] = 1e-6 * math.max(1, scale)
      epsused = true
      epscells[i - 1] = true
      notes[#notes + 1] = "Cero en 1a columna (fila s^" .. (n - (i - 2)) .. ") " .. U.arrow .. " se reemplaza por " .. U.eps .. "=0+"
    end
    local C = {}
    for j = 1, w do
      local a2 = A[j + 1] or 0
      local b2 = B[j + 1] or 0
      C[j] = (B[1] * a2 - A[1] * b2) / B[1]
      if abs(C[j]) <= tol then C[j] = 0 end
    end
    rows[i] = C
  end
  local ch = 0
  local prev = nil
  for i = 1, n + 1 do
    local v = rows[i][1]
    if v ~= 0 then
      local sgn = v > 0
      if prev ~= nil and sgn ~= prev then ch = ch + 1 end
      prev = sgn
    end
  end
  return rows, notes, ch, epsused, auxused, epscells
end

-- prueba rapida de estabilidad asintotica (Routh estricto)
function routhstable(pin)
  local p = ptrim(pcopy(pin))
  local n = #p - 1
  if n < 1 then return false end
  local a = {}
  for i = n + 1, 1, -1 do a[#a + 1] = p[i] end
  local sg = a[1] > 0 and 1 or -1
  for i = 1, #a do if a[i] * sg <= 0 then return false end end
  local w = floor(n / 2) + 1
  local A, B = {}, {}
  for j = 1, w do A[j] = a[2 * j - 1] or 0; B[j] = a[2 * j] or 0 end
  local scale = 0
  for i = 1, #a do scale = math.max(scale, abs(a[i])) end
  for i = 3, n + 1 do
    if B[1] * sg <= 1e-12 * scale then return false end
    local C = {}
    for j = 1, w do C[j] = (B[1] * (A[j + 1] or 0) - A[1] * (B[j + 1] or 0)) / B[1] end
    A, B = B, C
  end
  return B[1] * sg > 1e-12 * scale
end

-- rango de k estable; fpoly(k) -> polinomio caracteristico
function krange(fpoly)
  local ks = {}
  for e = 50, -40, -1 do ks[#ks + 1] = -(10 ^ (e / 10)) end
  ks[#ks + 1] = 0
  for e = -40, 50 do ks[#ks + 1] = 10 ^ (e / 10) end
  -- refinar la grilla
  local grid = {}
  for i = 1, #ks - 1 do
    local a, b = ks[i], ks[i + 1]
    for j = 0, 9 do grid[#grid + 1] = a + (b - a) * j / 10 end
  end
  grid[#grid + 1] = ks[#ks]
  local function st(k)
    local ok, p = pcall(fpoly, k)
    if not ok then return false end
    return routhstable(p)
  end
  local prev = st(grid[1])
  local bounds = {}
  local start = prev and grid[1] or nil
  local intervals = {}
  for i = 2, #grid do
    local cur = st(grid[i])
    if cur ~= prev then
      local a, b = grid[i - 1], grid[i]
      for it = 1, 70 do
        local c = (a + b) / 2
        if st(c) == prev then a = c else b = c end
      end
      local kb = (a + b) / 2
      if abs(kb) < 1e-8 then kb = 0 end
      bounds[#bounds + 1] = kb
      if cur then start = kb else intervals[#intervals + 1] = {start, kb}; start = nil end
    end
    prev = cur
  end
  if prev and start then intervals[#intervals + 1] = {start, INF} end
  -- extremos de la grilla
  for _, iv in ipairs(intervals) do
    if iv[1] == grid[1] then iv[1] = -INF end
  end
  return intervals, bounds
end

--------------------------------------------------------------------------
-- errores en estado estacionario (retro unitaria)
--------------------------------------------------------------------------
function sserrors(G0, A)
  A = A or 1
  local n, d = G0.n, G0.d
  local zn, zd = pzeroroots(n), pzeroroots(d)
  local N = zd - zn
  local nn, dd = pshiftdown(n, zn), pshiftdown(d, zd)
  local Kst = nn[1] / dd[1]
  local Kp, Kv, Ka
  if N < 0 then Kp, Kv, Ka = 0, 0, 0
  elseif N == 0 then Kp, Kv, Ka = Kst, 0, 0
  elseif N == 1 then Kp, Kv, Ka = INF, Kst, 0
  elseif N == 2 then Kp, Kv, Ka = INF, INF, Kst
  else Kp, Kv, Ka = INF, INF, INF end
  local function inv(x) if x == INF then return 0 end if x == 0 then return INF end return A / x end
  local ep = (Kp == INF) and 0 or A / (1 + Kp)
  local ev = inv(Kv)
  local ea = inv(Ka)
  if Kv == 0 then ev = INF end
  if Ka == 0 then ea = INF end
  return {N = N, Kst = Kst, Kp = Kp, Kv = Kv, Ka = Ka, ep = ep, ev = ev, ea = ea}
end

--------------------------------------------------------------------------
-- especificaciones de segundo orden
--------------------------------------------------------------------------
function zeta_from_Mp(Mp)
  local l = log(Mp / 100)
  return -l / sqrt(pi * pi + l * l)
end

function Mp_from_zeta(z)
  if z >= 1 then return 0 end
  if z <= 0 then return 100 end
  return 100 * exp(-pi * z / sqrt(1 - z * z))
end

--------------------------------------------------------------------------
-- fracciones parciales (residuos con polos multiples)
--------------------------------------------------------------------------
local function cmul(ar, ai, br, bi) return ar * br - ai * bi, ar * bi + ai * br end
local function cdiv(ar, ai, br, bi)
  local d = br * br + bi * bi
  return (ar * br + ai * bi) / d, (ai * br - ar * bi) / d
end

-- coeficientes de Taylor de un polinomio real p en s0 (complejo), m terminos
local function taylor(p, sr, si, m)
  local cr, ci = {}, {}
  for i = 1, #p do cr[i] = p[i]; ci[i] = 0 end
  local out = {}
  local n = #p
  for k = 1, m do
    -- division sintetica: el residuo es el coef. k-1
    if n < 1 then out[k] = {0, 0}
    else
      local br, bi = cr[n], ci[n]
      local qr, qi = {}, {}
      for i = n - 1, 1, -1 do
        qr[i], qi[i] = br, bi
        local tr, ti = cmul(br, bi, sr, si)
        br, bi = cr[i] + tr, ci[i] + ti
      end
      out[k] = {br, bi}
      cr, ci = qr, qi
      n = n - 1
    end
  end
  return out
end

-- Y = N/D (estrictamente propia). Devuelve grupos {p={re,im}, m, c={...}}
-- c[j] = coeficiente de 1/(s-p)^j
function partialfractions(N, D)
  if #N >= #D then error("Y(s) no es estrictamente propia") end
  local roots = proots(D)
  local groups = {}
  for _, r in ipairs(roots) do
    local found = nil
    for _, g in ipairs(groups) do
      local sc = math.max(1, abs(g.p[1]) + abs(g.p[2]))
      if abs(g.p[1] - r[1]) + abs(g.p[2] - r[2]) < 1e-6 * sc then found = g; break end
    end
    if found then found.m = found.m + 1 else groups[#groups + 1] = {p = {r[1], r[2]}, m = 1} end
  end
  local lead = plead(D)
  for _, g in ipairs(groups) do
    local m = g.m
    local pr, pi_ = g.p[1], g.p[2]
    -- serie del denominador sin (s-p)^m
    local dr, di = {lead}, {0}
    for k = 2, m do dr[k] = 0; di[k] = 0 end
    for _, h in ipairs(groups) do
      if h ~= g then
        local c_r, c_i = pr - h.p[1], pi_ - h.p[2]
        for rep = 1, h.m do
          -- multiplicar serie por (c + u)
          local nr, ni = {}, {}
          for k = 1, m do
            local ar, ai = cmul(dr[k], di[k], c_r, c_i)
            if k > 1 then ar = ar + dr[k - 1]; ai = ai + di[k - 1] end
            nr[k], ni[k] = ar, ai
          end
          dr, di = nr, ni
        end
      end
    end
    local nt = taylor(N, pr, pi_, m)
    local qr, qi = {}, {}
    for k = 1, m do
      local ar, ai = nt[k][1], nt[k][2]
      for i = 2, k do
        local tr, ti = cmul(dr[i], di[i], qr[k - i + 1], qi[k - i + 1])
        ar, ai = ar - tr, ai - ti
      end
      qr[k], qi[k] = cdiv(ar, ai, dr[1], di[1])
    end
    g.c = {}
    for k = 1, m do
      local j = m - k + 1
      g.c[j] = {qr[k], qi[k]}
      local sc = math.max(1, abs(qr[k]) + abs(qi[k]))
      if abs(g.c[j][2]) < 1e-9 * sc and g.p[2] == 0 then g.c[j][2] = 0 end
    end
  end
  return groups
end


--------------------------------------------------------------------------
-- respuesta temporal exacta a partir de fracciones parciales
--------------------------------------------------------------------------
-- y(t) = sum Re( c_j t^(j-1)/(j-1)! e^(p t) )  (todos los polos, incl. conjugados)
function evalpf(groups, t)
  local y = 0
  for _, g in ipairs(groups) do
    local er = exp(g.p[1] * t)
    local cr, ci = math.cos(g.p[2] * t), math.sin(g.p[2] * t)
    local tp, f = 1, 1
    for j = 1, g.m do
      if j > 1 then tp = tp * t; f = f * (j - 1) end
      local c = g.c[j]
      -- Re((cr_ + j ci_) * (cos + j sin)) * e^(sigma t)
      y = y + er * (c[1] * cr - c[2] * ci) * tp / f
    end
  end
  return y
end

-- respuesta al escalon de amplitud A de T (racional propia)
function stepresp(T, A, tmax, N)
  local groups = partialfractions(pscale(T.n, A), pmul(T.d, {0, 1}))
  local ts, ys = {}, {}
  for i = 0, N do
    local t = tmax * i / N
    ts[#ts + 1] = t
    ys[#ys + 1] = evalpf(groups, t)
  end
  return ts, ys
end

-- horizonte de simulacion razonable segun los polos
function horizon(poles)
  local dom = nil
  local fastest = 0
  for _, r in ipairs(poles) do
    if r[1] < -1e-9 then
      if not dom or -r[1] < dom then dom = -r[1] end
    end
    if abs(r[1]) > fastest then fastest = abs(r[1]) end
    if r[2] ~= 0 and abs(r[2]) > fastest then fastest = abs(r[2]) end
  end
  local unst = 0
  for _, r in ipairs(poles) do if r[1] > unst then unst = r[1] end end
  if unst > 1e-9 then return 4 / unst end
  if not dom then return 20 end
  return 1.6 * 4 / dom
end

--------------------------------------------------------------------------
-- simulacion en el tiempo con retardo (RK4, retenedor en cada paso)
--------------------------------------------------------------------------
-- realizacion canonica controlable de una racional propia
function ssreal(R)
  local d = R.d
  local m = #d - 1
  local l = d[#d]
  local a = {}
  for i = 1, m + 1 do a[i] = d[i] / l end
  local n = {}
  for i = 1, #R.n do n[i] = R.n[i] / l end
  if #n > m + 1 then error("Gc impropio: no se puede simular") end
  local D = n[m + 1] or 0
  local c = {}
  for i = 1, m do c[i] = (n[i] or 0) - D * a[i] end
  return {m = m, a = a, c = c, D = D}
end

local function ssderiv(S, x, u, dx)
  local m = S.m
  for i = 1, m - 1 do dx[i] = x[i + 1] end
  if m > 0 then
    local s_ = u
    for i = 1, m do s_ = s_ - S.a[i] * x[i] end
    dx[m] = s_
  end
end

local function ssout(S, x, u)
  local y = S.D * u
  for i = 1, S.m do y = y + S.c[i] * x[i] end
  return y
end

local function rk4(S, x, u, h)
  local m = S.m
  if m == 0 then return end
  if not S.k1 then S.k1, S.k2, S.k3, S.k4, S.xt = {}, {}, {}, {}, {} end
  local k1, k2, k3, k4, xt = S.k1, S.k2, S.k3, S.k4, S.xt
  ssderiv(S, x, u, k1)
  for i = 1, m do xt[i] = x[i] + 0.5 * h * k1[i] end
  ssderiv(S, xt, u, k2)
  for i = 1, m do xt[i] = x[i] + 0.5 * h * k2[i] end
  ssderiv(S, xt, u, k3)
  for i = 1, m do xt[i] = x[i] + h * k3[i] end
  ssderiv(S, xt, u, k4)
  for i = 1, m do x[i] = x[i] + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) end
end

-- lazo cerrado unitario: e = A - y ; u = Gc(e) ; y = Gp(u(t - L))
function simloop(Gc, Gp, L, A, tmax, nsteps)
  local Sc, Sp = ssreal(Gc), ssreal(Gp)
  local h = tmax / nsteps
  local nd = math.floor((L or 0) / h + 0.5)
  local xc, xp = {}, {}
  for i = 1, Sc.m do xc[i] = 0 end
  for i = 1, Sp.m do xp[i] = 0 end
  local buf = {}
  for i = 1, nd do buf[i] = 0 end
  local head = 1
  local ts, ys, us = {0}, {0}, {}
  local y = ssout(Sp, xp, 0)
  ys[1] = y
  for k = 1, nsteps do
    local e = A - y
    local u = ssout(Sc, xc, e)
    us[#us + 1] = u
    local ud = u
    if nd > 0 then ud = buf[head]; buf[head] = u; head = head % nd + 1 end
    rk4(Sc, xc, e, h)
    rk4(Sp, xp, ud, h)
    y = ssout(Sp, xp, ud)
    if y ~= y or abs(y) > 1e12 then break end
    ts[#ts + 1] = k * h
    ys[#ys + 1] = y
  end
  return ts, ys, us
end

-- metricas de una respuesta al escalon
function stepmetrics(ts, ys, yss)
  local n = #ys
  local res = {yss = yss}
  if not yss or yss ~= yss or abs(yss) < 1e-12 then return res end
  local ymax, imax = -INF, 1
  for i = 1, n do if ys[i] > ymax then ymax, imax = ys[i], i end end
  if yss > 0 then
    res.Mp = math.max(0, (ymax - yss) / abs(yss) * 100)
  else
    local ymin, imin = INF, 1
    for i = 1, n do if ys[i] < ymin then ymin, imin = ys[i], i end end
    res.Mp = math.max(0, (yss - ymin) / abs(yss) * 100)
    imax = imin
  end
  res.tp = (res.Mp > 0.05) and ts[imax] or nil
  for i = 2, n do
    if (ys[i] - yss) * (ys[1] - yss) <= 0 then res.tr = ts[i]; break end
  end
  for _, tol in ipairs({2, 5}) do
    local band = tol / 100 * abs(yss)
    local last = nil
    for i = n, 1, -1 do
      if abs(ys[i] - yss) > band then last = i; break end
    end
    if last == nil then res["ts" .. tol] = 0
    elseif last < n then res["ts" .. tol] = ts[last + 1]
    else res["ts" .. tol] = nil end
  end
  return res
end

-- submuestreo para graficar
function decimate(ts, ys, maxn)
  local n = #ts
  if n <= maxn then return ts, ys end
  local t2, y2 = {}, {}
  for i = 1, maxn do
    local k = math.floor((i - 1) * (n - 1) / (maxn - 1)) + 1
    t2[i] = ts[k]; y2[i] = ys[k]
  end
  return t2, y2
end

-- ===== modules.lua =====
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

-- graficos asociados al resultado (tecla tab)
local function view(L, v) L.views = L.views or {}; L.views[#L.views + 1] = v end

-- recuadro de respuestas al inicio del resultado
local function say(L, txt) L.say = L.say or {}; L.say[#L.say + 1] = txt end

-- resultados clave para comparar (tecla c)
local function metric(L, key, label, val, unit)
  L.met = L.met or {}
  for _, m in ipairs(L.met) do if m.key == key then m.val = val; return end end
  L.met[#L.met + 1] = {key = key, label = label, val = val, unit = unit or ""}
end

-- frases individuales de especificaciones de 2o orden
local function sayspecs(L, z, wn)
  local Zs, WNs, S_ = U.zeta, U.wn, U.sigma
  local sg = z * wn
  metric(L, "zeta", Zs, z)
  metric(L, "wn", WNs, wn, "rad/s")
  metric(L, "sigma", S_ .. " = " .. Zs .. WNs, sg)
  metric(L, "ts2", "ts(2%)", 4 / sg, "s")
  metric(L, "ts5", "ts(5%)", 3 / sg, "s")
  if z < 1 then
    local wd_ = wn * sqrt(1 - z * z)
    metric(L, "wd", U.w .. "d", wd_, "rad/s")
    metric(L, "tr", "tr", (pi - math.acos(z)) / wd_, "s")
    metric(L, "tp", "tp", pi / wd_, "s")
    metric(L, "Mp", "Mp", Mp_from_zeta(z), "%")
  else
    metric(L, "Mp", "Mp", 0, "%")
  end
  say(L, "El coeficiente de amortiguamiento es " .. Zs .. " = " .. fmt(z, 4) .. " y la frecuencia natural es " .. WNs .. " = " .. fmt(wn, 4) .. " rad/s.")
  say(L, "El tiempo de asentamiento al 2% es ts(2%) = 4/(" .. Zs .. WNs .. ") = 4/" .. fmt(sg, 4) .. " = " .. fmt(4 / sg, 4) .. " s.")
  say(L, "El tiempo de asentamiento al 5% es ts(5%) = 3/(" .. Zs .. WNs .. ") = 3/" .. fmt(sg, 4) .. " = " .. fmt(3 / sg, 4) .. " s.")
  if z < 1 then
    local wd = wn * sqrt(1 - z * z)
    say(L, "La frecuencia natural amortiguada es " .. U.w .. "d = " .. WNs .. U.sqrt .. "(1-" .. Zs .. U.sq .. ") = " .. fmt(wd, 4) .. " rad/s.")
    say(L, "El tiempo de subida es tr = (" .. U.pi .. " - acos " .. Zs .. ")/" .. U.w .. "d = (" .. U.pi .. " - " .. fmt(math.acos(z), 4) .. ")/" .. fmt(wd, 4) .. " = " .. fmt((pi - math.acos(z)) / wd, 4) .. " s.")
    say(L, "El tiempo de peak es tp = " .. U.pi .. "/" .. U.w .. "d = " .. U.pi .. "/" .. fmt(wd, 4) .. " = " .. fmt(pi / wd, 4) .. " s.")
    say(L, "El sobrepaso maximo es Mp = 100" .. "\194\183" .. "e^(-" .. Zs .. U.pi .. "/" .. U.sqrt .. "(1-" .. Zs .. U.sq .. ")) = " .. fmt(Mp_from_zeta(z), 4) .. " %.")
  else
    say(L, "Como " .. Zs .. " " .. U.ge .. " 1 la respuesta es sobreamortiguada: no hay sobrepaso (Mp = 0 %) ni tiempo de peak.")
  end
end

local function summary(L, items)
  local out = {"== RESPUESTAS =="}
  for _, it in ipairs(items) do out[#out + 1] = it end
  if L.say then
    out[#out + 1] = "== Respuesta redactada =="
    for _, t in ipairs(L.say) do out[#out + 1] = t end
  end
  if L.chainhelp then
    for _, h in ipairs(L.chainhelp) do out[#out + 1] = h end
  end
  if L.met then out[#out + 1] = "Tecla c: modificar datos y COMPARAR con este resultado" end
  out[#out + 1] = "== Desarrollo =="
  for i = #out, 1, -1 do table.insert(L, 1, out[i]) end
end

-- enviar datos a otro modulo con una tecla numerica
local function chain(L, key, target, vals, desc)
  L.chain = L.chain or {}
  L.chain[key] = {target = target, vals = vals}
  L.chainhelp = L.chainhelp or {}
  L.chainhelp[#L.chainhelp + 1] = "Tecla " .. key .. ": " .. desc
end

-- polinomio y racional en texto exacto (para reingresar en otro modulo)
local function pexact(p)
  local parts = {}
  for i = #p, 1, -1 do
    local c = p[i]
    if c ~= 0 then
      local e = i - 1
      local cs = string.format("%.7g", c)
      local term = cs .. ((e == 0) and "" or ((e == 1) and "*s" or ("*s^" .. e)))
      parts[#parts + 1] = term
    end
  end
  if #parts == 0 then return "0" end
  return (table.concat(parts, "+"):gsub("%+%-", "-"))
end
local function rexact(R) return "(" .. pexact(R.n) .. ")/(" .. pexact(R.d) .. ")" end

-- vista de respuesta al escalon de T (lazo ya cerrado), via fracciones parciales
local function timeview(T, A, title, name)
  local poles = proots(T.d)
  local tmax = horizon(poles)
  local ts, ys = stepresp(T, A, tmax, 300)
  local st = isstable(poles)
  local yss = (st and abs(T.d[1]) > 1e-12) and A * T.n[1] / T.d[1] or nil
  local m = stepmetrics(ts, ys, yss)
  local ymax = -INF
  for i = 1, #ys do if ys[i] > ymax then ymax = ys[i] end end
  local V = {kind = "time", title = title, ref = yss, band = yss ~= nil,
    curves = {{t = ts, y = ys, name = name or "y(t)"}},
    marks = {Mp = m.Mp, tp = m.tp, ymax = ymax, ts = m.ts2}}
  return V, m, st
end

local function metricstext(m)
  local t = {}
  if m.Mp then t[#t + 1] = "Mp = " .. fmt(m.Mp, 4) .. " %" end
  if m.tp then t[#t + 1] = "tp = " .. fmt(m.tp, 4) .. " s" end
  if m.ts2 then t[#t + 1] = "ts2% = " .. fmt(m.ts2, 4) .. " s" end
  if m.ts5 then t[#t + 1] = "ts5% = " .. fmt(m.ts5, 4) .. " s" end
  if m.tr then t[#t + 1] = "tr = " .. fmt(m.tr, 4) .. " s" end
  return table.concat(t, ", ")
end

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
    -- justificacion: cual polo domina y por que
    add(L, "-- Eleccion del polo dominante --")
    add(L, "Cada polo p aporta un modo e^(Re(p)" .. "\194\183" .. "t): mientras menor es |Re(p)|, mas lento decae.")
    add(L, "Tabla (ordenada por cercania al eje j" .. W .. "):")
    local seen = {}
    for _, q in ipairs(poles) do
      local key = fmt(q[1], 6) .. "," .. fmt(abs(q[2]), 6)
      if not seen[key] then
        seen[key] = true
        local name = (q[2] ~= 0) and (fmt(q[1], 4) .. " " .. PM .. " j" .. fmt(abs(q[2]), 4)) or fmt(q[1], 4)
        add(L, "  " .. name .. ":  |Re| = " .. fmt(abs(q[1]), 4) .. ",  " .. TAU .. " = 1/|Re| = " .. fmt(1 / abs(q[1]), 4) .. " s,  se extingue en " .. U.approx .. "4" .. TAU .. " = " .. fmt(4 / abs(q[1]), 4) .. " s")
      end
    end
    add(L, "El dominante es el de MENOR |Re| (mas cercano al eje j" .. W .. "): su modo")
    add(L, "es el ultimo en extinguirse, por eso gobierna la respuesta (ts, Mp).")
    add(L, "Los demas desaparecen antes; se desprecian si |Re| " .. U.ge .. " 5" .. "\194\183" .. "|Re dom|.")
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
  local S = {Z .. " = " .. fmt(z) .. "   Mp = " .. fmt(Mp_from_zeta(z), 4) .. " %"}
  if wn then
    S[#S + 1] = WN .. " = " .. fmt(wn) .. " rad/s,  " .. SIG .. " = " .. fmt(sig)
    if z < 1 and wd then
      S[#S + 1] = "pd = -" .. fmt(sig) .. " " .. PM .. " j" .. fmt(wd)
      S[#S + 1] = "ts2% = " .. fmt(4 / sig, 4) .. " s, ts5% = " .. fmt(3 / sig, 4) .. " s, tp = " .. fmt(pi / wd, 4) .. " s, tr = " .. fmt((pi - math.acos(z)) / wd, 4) .. " s"
    end
    local V = timeview(rat({wn * wn}, {wn * wn, 2 * sig, 1}), 1, "2o orden estandar", "y(t) escalon unitario")
    view(L, V)
    sayspecs(L, z, wn)
    if z < 1 and wd then
      say(L, "Los polos deseados de lazo cerrado son s = -" .. fmt(sig, 4) .. " " .. PM .. " j" .. fmt(wd, 4) .. ", es decir, el polinomio deseado es s" .. U.sq .. " + " .. fmt(2 * sig, 4) .. "s + " .. fmt(wn * wn, 4) .. ".")
    end
  else
    say(L, "El coeficiente de amortiguamiento es " .. Z .. " = " .. fmt(z, 4) .. ".")
  end
  summary(L, S)
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
  if lazo == 1 then
    say(L, "La funcion de transferencia de lazo cerrado en funcion de K es GLC(s) = Y(s)/R(s) = [" .. polyk_str(N0, N1) .. "] / [" .. polyk_str(C0, C1) .. "].")
    local n = math.max(#C0, #C1)
    if n == 3 and abs(C1[3] or 0) < 1e-12 and abs(C1[2] or 0) < 1e-12 and abs(C0[3] or 0) > 1e-12 then
      local l = C0[3]
      say(L, "Comparando con la forma estandar, " .. WN .. U.sq .. " = " .. linstr((C0[1] or 0) / l, (C1[1] or 0) / l) .. " y 2" .. Z .. WN .. " = " .. fmt((C0[2] or 0) / l) .. "; como " .. Z .. WN .. " = " .. fmt((C0[2] or 0) / l / 2) .. " no depende de K, el tiempo de asentamiento ts(2%) = 4/(" .. Z .. WN .. ") = " .. fmt(8 / ((C0[2] or 0) / l)) .. " s no cambia al aumentar K, mientras que " .. WN .. " aumenta, " .. Z .. " disminuye y el sobrepaso Mp aumenta (respuesta mas oscilatoria).")
    end
    summary(L, {"GLC(s) = [" .. polyk_str(N0, N1) .. "] / [" .. polyk_str(C0, C1) .. "]"})
  end
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
  do
    local V, m, st = timeview(T, A, (lazo == 1) and "Respuesta LC al escalon" or "Respuesta al escalon", "y(t), A = " .. fmt(A))
    if #poles >= 1 then
      local dom = poles[1]
      local nd = (dom[2] ~= 0) and 2 or 1
      local PV = {kind = "poles", poles = poles, zeros = zeros, dom = dom, nd = nd, legend = {}}
      if isstable(poles) then
        PV.legend[#PV.legend + 1] = {"azul: dominante(s) " .. fmtc(dom[1], abs(dom[2]), 4) .. ((nd == 2) and " (par)" or "") .. ", |Re| = " .. fmt(abs(dom[1]), 4), "pd"}
        PV.legend[#PV.legend + 1] = {"linea verde: Re = 5" .. "\194\183" .. "Re dom = " .. fmt(5 * dom[1], 4) .. " (a la izq. se desprecian)", "cz"}
        if poles[nd + 1] then PV.legend[#PV.legend + 1] = {"siguiente polo: |Re| = " .. fmt(abs(poles[nd + 1][1]), 4) .. "  " .. U.arrow .. " relacion " .. fmt(abs(poles[nd + 1][1]) / abs(dom[1]), 3), "txt"} end
      else
        PV.legend[#PV.legend + 1] = {"hay polos en el SPD o en el eje j" .. W .. ": sistema no estable", "pp"}
      end
      view(L, PV)
    end
    view(L, V)
    local S = {((lazo == 1) and "GLC = " or "G = ") .. rstr(T), "Polos: " .. rootsstr(poles, 4), (stabilitytext(poles))}
    if st and m.yss then
      add(L, "== Simulacion exacta (fracciones parciales) ==")
      add(L, "y(" .. U.inf .. ") = " .. fmt(m.yss) .. ";  " .. metricstext(m))
      add(L, "(valores reales; las formulas de 2o orden son aproximaciones)")
      S[#S + 1] = "Real: " .. metricstext(m)
    end
    do
      say(L, ((lazo == 1) and "La funcion de transferencia de lazo cerrado es GLC(s) = Y(s)/R(s) = " or "La funcion de transferencia es G(s) = ") .. rstr(T) .. ".")
      metric(L, "estable", "estabilidad", st and 1 or 0)
      metric(L, "polos", "polos LC", rootsstr(poles, 4))
      say(L, "Su ecuacion caracteristica es " .. pstr(T.d) .. " = 0, cuyos polos son " .. rootsstr(poles, 4) .. "; por lo tanto el sistema es " .. (st and "estable (todos los polos en el semiplano izquierdo)." or "NO estable."))
      if st and #poles >= 1 then
        local dom = poles[1]
        local nd = (dom[2] ~= 0) and 2 or 1
        if #poles > nd then
          local ratio = abs(poles[nd + 1][1]) / abs(dom[1])
          metric(L, "ratio", "separacion de polos", ratio)
          do
            local nm = (nd == 2) and (fmt(dom[1], 4) .. " " .. PM .. " j" .. fmt(abs(dom[2]), 4)) or fmt(dom[1], 4)
            local nn = poles[nd + 1]
            local nn_s = (nn[2] ~= 0) and (fmt(nn[1], 4) .. " " .. PM .. " j" .. fmt(abs(nn[2]), 4)) or fmt(nn[1], 4)
            say(L, "Se elige como dominante " .. ((nd == 2) and "el par " or "el polo ") .. nm .. " porque es el mas cercano al eje imaginario (|Re| = " .. fmt(abs(dom[1]), 4) .. "): su modo e^(" .. ((abs(dom[1] + 1) < 1e-12) and "-" or fmt(dom[1], 4)) .. "t) decae mas lento (" .. TAU .. " = " .. fmt(1 / abs(dom[1]), 4) .. " s) que el del polo " .. nn_s .. " (" .. TAU .. " = " .. fmt(1 / abs(nn[1]), 4) .. " s), que se extingue " .. fmt(ratio, 3) .. " veces mas rapido; por eso es el que gobierna la respuesta.")
          end
          if ratio >= 5 then
            say(L, ((nd == 2) and ("Los polos dominantes son " .. fmt(dom[1], 4) .. " " .. PM .. " j" .. fmt(abs(dom[2]), 4)) or ("El polo dominante es " .. fmt(dom[1], 4))) .. ", ya que el siguiente polo esta " .. fmt(ratio, 3) .. " veces mas lejos del eje imaginario (" .. U.ge .. " 5); por ello es razonable aproximar por un modelo de orden " .. nd .. " que conserva la ganancia estatica.")
            local dpp = {dom}
            if nd == 2 then dpp[2] = {dom[1], -dom[2]} end
            local den = pfromroots(dpp, 1)
            if abs(T.d[1]) > 1e-12 then say(L, "El modelo reducido es Gr(s) = " .. fmt(T.n[1] / T.d[1] * den[1], 5) .. " / (" .. pstr(den) .. "), con Gr(0) = GLC(0) = " .. fmt(T.n[1] / T.d[1], 4) .. ".") end
          else
            say(L, "El siguiente polo esta solo " .. fmt(ratio, 3) .. " veces mas lejos que el dominante (< 5), por lo que NO es razonable despreciarlo.")
          end
        end
        if nd == 2 then
          local sg, w = -dom[1], abs(dom[2])
          local wnn = sqrt(sg * sg + w * w)
          local zz = sg / wnn
          say(L, "Con los polos dominantes (" .. ((#poles > 2) and "modelo reducido" or "segundo orden") .. "):")
          sayspecs(L, zz, wnn)
        else
          local tau1 = 1 / abs(dom[1])
          metric(L, "tau", TAU .. " dominante", tau1, "s"); metric(L, "ts2", "ts(2%)", 4 * tau1, "s"); metric(L, "Mp", "Mp", 0, "%")
          say(L, "La constante de tiempo dominante es " .. TAU .. " = 1/" .. fmt(abs(dom[1]), 4) .. " = " .. fmt(tau1, 4) .. " s.")
          say(L, "El tiempo de asentamiento al 2% es ts(2%) = 4" .. TAU .. " = " .. fmt(4 * tau1, 4) .. " s, y al 5% es ts(5%) = 3" .. TAU .. " = " .. fmt(3 * tau1, 4) .. " s.")
          say(L, "Como el polo dominante es real, no hay sobrepaso: Mp = 0 %.")
        end
        if m.Mp and #poles > 2 then
          say(L, "(Comprobacion con la respuesta exacta del sistema completo: " .. metricstext(m) .. ".)")
        end
        if m.yss then
          say(L, "El valor final ante el escalon es y(" .. U.inf .. ") = " .. fmt(m.yss) .. ".")
          metric(L, "yss", "y(" .. U.inf .. ")", m.yss)
          metric(L, "ess", "error escalon", A - m.yss)
          if #poles > 2 and m.Mp then metric(L, "Mpreal", "Mp real (simulado)", m.Mp, "%") end
        end
      end
    end
    chain(L, "4", 4, {Gc = v.Gc or "1", Gp = v.Gp, A = v.A or "1", k = v.k or "", eobj = "", senal = "1"}, "error estacionario de este lazo")
    chain(L, "5", 5, {P = "", Gc = v.Gc or "1", Gp = v.Gp, k = v.k or ""}, "Routh / estabilidad")
    chain(L, "6", 6, {Gc = v.Gc or "1", Gp = v.Gp, k = v.k or "", L = "0"}, "margenes / Bode de este lazo")
    summary(L, S)
  end
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
  do
    local tmax = horizon(gp)
    local ts, ys = {}, {}
    for i = 0, 300 do
      local t = tmax * i / 300
      ts[#ts + 1] = t
      ys[#ys + 1] = evalpf(groups, t)
    end
    local yss = nil
    if not diverge and ent == 1 and abs(G.d[1]) > 1e-12 then yss = A * G.n[1] / G.d[1] end
    local m = stepmetrics(ts, ys, yss)
    local ymax = -INF
    for i = 1, #ys do if ys[i] > ymax then ymax = ys[i] end end
    view(L, {kind = "time", title = "y(t)", ref = yss, band = yss ~= nil,
      curves = {{t = ts, y = ys, name = "y(t) " .. (diverge and "(diverge)" or "")}},
      marks = yss and {Mp = m.Mp, tp = m.tp, ymax = ymax, ts = m.ts2} or nil})
    local S = {txt}
    metric(L, "estable", "estabilidad", diverge and 0 or 1)
    if slow then metric(L, "dom", "Re polo dominante", slow[1]); metric(L, "ts2", "ts(2%) " .. U.approx .. " 4/|Re dom|", 4 / abs(slow[1]), "s") end
    if fast then metric(L, "fast", "Re polo mas veloz", fast[1]) end
    if yss then metric(L, "yss", "y(" .. U.inf .. ")", yss); if m.Mp then metric(L, "Mp", "Mp (medido)", m.Mp, "%") end end
    say(L, "La respuesta temporal ante " .. ((ent == 0) and "un impulso" or ((ent == 2) and "una rampa" or "un escalon")) .. " de magnitud " .. fmt(A) .. " es " .. txt .. ".")
    if diverge then
      say(L, "El sistema es INESTABLE: existen modos con parte real positiva (o nula) que no decaen, por lo que y(t) crece sin limite y nunca alcanza el estado estable; no tiene sentido hablar del modo mas veloz para llegar al estado estable.")
    else
      if fast and slow then
        say(L, "Todos los modos decaen, por lo que el sistema es estable. El modo mas veloz es el asociado al polo " .. fmtc(fast[1], fast[2], 4) .. " (mayor |Re|, se extingue en " .. U.approx .. " " .. fmt(4 / abs(fast[1]), 4) .. " s); el modo dominante es el del polo " .. fmtc(slow[1], slow[2], 4) .. ", que fija el tiempo de asentamiento ts " .. U.approx .. " " .. fmt(4 / abs(slow[1]), 4) .. " s.")
      end
      if yss then
        say(L, "El valor final es y(" .. U.inf .. ") = " .. fmt(yss) .. ".")
        if m.Mp then say(L, "Medido sobre y(t): el sobrepaso maximo es Mp = " .. fmt(m.Mp, 4) .. " %" .. (m.tp and (", el tiempo de peak es tp = " .. fmt(m.tp, 4) .. " s") or "") .. (m.tr and (", el tiempo de subida es tr = " .. fmt(m.tr, 4) .. " s") or "") .. (m.ts2 and (" y el tiempo de asentamiento es ts(2%) = " .. fmt(m.ts2, 4) .. " s") or "") .. ".") end
      end
    end
    if diverge then S[#S + 1] = "INESTABLE: y(t) no alcanza estado estable"
    else
      if fast then S[#S + 1] = "Modo mas veloz: polo " .. fmtc(fast[1], fast[2], 4) .. ";  dominante: " .. fmtc(slow[1], slow[2], 4) end
      if yss then S[#S + 1] = "y(" .. U.inf .. ") = " .. fmt(yss) .. ";  " .. metricstext(m) end
    end
    summary(L, S)
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
    local S = {}
    for _, kk in ipairs(sols) do S[#S + 1] = "k = " .. fmt(kk, 6) end
    if #S == 0 then S[1] = "sin solucion" end
    if sols[1] then metric(L, "ksol", "k requerido", sols[#sols]) end
    local pos = {}
    for _, kk in ipairs(sols) do if kk > 0 then pos[#pos + 1] = fmt(kk, 6) end end
    if #pos > 0 then
      say(L, "Para obtener un error en estado estable de " .. fmt(eobj) .. " ante " .. ({"un escalon", "una rampa", "una parabola"})[sig] .. " se requiere K = " .. table.concat(pos, " o ") .. " (verificando que el lazo cerrado sea estable).")
    end
    summary(L, S)
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
  metric(L, "N", "tipo N", e.N); metric(L, "Kp", "Kp", e.Kp); metric(L, "Kv", "Kv", e.Kv); metric(L, "Ka", "Ka", e.Ka)
  metric(L, "ep", "ep (escalon)", e.ep); metric(L, "ev", "ev (rampa)", e.ev); metric(L, "ea", "ea (parabola)", e.ea)
  metric(L, "estable", "estabilidad", ok and 1 or 0)
  say(L, "La funcion de lazo abierto es G0(s) = " .. rstr(G0) .. ", que tiene " .. math.max(e.N, 0) .. " integrador(es): el sistema es tipo " .. e.N .. ".")
  if not ok then say(L, "OJO: el lazo cerrado no es estable, por lo que el teorema del valor final no es aplicable.") end
  say(L, "Las constantes de error son Kp = " .. fmt(e.Kp) .. ", Kv = " .. fmt(e.Kv) .. " y Ka = " .. fmt(e.Ka) .. ".")
  say(L, "Por lo tanto, el error en estado estable ante un escalon de magnitud " .. fmt(A) .. " es ep = " .. fmt(e.ep) .. ", ante una rampa es ev = " .. fmt(e.ev) .. " y ante una parabola es ea = " .. fmt(e.ea) .. ".")
  chain(L, "3", 3, {Gc = v.Gc or "1", Gp = v.Gp, lazo = "1", A = v.A or "1", k = v.k or ""}, "polos/respuesta LC")
  chain(L, "5", 5, {P = "", Gc = v.Gc or "1", Gp = v.Gp, k = v.k or ""}, "Routh / estabilidad")
  chain(L, "6", 6, {Gc = v.Gc or "1", Gp = v.Gp, L = "0", k = v.k or ""}, "margenes / Bode")
  summary(L, {
    "Tipo N = " .. e.N .. ";  lazo cerrado " .. (ok and "estable" or "NO estable"),
    "Kp = " .. fmt(e.Kp) .. "  " .. U.arrow .. " ep (escalon) = " .. fmt(e.ep),
    "Kv = " .. fmt(e.Kv) .. "  " .. U.arrow .. " ev (rampa) = " .. fmt(e.ev),
    "Ka = " .. fmt(e.Ka) .. "  " .. U.arrow .. " ea (parabola) = " .. fmt(e.ea),
  })
  return L
end

--------------------------------------------------------------------------
-- 4) Routh-Hurwitz
--------------------------------------------------------------------------

--------------------------------------------------------------------------
-- Routh-Hurwitz simbolico en k (P(s) = P0(s) + k P1(s)), paso a paso
--------------------------------------------------------------------------
local function ksimp(r)
  local n, d = ptrim(pcopy(r.n)), ptrim(pcopy(r.d))
  if piszero(n) then return rat({0}, {1}) end
  if #d == 1 then
    if d[1] < 0 then n = pscale(n, -1); d = {-d[1]} end
    -- reducir por el MCD si todos los coeficientes son enteros
    local function isint(x) return abs(x - math.floor(x + 0.5)) < 1e-9 * math.max(1, abs(x)) and abs(x) < 1e9 end
    local allint = isint(d[1])
    for i = 1, #n do if not isint(n[i]) then allint = false end end
    if allint then
      local function gcd(a, b) a, b = abs(a), abs(b) while b > 0.5 do a, b = b, a % b end return a end
      local g = math.floor(d[1] + 0.5)
      for i = 1, #n do g = gcd(g, math.floor(abs(n[i]) + 0.5)) end
      if g > 1 then n = pscale(n, 1 / g); d = {d[1] / g} end
    else
      n = pscale(n, 1 / d[1]); d = {1}
    end
    if abs(d[1] - 1) < 1e-12 then d = {1} end
    return rat(n, d)
  end
  local q = rsimplify(rat(n, d))
  if #q.d == 1 then return ksimp(q) end
  return q
end

-- polinomio en k en orden ascendente: 140 - 10k
local function kpoly(p)
  local parts = {}
  for i = 1, #p do
    local c = p[i]
    if abs(c) > 1e-12 or #p == 1 then
      local e = i - 1
      local mag = abs(c)
      local cs = (e > 0 and abs(mag - 1) < 1e-12) and "" or fmt(mag)
      local term = (e == 0) and fmt(mag) or (cs .. "k" .. ((e > 1) and ("^" .. e) or ""))
      if #parts == 0 then parts[1] = (c < 0 and "-" or "") .. term
      else parts[#parts + 1] = (c < 0 and " - " or " + ") .. term end
    end
  end
  if #parts == 0 then return "0" end
  return table.concat(parts)
end

local function kstr(r)
  local n, d = r.n, r.d
  if #d == 1 then
    if #n == 1 then return fmt(n[1] / d[1]) end
    if abs(d[1] - 1) < 1e-12 then return kpoly(n) end
    return "(" .. kpoly(n) .. ")/" .. fmt(d[1])
  end
  return "(" .. kpoly(n) .. ")/(" .. kpoly(d) .. ")"
end

-- con parentesis si es suma (para usar como factor)
local function kfac(r)
  local t = kstr(r)
  if t:find(" [%+%-] ") then
    if t:sub(1, 1) == "(" and t:find("%)/") then return "[" .. t .. "]" end
    return "(" .. t .. ")"
  end
  if t:find("/") then return "[" .. t .. "]" end
  return t
end

local function kzero(r) return piszero(ptrim(pcopy(r.n))) end

-- intervalos (abiertos) donde r(k) > 0
local function kpositive(r)
  local pts = {}
  local function addroots(p)
    if #p > 1 then
      for _, q in ipairs(proots(p)) do if q[2] == 0 then pts[#pts + 1] = q[1] end end
    end
  end
  addroots(ptrim(pcopy(r.n))); addroots(ptrim(pcopy(r.d)))
  table.sort(pts)
  local u = {}
  for _, x in ipairs(pts) do if #u == 0 or abs(x - u[#u]) > 1e-9 * math.max(1, abs(x)) then u[#u + 1] = x end end
  local edges = {-INF}
  for _, x in ipairs(u) do edges[#edges + 1] = x end
  edges[#edges + 1] = INF
  local out = {}
  for i = 1, #edges - 1 do
    local a, b = edges[i], edges[i + 1]
    local t
    if a == -INF and b == INF then t = 0
    elseif a == -INF then t = b - math.max(1, abs(b))
    elseif b == INF then t = a + math.max(1, abs(a))
    else t = (a + b) / 2 end
    local val = peval(r.n, t) / peval(r.d, t)
    if val > 0 then
      if #out > 0 and out[#out][2] == a and false then out[#out][2] = b else out[#out + 1] = {a, b} end
    end
  end
  return out
end

local function kintersect(A, B)
  local out = {}
  for _, a in ipairs(A) do
    for _, b in ipairs(B) do
      local lo, hi = math.max(a[1], b[1]), math.min(a[2], b[2])
      if lo < hi - 1e-12 then out[#out + 1] = {lo, hi} end
    end
  end
  return out
end

local function ivtext(iv)
  if #iv == 0 then return "ningun k" end
  local t = {}
  for _, x in ipairs(iv) do
    local a, b = x[1], x[2]
    if a == -INF and b == INF then t[#t + 1] = "todo k"
    elseif a == -INF then t[#t + 1] = "k < " .. fmt(b, 6)
    elseif b == INF then t[#t + 1] = "k > " .. fmt(a, 6)
    else t[#t + 1] = fmt(a, 6) .. " < k < " .. fmt(b, 6) end
  end
  return table.concat(t, "  o  ")
end

-- agrega el desarrollo a L; devuelve los intervalos estables o nil si no aplica
local function routh_k_steps(L, P0, P1)
  local n = math.max(#P0, #P1) - 1
  if n < 1 then return nil end
  local coef = {}
  for i = n + 1, 1, -1 do
    coef[#coef + 1] = rat({P0[i] or 0, P1[i] or 0}, {1})
  end
  local w = math.floor(n / 2) + 1
  local rows = {{}, {}}
  for j = 1, w do
    rows[1][j] = coef[2 * j - 1] or rat({0}, {1})
    rows[2][j] = coef[2 * j] or rat({0}, {1})
  end
  local formulas = {}
  for i = 3, n + 1 do
    local A, B = rows[i - 2], rows[i - 1]
    if kzero(B[1]) then return nil end
    local C = {}
    for j = 1, w do
      local a2 = A[j + 1] or rat({0}, {1})
      local b2 = B[j + 1] or rat({0}, {1})
      C[j] = ksimp(rdiv(rsub(rmul(B[1], a2), rmul(A[1], b2)), B[1]))
    end
    rows[i] = C
    formulas[i] = "s^" .. (n - i + 1) .. ": [" .. kfac(B[1]) .. "\194\183" .. kfac(A[2] or rat({0}, {1})) .. " - " .. kfac(A[1]) .. "\194\183" .. kfac(B[2] or rat({0}, {1})) .. "] / " .. kfac(B[1]) .. " = " .. kstr(C[1])
  end
  sep(L, "Paso 2: Tabla de Routh en funcion de k")
  for i, r in ipairs(rows) do
    local cells = {}
    for j = 1, w do
      if j == 1 or not kzero(r[j]) then cells[#cells + 1] = kstr(r[j]) end
    end
    add(L, "s^" .. (n - i + 1) .. " | " .. table.concat(cells, "    "))
  end
  add(L, "Cada fila nueva: c1 = (b1" .. "\194\183" .. "a2 - a1" .. "\194\183" .. "b2)/b1 (a = fila de 2 arriba, b = fila de arriba)")
  for i = 3, n + 1 do add(L, "  " .. formulas[i]) end
  sep(L, "Paso 3: Condicion: toda la 1a columna > 0")
  local acc = {{-INF, INF}}
  for i, r in ipairs(rows) do
    local iv = kpositive(r[1])
    local txt = "s^" .. (n - i + 1) .. ":  " .. kstr(r[1]) .. " > 0"
    if #r[1].n == 1 and #r[1].d == 1 then
      txt = txt .. ((#iv > 0) and "  (siempre)" or "  (NUNCA: inestable)")
    else
      txt = txt .. "  " .. U.arrow .. "  " .. ivtext(iv)
    end
    add(L, txt)
    acc = kintersect(acc, iv)
  end
  sep(L, "Paso 4: Interseccion de las condiciones")
  add(L, "Estable para: " .. ivtext(acc))
  local phys = kintersect(acc, {{0, INF}})
  if #phys > 0 then add(L, "Con k > 0 (ganancia fisica): " .. ivtext(phys)) end
  -- limites: fila nula y polinomio auxiliar
  local lims = {}
  for _, x in ipairs(acc) do
    if x[1] ~= -INF then lims[#lims + 1] = x[1] end
    if x[2] ~= INF then lims[#lims + 1] = x[2] end
  end
  if #lims > 0 then
    sep(L, "Paso 5: Limites (oscilacion sostenida)")
    for _, kb in ipairs(lims) do
      local found = false
      for i = 2, #rows do
        local v1 = peval(rows[i][1].n, kb) / peval(rows[i][1].d, kb)
        if abs(v1) < 1e-7 * math.max(1, abs(kb)) and i >= 2 then
          local pw = n - (i - 2) -- potencia de la fila de arriba
          local up = rows[i - 1]
          local aux = {}
          for j = 1, w do
            local e = pw - 2 * (j - 1)
            if e >= 0 then aux[e + 1] = peval(up[j].n, kb) / peval(up[j].d, kb) end
          end
          for e = 1, pw + 1 do aux[e] = aux[e] or 0 end
          aux = ptrim(aux)
          if i == #rows then
            add(L, "k = " .. fmt(kb, 6) .. ": se anula el termino s^0 " .. U.arrow .. " raiz en s = 0")
          else
            add(L, "k = " .. fmt(kb, 6) .. ": se anula la fila s^" .. (pw - 1) .. " " .. U.arrow .. " pol. auxiliar (fila s^" .. pw .. "):")
            add(L, "  A(s) = " .. pstr(aux) .. " = 0  " .. U.arrow .. "  s = " .. rootsstr(proots(aux), 5))
          end
          found = true
          break
        end
      end
      if not found then add(L, "k = " .. fmt(kb, 6) .. ": cambia el signo de la 1a columna") end
    end
  end
  return acc
end

local function charpoly(v, kk)
  if not blank(v.P) then
    local r = tf(v, "P", kk)
    return r.n
  end
  local G0 = rmul(tf(v, "Gc", kk, "1"), tf(v, "Gp", kk))
  return padd(G0.d, G0.n)
end

-- lugar de raices para P(s) = P0 + k P1, k > 0
local function rlocus_view(P0, P1, bounds, kmark)
  local olp, olz = proots(P0), proots(P1)
  local kpos = {}
  for _, b in ipairs(bounds or {}) do if b > 0 then kpos[#kpos + 1] = b end end
  local kref = 1
  for _, b in ipairs(kpos) do if b > kref then kref = b end end
  local sc = 0
  for i = 1, #P0 do sc = math.max(sc, abs(P0[i])) end
  local s1 = 0
  for i = 1, #P1 do s1 = math.max(s1, abs(P1[i])) end
  local kmax = math.max(20 * kref, 100 * sc / math.max(s1, 1e-12))
  local kmin = 1e-4 * math.min(kref, sc / math.max(s1, 1e-12))
  local pts = {}
  local prev = nil
  local N = 160
  for i = 0, N do
    local kk = kmin * (kmax / kmin) ^ (i / N)
    local r = proots(padd(P0, pscale(P1, kk)), prev)
    prev = r
    for _, q in ipairs(r) do pts[#pts + 1] = {q[1], q[2]} end
  end
  -- extension de la vista
  local Rr = 1
  for _, q in ipairs(olp) do Rr = math.max(Rr, abs(q[1]), abs(q[2])) end
  for _, q in ipairs(olz) do Rr = math.max(Rr, abs(q[1]), abs(q[2])) end
  local cross = {}
  for _, kb in ipairs(kpos) do
    local r = proots(padd(P0, pscale(P1, kb)))
    for _, q in ipairs(r) do
      if abs(q[1]) < 1e-4 * math.max(1, abs(q[2])) and q[2] > 1e-6 then
        cross[#cross + 1] = {k = kb, w = q[2]}; Rr = math.max(Rr, q[2]); break
      end
    end
  end
  local lim = 2.2 * Rr
  local xmin, xmax, ymax = -0.3 * Rr, 0.3 * Rr, 0.5 * Rr
  for _, q in ipairs(pts) do
    if abs(q[1]) <= lim and abs(q[2]) <= lim then
      xmin = math.min(xmin, q[1]); xmax = math.max(xmax, q[1]); ymax = math.max(ymax, abs(q[2]))
    end
  end
  local pad = 0.1 * Rr
  local V = {kind = "rlocus", pts = pts, poles = olp, zeros = olz, cross = cross,
    xmin = xmin - pad, xmax = xmax + pad, ymax = ymax + pad, R = Rr}
  V.legend = {{"x rojo: polos LA (k=0)   o: ceros (k" .. U.arrow .. U.inf .. ")", "pp"}}
  for _, c in ipairs(cross) do
    V.legend[#V.legend + 1] = {"cruce eje j" .. U.w .. ": k = " .. fmt(c.k, 5) .. ", s = " .. U.pm .. "j" .. fmt(c.w, 4), "cz"}
  end
  if kmark then
    V.kpoles = proots(padd(P0, pscale(P1, kmark)))
    V.legend[#V.legend + 1] = {"azul: polos LC con k = " .. fmt(kmark), "pd"}
  end
  return V
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
    sep(L, "Paso 1: Polinomio caracteristico en funcion de k")
    add(L, "P(s) = den(Gc" .. "\194\183" .. "Gp) + num(Gc" .. "\194\183" .. "Gp) = 0")
    if lin then
      add(L, "P(s) = " .. polyk_str(P0, d1))
    else
      add(L, "P(s) depende de k en forma no lineal")
      add(L, "P(s)|k=1 = " .. pstr(P1))
    end
    sep(L, "Criterio de Routh-Hurwitz")
    add(L, "Se exige 1a columna de Routh > 0; los limites son los k donde")
    add(L, "un elemento se anula (polos sobre el eje j" .. W .. ").")
    local ivs, bounds = krange(function(kk) return charpoly(v, kk) end)
    local symok = false
    if lin then
      local okc, acc = pcall(routh_k_steps, L, P0, d1)
      symok = okc and acc ~= nil
    end
    if not symok then
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
    end -- not symok
    add(L, "(ingrese un valor en k para ver la tabla numerica)")
    local S = {}
    for _, iv in ipairs(ivs) do
      local a = (iv[1] == -INF) and ("-" .. U.inf) or fmt(iv[1], 6)
      local b = (iv[2] == INF) and U.inf or fmt(iv[2], 6)
      S[#S + 1] = "Estable para " .. a .. " < k < " .. b
    end
    if #S == 0 then S[1] = "Ningun k estabiliza el lazo" end
    if lin then say(L, "La ecuacion caracteristica de lazo cerrado es " .. polyk_str(P0, d1) .. " = 0.") end
    do
      local kmax, kmin = nil, nil
      for _, iv in ipairs(ivs) do if iv[2] > 0 then kmax = iv[2]; kmin = math.max(iv[1], 0) end end
      if kmax then metric(L, "kmax", "k maximo estable", kmax); metric(L, "kmin", "k minimo estable", kmin) end
      for _, kb in ipairs(bounds) do
        if kb > 0 then
          for _, q in ipairs(proots(charpoly(v, kb))) do
            if abs(q[1]) < 1e-4 * math.max(1, abs(q[2])) and q[2] > 1e-6 then metric(L, "wosc", U.w .. " de oscilacion", q[2], "rad/s"); break end
          end
        end
      end
    end
    do
      local phys = {}
      for _, iv in ipairs(ivs) do
        local a, b = math.max(iv[1], 0), iv[2]
        if a < b then phys[#phys + 1] = {a, b} end
      end
      local txts = {}
      for _, iv in ipairs(phys) do txts[#txts + 1] = fmt(iv[1], 6) .. " < K < " .. ((iv[2] == INF) and U.inf or fmt(iv[2], 6)) end
      if #txts > 0 then
        say(L, "Aplicando el criterio de Routh-Hurwitz (toda la primera columna positiva), el sistema en lazo cerrado es estable para " .. table.concat(txts, " o ") .. ".")
      else
        say(L, "Aplicando el criterio de Routh-Hurwitz, no existe K > 0 que haga estable el lazo cerrado.")
      end
      for _, kb in ipairs(bounds) do
        if kb > 0 then
          local r = proots(charpoly(v, kb))
          for _, q in ipairs(r) do
            if abs(q[1]) < 1e-4 * math.max(1, abs(q[2])) and q[2] > 1e-6 then
              say(L, "Para K = " .. fmt(kb, 6) .. " el sistema es marginalmente estable: oscila con frecuencia " .. W .. " = " .. fmt(q[2], 5) .. " rad/s (polos en s = " .. PM .. "j" .. fmt(q[2], 5) .. ").")
              break
            end
          end
        end
      end
    end
    if blank(v.P) then
      chain(L, "3", 3, {Gc = v.Gc or "k", Gp = v.Gp, lazo = "1", A = "1", k = v.k or ""}, "polos/respuesta LC (GLC en funcion de k si k vacio)")
      chain(L, "4", 4, {Gc = v.Gc or "k", Gp = v.Gp, A = "1", k = v.k or "", eobj = "", senal = "1"}, "error estacionario")
      chain(L, "6", 6, {Gc = v.Gc or "k", Gp = v.Gp, L = "0", k = v.k or ""}, "margenes / Bode")
    end

    if lin then
      local okv, V = pcall(rlocus_view, P0, d1, bounds, nil)
      if okv then view(L, V); S[#S + 1] = "tab: lugar de raices" end
    end
    summary(L, S)
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
  do
    local S = {"Cambios de signo en 1a columna: " .. ch, (stabilitytext(roots))}
    metric(L, "cambios", "cambios de signo", ch)
    metric(L, "estable", "estabilidad", isstable(roots) and 1 or 0)
    metric(L, "polos", "polos LC", rootsstr(roots, 4))
    say(L, "La ecuacion caracteristica es " .. pstr(P) .. " = 0.")
    say(L, "En la primera columna de la tabla de Routh hay " .. ch .. " cambio(s) de signo, por lo que hay " .. ch .. " polo(s) en el semiplano derecho: el sistema en lazo cerrado es " .. ((ch == 0 and isstable(roots)) and "ESTABLE." or ((ch == 0) and "criticamente estable (polos sobre el eje j" .. W .. ")." or "INESTABLE.")))
    if blank(v.P) then
      chain(L, "3", 3, {Gc = v.Gc or "k", Gp = v.Gp, lazo = "1", A = "1", k = v.k or ""}, "polos/respuesta LC (GLC en funcion de k si k vacio)")
      chain(L, "4", 4, {Gc = v.Gc or "k", Gp = v.Gp, A = "1", k = v.k or "", eobj = "", senal = "1"}, "error estacionario")
      chain(L, "6", 6, {Gc = v.Gc or "k", Gp = v.Gp, L = "0", k = v.k or ""}, "margenes / Bode")
    end

    if hask then
      local P0, P1, P2 = charpoly(v, 0), charpoly(v, 1), charpoly(v, 2)
      local d1 = psub(P1, P0)
      local chk = psub(psub(P2, P0), pscale(d1, 2))
      local lin = true
      for i = 1, #chk do if abs(chk[i]) > 1e-9 * math.max(1, abs(P2[i] or 0)) then lin = false end end
      if lin then
        local _, bounds = krange(function(kk) return charpoly(v, kk) end)
        local okv, V = pcall(rlocus_view, P0, d1, bounds, k)
        if okv then view(L, V) end
      end
    end
    summary(L, S)
  end
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
  do
    -- rango de frecuencias para el Bode
    local ws = {}
    if m.gc[1] then ws[#ws + 1] = m.gc[1].w end
    if m.pc[1] then ws[#ws + 1] = m.pc[1].w end
    for _, r in ipairs(o.p) do local a = sqrt(r[1] ^ 2 + r[2] ^ 2); if a > 1e-9 then ws[#ws + 1] = a end end
    for _, r in ipairs(o.z) do local a = sqrt(r[1] ^ 2 + r[2] ^ 2); if a > 1e-9 then ws[#ws + 1] = a end end
    if #ws == 0 then ws[1] = 1 end
    local lo, hi = INF, 0
    for _, x in ipairs(ws) do lo = math.min(lo, x); hi = math.max(hi, x) end
    local l0, l1 = math.floor(math.log10(lo) - 1), math.ceil(math.log10(hi) + 1)
    if l1 - l0 > 6 then l0 = l1 - 6 end
    local lw, mg_, ph_ = {}, {}, {}
    for i = 0, 200 do
      local lwv = l0 + (l1 - l0) * i / 200
      local mm, pp = freqat(o, 10 ^ lwv)
      lw[#lw + 1] = lwv; mg_[#mg_ + 1] = (mm > 0) and db(mm) or -200; ph_[#ph_ + 1] = pp
    end
    local V = {kind = "bode", lw = lw, mag = mg_, ph = ph_, legend = {}}
    if m.gc[1] then
      V.wcg, V.phcg, V.MF = m.gc[1].w, m.gc[1].ph, m.gc[1].MF
      V.legend[#V.legend + 1] = {"verde: " .. W .. "cg = " .. fmt(V.wcg, 4) .. " (|G0|=0 dB), MF = " .. fmt(V.MF, 4) .. DEG, "cz"}
    end
    if m.pc[1] then
      V.wcf, V.magcf = m.pc[1].w, db(m.pc[1].mag)
      V.legend[#V.legend + 1] = {"rojo: " .. W .. "cf = " .. fmt(V.wcf, 4) .. " (fase -180" .. DEG .. "), MG = " .. fmt(m.pc[1].MG, 4) .. " = " .. fmt(-V.magcf, 4) .. " dB", "pp"}
    end
    view(L, V)
    -- respuesta al escalon del lazo cerrado
    local S = {}
    S[#S + 1] = m.gc[1] and (W .. "cg = " .. fmt(m.gc[1].w) .. " rad/s,  MF = " .. fmt(m.gc[1].MF, 4) .. DEG) or ("MF = " .. U.inf)
    S[#S + 1] = m.pc[1] and (W .. "cf = " .. fmt(m.pc[1].w) .. " rad/s,  MG = " .. fmt(m.pc[1].MG, 4) .. " (" .. fmt(db(m.pc[1].MG), 4) .. " dB)") or ("MG = " .. U.inf)
    S[#S + 1] = (mf > 0 and mg > 1) and "Lazo cerrado ESTABLE" or "Lazo cerrado INESTABLE"
    if m.gc[1] then metric(L, "wcg", U.w .. "cg", m.gc[1].w, "rad/s"); metric(L, "MF", "MF", m.gc[1].MF, DEG) end
    if m.pc[1] then metric(L, "wcf", U.w .. "cf", m.pc[1].w, "rad/s"); metric(L, "MG", "MG", m.pc[1].MG); metric(L, "MGdB", "MG [dB]", db(m.pc[1].MG), "dB") end
    metric(L, "estable", "estabilidad", (mf > 0 and mg > 1) and 1 or 0)
    if m.gc[1] then say(L, "La frecuencia de cruce de ganancia (|G0(j" .. W .. ")| = 1) es " .. W .. "cg = " .. fmt(m.gc[1].w, 4) .. " rad/s; ahi la fase es " .. fmt(m.gc[1].ph, 4) .. DEG .. ", por lo que el margen de fase es MF = 180" .. DEG .. " + (" .. fmt(m.gc[1].ph, 4) .. DEG .. ") = " .. fmt(m.gc[1].MF, 4) .. DEG .. ".")
    else say(L, "La magnitud no cruza 1 (0 dB), por lo que el margen de fase es infinito.") end
    if m.pc[1] then say(L, "La frecuencia de cruce de fase (" .. U.ang .. "G0 = -180" .. DEG .. ") es " .. W .. "cf = " .. fmt(m.pc[1].w, 4) .. " rad/s; ahi |G0| = " .. fmt(m.pc[1].mag, 4) .. ", por lo que el margen de ganancia es MG = 1/" .. fmt(m.pc[1].mag, 4) .. " = " .. fmt(m.pc[1].MG, 4) .. " (" .. fmt(db(m.pc[1].MG), 4) .. " dB).")
    else say(L, "La fase nunca alcanza -180" .. DEG .. ", por lo que el margen de ganancia es infinito.") end
    say(L, (mf > 0 and mg > 1) and ("Como MF > 0" .. DEG .. " y MG > 1, el sistema en lazo cerrado es ESTABLE.") or ("Como MF " .. U.le .. " 0" .. DEG .. " o MG " .. U.le .. " 1, el sistema en lazo cerrado es INESTABLE."))
    if Ld == 0 then
      chain(L, "3", 3, {Gc = v.Gc or "1", Gp = v.Gp, lazo = "1", A = "1", k = v.k or ""}, "polos/respuesta LC")
      chain(L, "4", 4, {Gc = v.Gc or "1", Gp = v.Gp, A = "1", k = v.k or "", eobj = "", senal = "1"}, "error estacionario")
      chain(L, "5", 5, {P = "", Gc = v.Gc or "1", Gp = v.Gp, k = v.k or ""}, "Routh / estabilidad")
    end
    if Ld == 0 then
      local Vt, mt, st = timeview(closedloop(G0), 1, "Respuesta LC (escalon unitario)")
      view(L, Vt)
      if st and mt.yss then add(L, "Simulacion LC: y(" .. U.inf .. ") = " .. fmt(mt.yss) .. "; " .. metricstext(mt)) end
    elseif #Gc.n <= #Gc.d and mf > 0 and mg > 1 then
      local wc = m.gc[1] and m.gc[1].w or 1
      local tmax = 20 / wc
      local okc, ts, ys = pcall(simloop, Gc, Gp, Ld, 1, tmax, 3000)
      if okc then
        local yss = nil
        local e = sserrors(G0, 1)
        if e.N >= 1 then yss = 1 elseif e.N == 0 then yss = e.Kp / (1 + e.Kp) end
        local mt = stepmetrics(ts, ys, yss)
        local t2, y2 = decimate(ts, ys, 300)
        local ymax = -INF
        for i = 1, #ys do if ys[i] > ymax then ymax = ys[i] end end
        view(L, {kind = "time", title = "Respuesta LC con retardo", ref = yss, band = yss ~= nil,
          curves = {{t = t2, y = y2, name = "y(t) simulada (RK4)"}}, marks = {Mp = mt.Mp, tp = mt.tp, ymax = ymax, ts = mt.ts2}})
        add(L, "Simulacion LC con retardo: " .. metricstext(mt))
      end
    end
    summary(L, S)
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
  add(L, U.ang .. "G(pd) = " .. "\206\163" .. "ang. ceros - " .. "\206\163" .. "ang. polos; cada angulo es el del vector raiz" .. U.arrow .. "pd")
  add(L, "medido desde el eje real positivo (ver grafico: tecla tab)")
  local mag, ph = gainphase_at(zp, pp, lead, -sig, wd)
  local sump, sumz = 0, 0
  local dp, dz_ = {}, {}
  local function th(r) return atan2(wd - r[2], -sig - r[1]) * R2D end
  local function dist(r) return sqrt((-sig - r[1]) ^ 2 + (wd - r[2]) ^ 2) end
  -- grafico del plano s (pantalla 2, tecla tab)
  local GREEK = {"\206\177", "\206\178", "\206\179", "\206\180", "\206\181", "\206\183", "\206\186", "\206\187", "\206\188", "\206\189"}
  local plot = {kind = "splane", sig = sig, wd = wd, items = {}, legend = {}}
  view(L, plot)
  local gi = 0
  local pletters, zletters = {}, {}
  local function angitem(r, kind, desc)
    gi = gi + 1
    local g = GREEK[gi] or ("a" .. gi)
    local t = th(r)
    plot.items[#plot.items + 1] = {x = r[1], y = r[2], kind = kind, label = g, ang = t}
    plot.legend[#plot.legend + 1] = {g .. " = " .. desc .. " " .. fmtc(r[1], r[2], 3) .. U.arrow .. "pd: " .. fmt(t, 4) .. DEG, kind}
    return g, t
  end
  for _, r in ipairs(pp) do
    local g, t = angitem(r, "pp", "polo planta")
    sump = sump + t
    pletters[#pletters + 1] = g
    dp[#dp + 1] = dist(r)
    add(L, "  " .. g .. " (polo " .. fmtc(r[1], r[2], 4) .. ") = atan2(" .. fmt(wd - r[2], 4) .. ", " .. fmt(-sig - r[1], 4) .. ") = " .. fmt(t) .. DEG)
  end
  if useI then
    local g, t = angitem({0, 0}, "ci", "integrador (polo en 0)")
    sump = sump + t
    pletters[#pletters + 1] = g
    dp[#dp + 1] = sqrt(sig * sig + wd * wd)
    add(L, "  " .. g .. " (integrador, polo 0) = atan2(" .. fmt(wd, 4) .. ", " .. fmt(-sig, 4) .. ") = " .. fmt(t) .. DEG)
    ph = ph - t
    mag = mag / sqrt(sig * sig + wd * wd)
  end
  for _, r in ipairs(zp) do
    local g, t = angitem(r, "pz", "cero planta")
    sumz = sumz + t
    zletters[#zletters + 1] = g
    dz_[#dz_ + 1] = dist(r)
    add(L, "  " .. g .. " (cero " .. fmtc(r[1], r[2], 4) .. ") = " .. fmt(t) .. DEG)
  end
  local zs = (#zletters > 0) and ("(" .. table.concat(zletters, "+") .. ")") or "0"
  local ps = (#pletters > 0) and ("(" .. table.concat(pletters, "+") .. ")") or "0"
  add(L, U.ang .. "G(pd) = " .. zs .. " - " .. ps .. " = " .. fmt(sumz) .. DEG .. " - " .. fmt(sump) .. DEG .. ((lead < 0) and " - 180" .. DEG .. " (K<0)" or "") .. " = " .. fmt(ph) .. DEG)
  plot.legend[#plot.legend + 1] = {U.ang .. "G(pd) = " .. zs .. " - " .. ps .. " = " .. fmt(ph, 4) .. DEG, "txt"}
  local phi = -180 - ph
  while phi > 180 do phi = phi - 360 end
  while phi <= -180 do phi = phi + 360 end
  add(L, "Deficiencia: " .. PHI .. "d = -180" .. DEG .. " - (" .. fmt(ph) .. DEG .. ") = " .. fmt(phi) .. DEG)
  add(L, "(" .. PHI .. "d = fase que deben aportar los ceros del controlador)")
  plot.legend[#plot.legend + 1] = {PHI .. "d = -180" .. DEG .. " - " .. U.ang .. "G(pd) = " .. fmt(phi, 4) .. DEG .. " (lo aporta el cero verde)", "txt"}
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
    plot.items[#plot.items + 1] = {x = -zc, y = 0, kind = "cz", label = PHI, ang = phz}
    plot.legend[#plot.legend + 1] = {PHI .. " = cero controlador -" .. fmt(zc, 4) .. U.arrow .. "pd: " .. fmt(phz, 4) .. DEG .. (nz == 2 and " (x2)" or ""), "cz"}
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
    plot.items[#plot.items + 1] = {x = q[1], y = q[2], kind = "lc"}
    add(L, "Polo adicional " .. fmtc(q[1], q[2], 4) .. ": |Re|/" .. SIG .. " = " .. fmt(abs(q[1]) / sig, 3) .. ((abs(q[1]) / sig >= 5) and " (" .. U.ge .. "5, despreciable)" or " (<5, influye)"))
  end
  add(L, (stabilitytext(r)))
  sep(L, "Paso 10: Error estacionario (A=1)")
  local e = sserrors(rsimplify(G0), 1)
  errlines(L, e, 1)
  if e.ep > 0 and e.ep < INF then
    add(L, "y(" .. U.inf .. ") = 1 - ep = " .. fmt(1 - e.ep) .. " (escalon unitario)")
  end
  do
    local T = rat(pcopy(G0.n), padd(G0.d, G0.n))
    local V, m, st = timeview(T, 1, "Respuesta LC con Gc (escalon unitario)")
    view(L, V)
    sep(L, "Paso 11: Comprobacion por simulacion")
    add(L, "Respuesta exacta del lazo cerrado (incluye el efecto de los ceros):")
    if st and m.yss then
      add(L, "y(" .. U.inf .. ") = " .. fmt(m.yss) .. ";  " .. metricstext(m))
      if Mp then add(L, "Mp real " .. fmt(m.Mp, 4) .. " % vs pedido " .. fmt(Mp) .. " %" .. ((m.Mp > Mp * 1.15 + 0.5) and "  (mayor: efecto del cero)" or "  " .. U.arrow .. " cumple aprox.")) end
      if ts and m.ts2 then
        local tsr = (crit == 5) and m.ts5 or m.ts2
        if tsr then add(L, "ts real " .. fmt(tsr, 4) .. " s vs pedido " .. fmt(ts) .. " s") end
      end
    end
    local gcs = fmt(Kc)
    if nz == 1 then gcs = gcs .. "(s + " .. fmt(zc) .. ")" elseif nz == 2 then gcs = gcs .. "(s + " .. fmt(zc) .. ")" .. U.sq end
    if useI then gcs = gcs .. "/s" end
    local S = {"pd = -" .. fmt(sig) .. " " .. PM .. " j" .. fmt(wd) .. ";  " .. PHI .. "d = " .. fmt(phi, 4) .. DEG,
      "Estructura: " .. names[nz] .. ";  Gc(s) = " .. gcs}
    if useI and nz == 1 then S[#S + 1] = "Kp = " .. fmt(Kc) .. ",  Ti = " .. fmt(1 / zc) .. " s"
    elseif useI and nz == 2 then S[#S + 1] = "Kp = " .. fmt(2 * Kc * zc) .. ",  Ti = " .. fmt(2 / zc) .. " s,  Td = " .. fmt(1 / (2 * zc)) .. " s"
    elseif (not useI) and nz == 1 then S[#S + 1] = "Kp = " .. fmt(Kc * zc) .. ",  Td = " .. fmt(1 / zc) .. " s" end
    S[#S + 1] = "Polos LC: " .. rootsstr(r, 4)
    S[#S + 1] = "ep = " .. fmt(e.ep) .. ",  ev = " .. fmt(e.ev)
    if st and m.yss then S[#S + 1] = "Simulacion: " .. metricstext(m) end
    S[#S + 1] = "tab: plano s y respuesta y(t)"
    L.res = {phi = phi, a = zc or 0, Kc = Kc}
    do
      local wnn = sqrt(sig * sig + wd * wd)
      metric(L, "zeta", U.zeta, sig / wnn); metric(L, "wn", U.wn, wnn, "rad/s")
      metric(L, "sigma", U.sigma .. " (|Re pd|)", sig); metric(L, "wd", U.w .. "d (Im pd)", wd)
      metric(L, "phi", PHI .. "d", phi, DEG)
      metric(L, "estr", "estructura", names[nz])
      if zc then metric(L, "a", "cero a", zc) end
      metric(L, "Kc", "Kc", Kc)
      metric(L, "ep", "ep (escalon)", e.ep)
      if st and m.Mp then metric(L, "Mpreal", "Mp real", m.Mp, "%") end
      if st and m.ts2 then metric(L, "tsreal", "ts(2%) real", m.ts2, "s") end
    end
    do
      local wnn = sqrt(sig * sig + wd * wd)
      local zz = sig / wnn
      say(L, "De Mp y ts se obtiene " .. Z .. " = " .. fmt(zz, 4) .. " y " .. WN .. " = " .. fmt(wnn, 4) .. " rad/s; los polos deseados de lazo cerrado son pd = -" .. fmt(sig, 4) .. " " .. PM .. " j" .. fmt(wd, 4) .. ".")
    end
    if useI then say(L, "Como se exige error nulo ante escalon y la planta es tipo 0, se agrega accion integral (polo en el origen).")
    else say(L, "Como no se exige seguimiento exacto en regimen permanente, no se requiere accion integral.") end
    say(L, "La condicion de angulo en pd entrega una deficiencia angular " .. PHI .. "d = " .. fmt(phi, 4) .. DEG .. "; por lo tanto, la estructura minima es un controlador " .. names[nz] .. ((nz == 1) and " (un cero aporta entre 0" .. DEG .. " y 90" .. DEG .. ")." or ((nz == 2) and " (se requieren dos ceros, " .. PHI .. "d/2 cada uno)." or ".")))
    if nz > 0 then say(L, "El cero del controlador se ubica en s = -" .. fmt(zc, 5) .. " y, por la condicion de magnitud, Kc = " .. fmt(Kc, 5) .. ".") end
    say(L, "El controlador obtenido es Gc(s) = " .. gcs .. ".")
    say(L, "Se verifica que los polos de lazo cerrado son " .. rootsstr(r, 4) .. ", que incluyen a pd; el sistema es " .. ((isstable(r)) and "estable." or "inestable."))
    if st and m.yss then
      if m.Mp then say(L, "Con el controlador, el sobrepaso maximo real es Mp = " .. fmt(m.Mp, 4) .. " %" .. (Mp and (" (pedido: " .. fmt(Mp) .. " %)") or "") .. ".") end
      if m.tp then say(L, "El tiempo de peak real es tp = " .. fmt(m.tp, 4) .. " s.") end
      if m.tr then say(L, "El tiempo de subida real es tr = " .. fmt(m.tr, 4) .. " s.") end
      if m.ts2 then say(L, "El tiempo de asentamiento real es ts(2%) = " .. fmt(m.ts2, 4) .. " s" .. (ts and (" (pedido: " .. fmt(ts) .. " s)") or "") .. ".") end
    end
    say(L, "El error en estado estable ante escalon unitario es ep = " .. fmt(e.ep, 4) .. ((e.ep == 0) and " (seguimiento exacto)." or (" (la salida se estabiliza en " .. fmt(1 - e.ep, 4) .. ").")))
    local gtxt = rexact(Gc)
    chain(L, "3", 3, {Gc = gtxt, Gp = v.Gp, lazo = "1", A = "1", k = ""}, "polos y respuesta LC con este Gc")
    chain(L, "6", 6, {Gc = gtxt, Gp = v.Gp, L = "0", k = ""}, "margenes / Bode con este Gc")
    summary(L, S)
  end
  return L
end

--------------------------------------------------------------------------
-- 7) Sintonia empirica (ZN, rele, Ho-Hang-Cao)
--------------------------------------------------------------------------
-- PID con derivada filtrada (Td/10) para poder simular
local function pidtf_sim(Kp, Ti, Td)
  if not Td or Td == 0 then
    if not Ti then return rconst(Kp) end
    return rat({Kp, Kp * Ti}, {0, Ti})
  end
  local a = Td / 10
  -- Kp(1 + 1/(Ti s) + Td s/(a s + 1))
  local base = rat({Kp, Kp * Ti}, {0, Ti})
  return radd(base, rat({0, Kp * Td}, {1, a}))
end

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
  if haveModel and #verify > 0 then
    sep(L, "Simulacion LC (escalon unitario, retardo exacto)")
    add(L, "(PID con derivada filtrada Td/10 para simular)")
    local curves, S = {}, {}
    local wc = 1
    do
      local G0 = rmul(pidtf(verify[1][2], verify[1][3], verify[1][4]), Gp)
      local mm = margins(freqobj(G0, Ld))
      if mm.gc[1] then wc = mm.gc[1].w end
    end
    local tmax = math.max(25 / wc, 6 * Ld)
    local function simc(name, Kp, Ti, Td)
      local okc, ts, ys = pcall(simloop, pidtf_sim(Kp, Ti, Td), Gp, Ld, 1, tmax, 1500)
      if not okc then return end
      local mt = stepmetrics(ts, ys, 1)
      local t2, y2 = decimate(ts, ys, 250)
      local div = abs(ys[#ys] - 1) > 0.5
      add(L, name .. ": " .. (div and "no converge (inestable)" or metricstext(mt)))
      if #curves < 3 then curves[#curves + 1] = {t = t2, y = y2, name = name .. (div and " (inestable)" or (", Mp=" .. fmt(mt.Mp or 0, 3) .. "%"))} end
      return mt, ys
    end
    local first
    for i, c in ipairs(verify) do
      local mt = simc(c[1], c[2], c[3], c[4])
      if i == 1 then first = mt end
    end
    if fK and hhc then simc("HHC Kp x " .. fmt(fK), fK * hhc[2], hhc[3], hhc[4]) end
    if #curves > 0 then
      local ymax = -INF
      for _, yv in ipairs(curves[1].y) do if yv > ymax then ymax = yv end end
      view(L, {kind = "time", title = "Lazo cerrado con retardo", ref = 1, band = true, curves = curves,
        marks = first and {Mp = first.Mp, tp = first.tp, ymax = ymax, ts = first.ts2} or nil})
    end
  end
  do
    local S = {}
    if haveModel then S[#S + 1] = "ZN-LA PI: Kp = " .. fmt(0.9 * tau / (K * Ld)) .. ", Ti = " .. fmt(3.33 * Ld) .. " s" end
    if Ku and Pu then S[#S + 1] = "ZN-LC PI: Kp = " .. fmt(0.45 * Ku) .. ", Ti = " .. fmt(Pu / 1.2) .. " s  (Ku=" .. fmt(Ku, 4) .. ", Pu=" .. fmt(Pu, 4) .. ")" end
    if hhc then
      S[#S + 1] = "HHC: Kp = " .. fmt(hhc[2]) .. ", Ti = " .. fmt(hhc[3]) .. " s" .. (hhc[4] and (", Td = " .. fmt(hhc[4])) or "")
      local G0 = rmul(pidtf(hhc[2], hhc[3], hhc[4]), Gp)
      local mm = margins(freqobj(G0, Ld))
      if mm.gc[1] then S[#S + 1] = "HHC verif.: " .. W .. "cg = " .. fmt(mm.gc[1].w, 4) .. ", MF = " .. fmt(mm.gc[1].MF, 4) .. DEG end
      if mm.pc[1] then S[#S + 1] = "            " .. W .. "cf = " .. fmt(mm.pc[1].w, 4) .. ", MG = " .. fmt(mm.pc[1].MG, 4) end
      local gtxt = string.format("%.7g*(1+1/(%.7g*s))", hhc[2], hhc[3])
      if hhc[4] then gtxt = string.format("%.7g*(1+1/(%.7g*s)+%.7g*s)", hhc[2], hhc[3], hhc[4]) end
      local gp = string.format("%.7g/(%.7g*s+1)", K, tau)
      if tau2 then gp = string.format("%.7g/((%.7g*s+1)*(%.7g*s+1))", K, tau, tau2) end
      chain(L, "6", 6, {Gc = gtxt, Gp = gp, L = string.format("%.7g", Ld), k = ""}, "Bode / margenes del lazo HHC")
    end
    if haveModel then metric(L, "znKp", "ZN-LA PI Kp", 0.9 * tau / (K * Ld)); metric(L, "znTi", "ZN-LA PI Ti", 3.33 * Ld, "s") end
    if Ku then metric(L, "Ku", "Ku", Ku) end
    if hhc then
      metric(L, "hKp", "HHC Kp", hhc[2]); metric(L, "hTi", "HHC Ti", hhc[3], "s")
      local G0m = rmul(pidtf(hhc[2], hhc[3], hhc[4]), Gp)
      local mm0 = margins(freqobj(G0m, Ld))
      if mm0.gc[1] then metric(L, "wcg", U.w .. "cg (HHC)", mm0.gc[1].w, "rad/s"); metric(L, "MF", "MF (HHC)", mm0.gc[1].MF, DEG) end
      if mm0.pc[1] then metric(L, "MG", "MG (HHC)", mm0.pc[1].MG) end
    end
    if haveModel then say(L, "Por Ziegler-Nichols (curva de reaccion), el controlador PI es Kp = 0.9" .. TAU .. "/(KL) = " .. fmt(0.9 * tau / (K * Ld), 4) .. " y Ti = 3.33L = " .. fmt(3.33 * Ld, 4) .. " s.") end
    if Ku and Pu then say(L, "Por Ziegler-Nichols en lazo cerrado (Ku = " .. fmt(Ku, 4) .. ", Pu = " .. fmt(Pu, 4) .. " s), el controlador PI es Kp = 0.45Ku = " .. fmt(0.45 * Ku, 4) .. " y Ti = Pu/1.2 = " .. fmt(Pu / 1.2, 4) .. " s.") end
    if hhc then
      say(L, "Por Ho-Hang-Cao (MG = " .. fmt(MG) .. ", MF = " .. fmt(MF) .. DEG .. "), el controlador es Kp = " .. fmt(hhc[2], 4) .. " y Ti = " .. fmt(hhc[3], 4) .. " s" .. (hhc[4] and (", Td = " .. fmt(hhc[4], 4) .. " s") or "") .. ".")
      local G0 = rmul(pidtf(hhc[2], hhc[3], hhc[4]), Gp)
      local mm = margins(freqobj(G0, Ld))
      if mm.gc[1] and mm.pc[1] then
        say(L, "Verificando con el retardo exacto: " .. W .. "cg = " .. fmt(mm.gc[1].w, 4) .. " rad/s con MF = " .. fmt(mm.gc[1].MF, 4) .. DEG .. ", y " .. W .. "cf = " .. fmt(mm.pc[1].w, 4) .. " rad/s con MG = " .. fmt(mm.pc[1].MG, 4) .. "; los margenes coinciden con los de diseno y el lazo cerrado es estable.")
      end
      if fK then
        local G1 = rmul(pidtf(fK * hhc[2], hhc[3], hhc[4]), Gp)
        local m1 = margins(freqobj(G1, Ld))
        if mm.gc[1] and m1.gc[1] and mm.pc[1] and m1.pc[1] then
          say(L, "Al multiplicar Kp por " .. fmt(fK) .. " (Ti fijo), " .. W .. "cg sube de " .. fmt(mm.gc[1].w, 4) .. " a " .. fmt(m1.gc[1].w, 4) .. " rad/s (respuesta mas rapida), pero el MF baja de " .. fmt(mm.gc[1].MF, 4) .. DEG .. " a " .. fmt(m1.gc[1].MF, 4) .. DEG .. " y el MG de " .. fmt(mm.pc[1].MG, 4) .. " a " .. fmt(m1.pc[1].MG, 4) .. ": disminuye la estabilidad relativa (mas oscilacion y sobrepaso)" .. ((m1.gc[1].MF <= 0 or m1.pc[1].MG <= 1) and ", llegando a ser INESTABLE." or ", aunque sigue estable."))
        end
      end
    end
    if #S > 0 then
      if L.views then S[#S + 1] = "tab: respuesta simulada" end
      summary(L, S)
    end
  end
  return L
end


--------------------------------------------------------------------------
-- 9) Identificacion de un modelo FOPDT desde la curva de reaccion
--------------------------------------------------------------------------
local function m_ident(v)
  local L = {}
  local du = num(v, "du") or 1
  local y0 = num(v, "y0") or 0
  local yf = req(num(v, "yf"), "y final")
  local t28, t63 = num(v, "t28"), num(v, "t63")
  local Lt, taut = num(v, "Lt"), num(v, "taut")
  local dy = yf - y0
  if du == 0 then error("du no puede ser 0") end
  sep(L, "Paso 1: Ganancia")
  local K = dy / du
  add(L, "K = (y final - y inicial)/" .. "\206\148" .. "u = (" .. fmt(yf) .. " - " .. fmt(y0) .. ")/" .. fmt(du) .. " = " .. fmt(K))
  local tau, Ld
  if Lt and taut then
    sep(L, "Paso 2: Metodo de la tangente (datos leidos)")
    tau, Ld = taut, Lt
    add(L, "L = " .. fmt(Ld) .. " s (corte de la tangente con y inicial)")
    add(L, TAU .. " = " .. fmt(tau) .. " s (desde L hasta el corte con y final)")
  else
    req(t28, "t28"); req(t63, "t63")
    sep(L, "Paso 2: Metodo de dos puntos (Smith)")
    add(L, "y(t28) = y0 + 0.283" .. "\206\148" .. "y = " .. fmt(y0 + 0.283 * dy) .. ",  y(t63) = y0 + 0.632" .. "\206\148" .. "y = " .. fmt(y0 + 0.632 * dy))
    tau = 1.5 * (t63 - t28)
    Ld = t63 - tau
    add(L, TAU .. " = 1.5(t63 - t28) = 1.5(" .. fmt(t63) .. " - " .. fmt(t28) .. ") = " .. fmt(tau) .. " s")
    add(L, "L = t63 - " .. TAU .. " = " .. fmt(t63) .. " - " .. fmt(tau) .. " = " .. fmt(Ld) .. " s")
    if Ld < 0 then add(L, "OJO: L < 0 (revise los tiempos leidos)") end
  end
  sep(L, "Modelo")
  add(L, "Gp(s) = " .. fmt(K) .. "e^(-" .. fmt(Ld) .. "s)/(" .. fmt(tau) .. "s + 1)")
  add(L, "Retardo normalizado L/" .. TAU .. " = " .. fmt(Ld / tau, 4))
  sep(L, "Sintonia ZN (curva de reaccion)")
  if Ld > 0 then
    local P = tau / (K * Ld)
    add(L, "PI : Kp = 0.9" .. TAU .. "/(KL) = " .. fmt(0.9 * P) .. "   Ti = 3.33L = " .. fmt(3.33 * Ld))
    add(L, "PID: Kp = 1.2" .. TAU .. "/(KL) = " .. fmt(1.2 * P) .. "  Ti = 2L = " .. fmt(2 * Ld) .. "  Td = 0.5L = " .. fmt(0.5 * Ld))
  end
  -- curva del modelo
  local ts, ys = {}, {}
  local tmax = Ld + 5 * tau
  for i = 0, 200 do
    local t = tmax * i / 200
    ts[#ts + 1] = t
    ys[#ys + 1] = y0 + ((t > Ld) and dy * (1 - exp(-(t - Ld) / tau)) or 0)
  end
  view(L, {kind = "time", title = "Modelo identificado", ref = yf, band = false,
    curves = {{t = ts, y = ys, name = "modelo FOPDT"}},
    legend = {"K = " .. fmt(K, 4) .. ",  " .. TAU .. " = " .. fmt(tau, 4) .. " s,  L = " .. fmt(Ld, 4) .. " s"}})
  chain(L, "8", 8, {K = string.format("%.7g", K), tau = string.format("%.7g", tau), L = string.format("%.7g", Ld), tau2 = ""}, "llevar el modelo a Sintonia (ZN / HHC)")
  metric(L, "K", "K", K); metric(L, "tau", TAU, tau, "s"); metric(L, "L", "L", Ld, "s"); metric(L, "Lt", "L/" .. TAU, Ld / tau)
  say(L, "A partir de la curva de reaccion, la ganancia es K = " .. fmt(K, 4) .. ", la constante de tiempo es " .. TAU .. " = " .. fmt(tau, 4) .. " s y el retardo es L = " .. fmt(Ld, 4) .. " s.")
  say(L, "El modelo identificado es Gp(s) = " .. fmt(K, 4) .. "e^(-" .. fmt(Ld, 4) .. "s)/(" .. fmt(tau, 4) .. "s + 1).")
  summary(L, {"K = " .. fmt(K) .. ",  " .. TAU .. " = " .. fmt(tau) .. " s,  L = " .. fmt(Ld) .. " s",
    "Gp = " .. fmt(K) .. "e^(-" .. fmt(Ld) .. "s)/(" .. fmt(tau) .. "s+1)"})
  return L
end

--------------------------------------------------------------------------
-- 10) Entrenamiento: ejercicios aleatorios con correccion
--------------------------------------------------------------------------
local function rng(seed)
  local x = (seed * 7919 + 104729) % 2147483647
  return function(t)
    x = (x * 48271) % 2147483647
    return t[(x % #t) + 1]
  end
end

local TRAIN = {}

TRAIN[1] = function(R)
  local Mp, ts = R({5, 9, 10, 15, 16.3, 20, 25, 30}), R({1, 2, 3, 4, 5, 8})
  local z = zeta_from_Mp(Mp)
  local sg = 4 / ts
  local wn = sg / z
  return {txt = {"Un sistema de 2o orden debe tener Mp = " .. fmt(Mp) .. " % y ts(2%) = " .. fmt(ts) .. " s.",
    "Determine " .. U.zeta .. ", " .. U.wn .. " y " .. U.w .. "d de los polos deseados."},
    ask = {U.zeta, U.wn .. " [rad/s]", U.w .. "d [rad/s]"}, ans = {z, wn, wn * sqrt(1 - z * z)},
    mod = 1, v = {Mp = tostring(Mp), ts = tostring(ts), crit = "2"}}
end

local function placeprob(R, enulo)
  for _ = 1, 30 do
    local p1, p2, b = R({1, 2, 3}), R({4, 5, 6, 8, 10}), R({2, 3, 4, 5, 10})
    local ts, Mp = R({1, 1.5, 2, 3, 4, 5}), R({8, 10, 15, 16.3, 20})
    local gp = string.format("%g/((s+%g)(s+%g))", b, p1, p2)
    local ok, L = pcall(MODULES[7].run, {Gp = gp, ts = tostring(ts), Mp = tostring(Mp), crit = "2", enulo = tostring(enulo), estr = "0"})
    if ok and L.res and L.res.phi > 2 and L.res.phi < 170 then
      return {txt = {"Gp(s) = " .. gp .. ", retroalimentacion unitaria.",
        "Se pide ts(2%) = " .. fmt(ts) .. " s y Mp = " .. fmt(Mp) .. " %" .. (enulo == 1 and ", con error nulo ante escalon." or ", sin exigir error nulo."),
        "Disene el controlador de estructura minima por localizacion de polos.",
        "r1 = deficiencia " .. U.phi .. "d [" .. U.deg .. "],  r2 = cero a (s = -a; 0 si no hay),  r3 = Kc"},
        ask = {U.phi .. "d [" .. U.deg .. "]", "cero a", "Kc"}, ans = {L.res.phi, L.res.a, L.res.Kc},
        mod = 7, v = {Gp = gp, ts = tostring(ts), Mp = tostring(Mp), crit = "2", enulo = tostring(enulo), estr = "0"}}
    end
  end
  error("no se pudo generar el ejercicio, cambie el numero")
end
TRAIN[2] = function(R) return placeprob(R, 0) end
TRAIN[3] = function(R) return placeprob(R, 1) end

TRAIN[4] = function(R)
  local a, b, c = R({1, 2, 3}), R({4, 5, 6, 8, 10}), R({1, 2, 5, 10})
  local gp = string.format("%g/(s(s+%g)(s+%g))", c, a, b)
  local km = (a + b) * a * b / c
  return {txt = {"Con Gc(s) = K y Gp(s) = " .. gp .. " (retro unitaria):",
    "Determine el rango de K > 0 estable (Routh) y la frecuencia de oscilacion en el limite."},
    ask = {"K maximo", U.w .. " cruce [rad/s]", "(vacio)"}, ans = {km, sqrt(a * b)},
    mod = 5, v = {P = "", Gc = "k", Gp = gp, k = ""}}
end

TRAIN[5] = function(R)
  local kind = R({0, 1})
  local a, b, K = R({1, 2, 4, 5}), R({3, 6, 10}), R({2, 5, 8, 10, 20})
  if kind == 1 then
    local gp = string.format("%g/(s(s+%g))", K, a)
    return {txt = {"G0(s) = " .. gp .. " (retro unitaria).", "Determine el tipo, Kv y el error ante rampa unitaria."},
      ask = {"tipo N", "Kv", "ev"}, ans = {1, K / a, a / K}, mod = 4, v = {Gc = "1", Gp = gp, A = "1", k = ""}}
  end
  local gp = string.format("%g/((s+%g)(s+%g))", K, a, b)
  local kp = K / (a * b)
  return {txt = {"G0(s) = " .. gp .. " (retro unitaria).", "Determine el tipo, Kp y el error ante escalon unitario."},
    ask = {"tipo N", "Kp", "ep"}, ans = {0, kp, 1 / (1 + kp)}, mod = 4, v = {Gc = "1", Gp = gp, A = "1", k = ""}}
end

TRAIN[6] = function(R)
  for _ = 1, 30 do
    local a, b, c, K = R({1, 2}), R({3, 4, 5}), R({6, 8, 10}), R({20, 40, 60, 100, 150})
    local gp = string.format("%g/((s+%g)(s+%g)(s+%g))", K, a, b, c)
    local G = parse(gp)
    local m = margins(freqobj(G, 0))
    if m.gc[1] and m.pc[1] and m.gc[1].MF > 5 then
      return {txt = {"G0(s) = " .. gp .. ".", "Determine " .. U.w .. "cg, el margen de fase y el margen de ganancia (lineal)."},
        ask = {U.w .. "cg [rad/s]", "MF [" .. U.deg .. "]", "MG"}, ans = {m.gc[1].w, m.gc[1].MF, m.pc[1].MG},
        mod = 6, v = {Gc = "1", Gp = gp, L = "0", w = "", A = "", k = ""}}
    end
  end
  error("no se pudo generar el ejercicio, cambie el numero")
end

TRAIN[7] = function(R)
  local K, tau, Ld = R({1, 1.5, 2, 3}), R({5, 8, 10, 12, 20}), R({0.5, 1, 1.5, 2})
  local MG, MF = R({2, 3, 4}), R({45, 50, 60})
  local mf = MF * D2R
  local wp = (MG * (mf - pi / 2) + (pi / 2) * MG * MG) / ((MG * MG - 1) * Ld)
  local Kp = wp * tau / (MG * K)
  local Ti = 1 / (2 * wp - 4 * wp * wp * Ld / pi + 1 / tau)
  return {txt = {"Gp(s) = " .. fmt(K) .. "e^(-" .. fmt(Ld) .. "s)/(" .. fmt(tau) .. "s+1).",
    "Disene un PI por Ho-Hang-Cao con MG = " .. MG .. " y MF = " .. MF .. DEG .. "."},
    ask = {U.w .. "p [rad/s]", "Kp", "Ti [s]"}, ans = {wp, Kp, Ti},
    mod = 8, v = {K = tostring(K), tau = tostring(tau), L = tostring(Ld), MG = tostring(MG), MF = tostring(MF)}}
end

TRAIN[8] = function(R)
  local K, tau, Ld = R({1, 1.5, 2, 3}), R({10, 20, 50, 100, 200}), R({2, 5, 10, 15})
  local Kp, Ti = 0.9 * tau / (K * Ld), 3.33 * Ld
  local G0 = rmul(rat({Kp, Kp * Ti}, {0, Ti}), rat({K}, {1, tau}))
  local m = margins(freqobj(G0, Ld))
  return {txt = {"Gp(s) = " .. fmt(K) .. "e^(-" .. fmt(Ld) .. "s)/(" .. fmt(tau) .. "s+1).",
    "a) Estime un PI por Ziegler-Nichols (lazo abierto).  b) Margen de fase del sistema controlado."},
    ask = {"Kp", "Ti [s]", "MF [" .. U.deg .. "]"}, ans = {Kp, Ti, m.gc[1] and m.gc[1].MF or INF},
    mod = 8, v = {K = tostring(K), tau = tostring(tau), L = tostring(Ld)}}
end

local TRAIN_NAMES = {"Especificaciones (" .. U.zeta .. ", " .. U.wn .. ", " .. U.w .. "d)", "Diseno P/PD por polos", "Diseno PI/PID por polos",
  "Routh: K maximo", "Error estacionario", "Margenes MF/MG", "PI Ho-Hang-Cao", "PI Ziegler-Nichols + MF"}

local function m_train(v)
  local L = {}
  local tipo = num(v, "tipo") or 1
  tipo = math.floor(tipo)
  if not TRAIN[tipo] then error("tipo debe ser 1 a " .. #TRAIN) end
  local nro = math.floor(num(v, "nro") or 1)
  local P = TRAIN[tipo](rng(nro * 31 + tipo))
  sep(L, "Ejercicio " .. tipo .. "." .. nro .. ": " .. TRAIN_NAMES[tipo])
  for _, t in ipairs(P.txt) do add(L, t) end
  local given = {}
  local anyg, wantsol = false, false
  for i = 1, 3 do
    local s_ = v["r" .. i]
    if s_ and s_:match("%?") then wantsol = true
    elseif not blank(s_) then
      local ok, r = pcall(parse, s_, nil)
      if ok and risconst(r) then given[i] = rval(r); anyg = true end
    end
  end
  sep(L, "Se pide")
  for i = 1, 3 do
    if P.ans[i] then add(L, "r" .. i .. " = " .. P.ask[i]) end
  end
  if not anyg and not wantsol then
    add(L, "Resuelva a mano, escriba r1..r3 y presione enter.")
    add(L, "(escriba ? en r1 para ver la solucion; cambie 'nro' para otro ejercicio)")
    summary(L, {"Ejercicio " .. tipo .. "." .. nro .. " listo: lea el enunciado abajo"})
    return L
  end
  local S = {}
  local okc = 0
  local tot = 0
  for i = 1, 3 do
    local a = P.ans[i]
    if a then
      tot = tot + 1
      local g = given[i]
      local line
      if g == nil then line = "r" .. i .. " (" .. P.ask[i] .. "): sin respuesta; correcto = " .. fmt(a, 4)
      else
        local tol = math.max(0.02 * abs(a), 0.01)
        if P.ask[i]:find(U.deg) then tol = math.max(tol, 0.5) end
        if abs(g - a) <= tol then line = "r" .. i .. " (" .. P.ask[i] .. "): " .. fmt(g, 4) .. "  CORRECTO"; okc = okc + 1
        else line = "r" .. i .. " (" .. P.ask[i] .. "): " .. fmt(g, 4) .. "  INCORRECTO, correcto = " .. fmt(a, 4) end
      end
      S[#S + 1] = line
    end
  end
  table.insert(S, 1, "Puntaje: " .. okc .. "/" .. tot)
  sep(L, "Solucion paso a paso (modulo " .. P.mod .. ")")
  local ok, sol = pcall(MODULES[P.mod].run, P.v)
  if ok then
    for i = 1, #sol do add(L, sol[i]) end
    if sol.views then for _, vw in ipairs(sol.views) do view(L, vw) end end
  else
    add(L, "(no se pudo generar la solucion: " .. tostring(sol) .. ")")
  end
  summary(L, S)
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
  "'Curva de reaccion / identificar K, tau, L' " .. U.arrow .. " 9",
  "Practicar con ejercicios al azar " .. U.arrow .. " 10 (tecla 0 en el menu)",
  "== En los resultados ==",
  "Arriba: RESPUESTAS (lo que se escribe en la hoja). Abajo: desarrollo.",
  "tab: graficos (plano s, y(t), Bode, lugar de raices).",
  "Teclas 3/4/6/8: enviar el resultado a otro modulo (se indica en RESPUESTAS).",
  "Tecla c: guarda el resultado como caso A, abre los datos para modificarlos",
  "  y al calcular muestra la COMPARACION A " .. U.arrow .. " B con razonamiento.",
  "  (menu > Editar > Quitar comparacion para volver a lo normal)",
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
-- comparacion A -> B (tecla c)
--------------------------------------------------------------------------
local function rel(a, b)
  if a == b then return true end
  if a ~= a or b ~= b or abs(a) == INF or abs(b) == INF then return false end
  return abs(b - a) <= 1e-3 * math.max(abs(a), abs(b), 1e-9)
end

local function chg(a, b)
  if type(a) ~= "number" or type(b) ~= "number" then return "" end
  if a == INF or b == INF or a == -INF or b == -INF then return "" end
  if rel(a, b) then return "  (igual)" end
  if abs(a) > 1e-12 then
    local pct = (b - a) / abs(a) * 100
    return string.format("  (%s%.3g%%)", pct > 0 and "+" or "", pct)
  end
  return ""
end

local function vtxt(m)
  if type(m.val) ~= "number" then return tostring(m.val) end
  if m.key == "estable" then return (m.val == 1) and "estable" or "INESTABLE" end
  return fmt(m.val, 4) .. ((m.unit ~= "") and (" " .. m.unit) or "")
end

-- razonamiento por resultado
local function reason(key, a, b, MA, MB)
  local up = b > a
  local function g(k, M) return M[k] and M[k].val end
  local zA, zB = g("zeta", MA), g("zeta", MB)
  local sA, sB = g("sigma", MA), g("sigma", MB)
  local wdA, wdB = g("wd", MA), g("wd", MB)
  if key == "ts2" or key == "ts5" or key == "tsreal" then
    if sA and sB and rel(sA, sB) then return "ts no cambia: depende solo de " .. U.sigma .. " = " .. U.zeta .. U.wn .. " (parte real de los polos dominantes), que se mantuvo en " .. fmt(sA, 4) .. "." end
    if sA and sB then return "ts " .. (up and "aumenta" or "disminuye") .. " porque " .. U.sigma .. " = " .. U.zeta .. U.wn .. " pasa de " .. fmt(sA, 4) .. " a " .. fmt(sB, 4) .. ": los polos dominantes se " .. (up and "acercan al" or "alejan del") .. " eje j" .. U.w .. " (respuesta mas " .. (up and "lenta" or "rapida") .. ")." end
    return "El tiempo de asentamiento " .. (up and "aumenta: respuesta mas lenta." or "disminuye: respuesta mas rapida.")
  elseif key == "Mp" or key == "Mpreal" then
    if zA and zB and not rel(zA, zB) then return "Mp " .. (up and "aumenta" or "disminuye") .. " porque " .. U.zeta .. " pasa de " .. fmt(zA, 4) .. " a " .. fmt(zB, 4) .. " (" .. ((zB < zA) and "menos amortiguamiento: mas oscilacion" or "mas amortiguamiento: menos oscilacion") .. ")." end
    return "El sobrepaso " .. (up and "aumenta: respuesta mas oscilatoria." or "disminuye: respuesta mas amortiguada.")
  elseif key == "tp" or key == "tr" then
    if wdA and wdB then return key .. " " .. (up and "aumenta" or "disminuye") .. " porque depende de 1/" .. U.w .. "d y " .. U.w .. "d pasa de " .. fmt(wdA, 4) .. " a " .. fmt(wdB, 4) .. " rad/s." end
  elseif key == "wn" then
    return U.wn .. " (distancia de los polos al origen) " .. (up and "aumenta" or "disminuye") .. "."
  elseif key == "zeta" then
    return U.zeta .. " = cos(angulo del polo): " .. (up and "aumenta (polos mas cerca del eje real, menos oscilacion)." or "disminuye (polos mas cerca del eje j" .. U.w .. ", mas oscilacion).")
  elseif key == "wd" then
    return U.w .. "d (parte imaginaria de los polos) " .. (up and "aumenta: oscila mas rapido." or "disminuye: oscila mas lento.")
  elseif key == "ep" or key == "ev" or key == "ea" or key == "ess" then
    if a == 0 and b == 0 then return "El error sigue siendo nulo (el tipo del sistema lo anula)." end
    if b == INF then return "El error pasa a ser infinito (el tipo no alcanza para esa senal)." end
    local cst = ({ep = "Kp", ev = "Kv", ea = "Ka"})[key]
    local cA, cB = cst and g(cst, MA), cst and g(cst, MB)
    if cA and cB then return "El error " .. (up and "aumenta" or "disminuye") .. " porque " .. cst .. " pasa de " .. fmt(cA, 4) .. " a " .. fmt(cB, 4) .. " (" .. key .. " es inversamente proporcional a " .. cst .. ")." end
    return "El error en estado estable " .. (up and "aumenta: menor precision." or "disminuye: mayor precision.")
  elseif key == "MF" then
    return "MF " .. (up and "aumenta: mayor estabilidad relativa (menos sobrepaso)." or "disminuye: menor estabilidad relativa (mas oscilacion y sobrepaso)" .. ((b <= 0) and "; con MF " .. U.le .. " 0 el lazo es INESTABLE." or "."))
  elseif key == "MG" or key == "MGdB" then
    return "MG " .. (up and "aumenta: se puede subir mas la ganancia antes de inestabilizar." or "disminuye: el lazo queda mas cerca de la inestabilidad" .. ((key == "MG" and b <= 1) and "; con MG " .. U.le .. " 1 es INESTABLE." or "."))
  elseif key == "wcg" then
    return U.w .. "cg " .. (up and "aumenta: mayor ancho de banda, respuesta mas rapida." or "disminuye: respuesta mas lenta.")
  elseif key == "estable" then
    return (b == 1) and "El sistema pasa a ser ESTABLE." or "El sistema pasa a ser INESTABLE."
  elseif key == "ratio" then
    return "La separacion entre polos " .. (up and "aumenta" or "disminuye") .. ": la aproximacion de 2o orden es " .. ((b >= 5) and "valida (" .. U.ge .. " 5)." or "menos valida (< 5).")
  elseif key == "kmax" then
    return "El limite de estabilidad cambia: ahora se admite hasta k = " .. fmt(b, 5) .. "."
  elseif key == "znKp" then
    return "ZN: Kp = 0.9" .. U.tau .. "/(KL) es inversamente proporcional a K y a L: pasa de " .. fmt(a, 4) .. " a " .. fmt(b, 4) .. "."
  elseif key == "znTi" then
    return "ZN: Ti = 3.33L es proporcional al retardo: pasa de " .. fmt(a, 4) .. " a " .. fmt(b, 4) .. " s."
  elseif key == "hKp" then
    return "HHC: Kp = " .. U.w .. "p" .. U.tau .. "/(MG K) con " .. U.w .. "p inversamente proporcional a td: Kp pasa de " .. fmt(a, 4) .. " a " .. fmt(b, 4) .. "."
  elseif key == "hTi" then
    return "HHC: Ti cambia de " .. fmt(a, 4) .. " a " .. fmt(b, 4) .. " s."
  elseif key == "Ku" then
    return "La ganancia ultima Ku " .. (up and "aumenta: el proceso tolera mas ganancia." or "disminuye: con mas retardo/ganancia el lazo se inestabiliza antes.")
  elseif key == "tau" or key == "L" or key == "K" then
    return key .. " del modelo cambia de " .. fmt(a, 4) .. " a " .. fmt(b, 4) .. "."
  elseif key == "Kc" then
    return "La ganancia del controlador " .. (up and "aumenta: polos deseados mas exigentes (mas lejos) requieren mas ganancia." or "disminuye.")
  elseif key == "phi" then
    return "La deficiencia angular " .. (up and "aumenta: el controlador debe aportar mas fase." or "disminuye: se necesita menos fase.")
  end
  return nil
end

function compare_results(mod, VA, VB, MA, MB)
  local out = {}
  local m = MODULES[mod]
  out[#out + 1] = "== COMPARACION  A " .. U.arrow .. " B =="
  -- datos cambiados
  local changed = {}
  for _, f in ipairs(m.fields) do
    local a, b = VA[f[1]] or "", VB[f[1]] or ""
    if a ~= b then
      local txt = f[2] .. ": " .. ((a == "") and "(vacio)" or a) .. "  " .. U.arrow .. "  " .. ((b == "") and "(vacio)" or b)
      local na, nb = tonumber(a), tonumber(b)
      if na and nb and na ~= 0 then txt = txt .. "  (x" .. fmt(nb / na, 4) .. ")" end
      changed[#changed + 1] = txt
    end
  end
  out[#out + 1] = "Datos modificados:"
  if #changed == 0 then out[#out + 1] = "  (ninguno)" end
  for _, c in ipairs(changed) do out[#out + 1] = "  " .. c end
  -- resultados
  local MAk, MBk = {}, {}
  for _, x in ipairs(MA or {}) do MAk[x.key] = x end
  for _, x in ipairs(MB or {}) do MBk[x.key] = x end
  out[#out + 1] = "Resultados (A " .. U.arrow .. " B):"
  local reasons = {}
  for _, x in ipairs(MB or {}) do
    local a = MAk[x.key]
    if a then
      out[#out + 1] = "  " .. x.label .. ": " .. vtxt(a) .. "  " .. U.arrow .. "  " .. vtxt(x) .. chg(a.val, x.val)
      if type(a.val) == "number" and type(x.val) == "number" and not rel(a.val, x.val) then
        local r = reason(x.key, a.val, x.val, MAk, MBk)
        if r then reasons[#reasons + 1] = r end
      elseif type(a.val) == "number" and type(x.val) == "number" and (x.key == "ts2") then
        local r = reason(x.key, a.val, x.val, MAk, MBk)
        if r then reasons[#reasons + 1] = r end
      end
    else
      out[#out + 1] = "  " .. x.label .. ": (nuevo) " .. vtxt(x)
    end
  end
  if #reasons > 0 then
    out[#out + 1] = "== Razonamiento =="
    for _, r in ipairs(reasons) do out[#out + 1] = r end
  end
  -- conclusion
  local con = {}
  local function d(k) local a, b = MAk[k], MBk[k]; if a and b and type(a.val) == "number" and type(b.val) == "number" and not rel(a.val, b.val) then return b.val - a.val end end
  local dts = d("ts2") or d("tsreal")
  if dts then con[#con + 1] = (dts < 0) and "mas rapido" or "mas lento"
  elseif MAk["ts2"] and MBk["ts2"] then con[#con + 1] = "con el mismo tiempo de asentamiento" end
  local dmp = d("Mp") or d("Mpreal")
  if dmp then con[#con + 1] = (dmp > 0) and "mas oscilatorio (mayor Mp)" or "menos oscilatorio (menor Mp)" end
  local dmf = d("MF")
  if dmf then con[#con + 1] = (dmf < 0) and "menos estable relativamente (menor MF)" or "mas estable relativamente (mayor MF)" end
  local der = d("ev") or d("ep") or d("ess")
  if der then con[#con + 1] = (der < 0) and "mas preciso (menor error)" or "menos preciso (mayor error)" end
  local dst = d("estable")
  if dst then con[#con + 1] = (dst > 0) and "pasa a ser ESTABLE" or "pasa a ser INESTABLE" end
  local dwc = d("wcg")
  if dwc and not dts then con[#con + 1] = (dwc > 0) and "mas rapido (mayor " .. U.w .. "cg)" or "mas lento (menor " .. U.w .. "cg)" end
  if #con > 0 then
    out[#out + 1] = "== Conclusion =="
    out[#out + 1] = "Con el cambio" .. ((#changed > 0) and (" (" .. changed[1]:gsub("%s+", " ") .. ")") or "") .. ", el sistema resulta " .. table.concat(con, ", ") .. "."
  end
  out[#out + 1] = "== Resultado B completo =="
  return out
end

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
  {title = "Identificacion (curva)", short = "K, " .. U.tau .. ", L desde la curva de reaccion", run = m_ident, fields = {
    {"du", U.arrow .. "u escalon", "Tamano del escalon aplicado a la entrada", "1"},
    {"y0", "y inicial", "Valor de la salida antes del escalon", "0"},
    {"yf", "y final", "Valor final de la salida (estado estable)"},
    {"t28", "t 28.3%", "Tiempo en que y llega al 28.3% del cambio"},
    {"t63", "t 63.2%", "Tiempo en que y llega al 63.2% del cambio"},
    {"Lt", "L tangente", "(opcional) L leido con la tangente"},
    {"taut", U.tau .. " tangente", "(opcional) " .. U.tau .. " leido con la tangente"},
  }},
  {title = "Entrenamiento", short = "ejercicios al azar con correccion", run = m_train, fields = {
    {"tipo", "tipo", "1 espec. 2 PD 3 PI/PID 4 Routh 5 error 6 MF/MG 7 HHC 8 ZN", "1"},
    {"nro", "n" .. "\194\186" .. " ejercicio", "Cambie el numero para otro ejercicio", "1"},
    {"r1", "r1", "Su respuesta 1 (? = ver solucion)"},
    {"r2", "r2", "Su respuesta 2"},
    {"r3", "r3", "Su respuesta 3"},
  }},
  {title = "Ayuda", short = "teclas y sintaxis", help = HELP, fields = {}},
}

-- ===== ui.lua =====
--------------------------------------------------------------------------
-- ui.lua : interfaz TI-Nspire (Lua API 2.0)
--------------------------------------------------------------------------
local screen = "home"      -- home | form | result
local hsel = 1             -- modulo seleccionado en el menu
local cur = nil            -- modulo activo
local fsel = 1             -- campo seleccionado
local cpos = 0             -- posicion del cursor (bytes) en el campo
local fscroll = 0          -- desplazamiento del formulario
local rlines = {}          -- lineas de resultado (sin ajustar)
local rwrap, rwrapW = nil, -1
local rscroll = 0
local rerror = false
local rviews = nil         -- graficos del resultado (tab los recorre)
local vidx = 1
local rchain = nil         -- envio de datos a otros modulos (teclas)
local rmet = nil           -- resultados clave del ultimo calculo
local cmpbase = nil        -- caso A para comparar {mod, vals, met}
local lscroll = 0          -- desplazamiento de la leyenda del grafico
local VALS = {}            -- VALS[i][key] = texto

local C = {
  title = {0x1F, 0x3A, 0x68}, titletxt = {255, 255, 255},
  sel = {0xD8, 0xE6, 0xF7}, grey = {0x70, 0x70, 0x70}, black = {0, 0, 0},
  box = {0xA0, 0xA0, 0xA0}, focus = {0x1F, 0x5F, 0xC8}, hdr = {0xE3, 0xEB, 0xF7},
  hdrtxt = {0x1F, 0x3A, 0x68}, err = {0xC0, 0x10, 0x10}, white = {255, 255, 255},
  hint = {0xF4, 0xF4, 0xEC},
}
local function col(gc, c) gc:setColorRGB(c[1], c[2], c[3]) end

local function initvals()
  for i, m in ipairs(MODULES) do
    VALS[i] = VALS[i] or {}
    for _, f in ipairs(m.fields) do
      if VALS[i][f[1]] == nil then VALS[i][f[1]] = f[4] or "" end
    end
  end
end
initvals()

local function inval() platform.window:invalidate() end

local function W() return platform.window:width() end
local function H() return platform.window:height() end

-- UTF-8: saltos de caracter
local function prevchar(s, p)
  if p <= 0 then return 0 end
  p = p - 1
  while p > 0 do
    local b = s:byte(p + 1)
    if b < 0x80 or b >= 0xC0 then break end
    p = p - 1
  end
  return p
end
local function nextchar(s, p)
  if p >= #s then return #s end
  p = p + 1
  while p < #s do
    local b = s:byte(p + 1)
    if b < 0x80 or b >= 0xC0 then break end
    p = p + 1
  end
  return p
end

local function curfield() return cur and MODULES[cur].fields[fsel] end
local function getval() local f = curfield(); return f and VALS[cur][f[1]] or "" end
local function setval(s) local f = curfield(); if f then VALS[cur][f[1]] = s end end

local function openmodule(i)
  if cmpbase and cmpbase.mod ~= i then cmpbase = nil end
  cur = i
  local m = MODULES[i]
  if m.help then
    rlines = m.help; rwrap = nil; rscroll = 0; rerror = false
    screen = "result"
  else
    screen = "form"
    fsel = 1; fscroll = 0
    cpos = #getval()
  end
  inval()
end

local function compute()
  local m = MODULES[cur]
  local ok, res = pcall(m.run, VALS[cur])
  if ok then
    rlines = res; rerror = false; rviews = res.views; rchain = res.chain; vidx = 1; lscroll = 0
    rmet = res.met
    if cmpbase and cmpbase.mod == cur and res.met then
      local okc, cl = pcall(compare_results, cur, cmpbase.vals, VALS[cur], cmpbase.met, res.met)
      if okc then
        local merged = {}
        for _, l in ipairs(cl) do merged[#merged + 1] = l end
        for _, l in ipairs(res) do merged[#merged + 1] = l end
        merged.views = res.views; merged.chain = res.chain; merged.met = res.met
        rlines = merged
      end
    end
    if rviews and #rviews == 0 then rviews = nil end
  else
    local msg = tostring(res):gsub("^.-:%d+: ", "")
    rlines = {"ERROR", msg, "", "Revise los datos (esc para volver)."}
    rviews = nil; rchain = nil
    rerror = true
  end
  rwrap = nil; rscroll = 0
  screen = "result"
  inval()
end

--------------------------------------------------------------------------
-- dibujo
--------------------------------------------------------------------------
local function titlebar(gc, txt, right)
  col(gc, C.title); gc:fillRect(0, 0, W(), 17)
  col(gc, C.titletxt); gc:setFont("sansserif", "b", 10)
  gc:drawString(txt, 4, 1, "top")
  if right then
    gc:setFont("sansserif", "r", 7)
    local w = gc:getStringWidth(right)
    gc:drawString(right, W() - w - 4, 3, "top")
  end
end

local function wraptext(gc, s, maxw)
  local out = {}
  if gc:getStringWidth(s) <= maxw then out[1] = s; return out end
  local indent = s:match("^(%s*)") or ""
  local line = indent
  for word in s:gmatch("%S+%s*") do
    local test = line .. word
    if gc:getStringWidth(test) > maxw and line ~= indent then
      out[#out + 1] = line
      line = indent .. "  " .. word
    else
      line = test
    end
    -- palabra mas larga que la linea: cortar
    while gc:getStringWidth(line) > maxw do
      local cut = #line
      while cut > 1 and gc:getStringWidth(line:sub(1, cut)) > maxw do cut = prevchar(line, cut) end
      out[#out + 1] = line:sub(1, cut)
      line = indent .. "  " .. line:sub(cut + 1)
    end
  end
  if line ~= indent then out[#out + 1] = line end
  return out
end

local HOME = {y0 = 19, rh = 18, vis = 9, scroll = 0}

local function paint_home(gc)
  titlebar(gc, "Control Automatico", "v11 - APN")
  -- alturas reales de las fuentes (evita superposiciones)
  gc:setFont("sansserif", "b", 10)
  local th = gc:getStringHeight("Ay")
  gc:setFont("sansserif", "r", 9)
  local fh = gc:getStringHeight("Ay")
  local footH = fh + 4
  local y0 = 19
  local rh = th + 2
  local n = #MODULES
  local vis = math.max(1, math.floor((H() - y0 - footH) / rh))
  if vis > n then vis = n end
  local sc = HOME.scroll
  if hsel > sc + vis then sc = hsel - vis end
  if hsel <= sc then sc = hsel - 1 end
  HOME.y0, HOME.rh, HOME.vis, HOME.scroll = y0, rh, vis, sc
  for i = sc + 1, math.min(n, sc + vis) do
    local m = MODULES[i]
    local y = y0 + (i - sc - 1) * rh
    if i == hsel then col(gc, C.sel); gc:fillRect(0, y, W(), rh) end
    gc:setFont("sansserif", "b", 10)
    col(gc, C.focus)
    gc:drawString(tostring(i), 6, y + 1, "top")
    col(gc, C.black)
    gc:drawString(m.title, 22, y + 1, "top")
  end
  if n > vis then
    local tot = vis * rh
    local bh = math.max(8, tot * vis / n)
    local by = y0 + (tot - bh) * sc / (n - vis)
    col(gc, C.focus); gc:fillRect(W() - 3, by, 2, bh)
  end
  -- pie: descripcion del modulo seleccionado
  local fy = H() - footH
  col(gc, C.hint); gc:fillRect(0, fy, W(), footH)
  col(gc, C.box); gc:drawLine(0, fy, W(), fy)
  gc:setFont("sansserif", "r", 9)
  col(gc, C.grey)
  local txt = MODULES[hsel].short .. "   (enter)"
  while gc:getStringWidth(txt) > W() - 8 and #txt > 1 do txt = txt:sub(1, prevchar(txt, #txt)) end
  gc:drawString(txt, 4, fy + 2, "top")
end

local FH = 18
local function formlayout()
  local top = 19
  local hintH = 30
  local vis = math.floor((H() - top - hintH) / FH)
  return top, hintH, vis
end

local function paint_form(gc)
  local m = MODULES[cur]
  titlebar(gc, cur .. ". " .. m.title, (cmpbase and cmpbase.mod == cur) and "enter: comparar con A" or "enter: calcular")
  local top, hintH, vis = formlayout()
  local fields = m.fields
  if fsel > fscroll + vis then fscroll = fsel - vis end
  if fsel <= fscroll then fscroll = fsel - 1 end
  local lw = 74
  for idx = fscroll + 1, math.min(#fields, fscroll + vis) do
    local f = fields[idx]
    local y = top + (idx - fscroll - 1) * FH
    gc:setFont("sansserif", idx == fsel and "b" or "r", 9)
    col(gc, idx == fsel and C.focus or C.black)
    local lab = f[2]
    while gc:getStringWidth(lab) > lw - 4 and #lab > 1 do lab = lab:sub(1, prevchar(lab, #lab)) end
    gc:drawString(lab, 3, y + 2, "top")
    local bx, bw = lw, W() - lw - 4
    col(gc, C.white); gc:fillRect(bx, y + 1, bw, FH - 3)
    if idx == fsel then
      col(gc, C.focus); gc:drawRect(bx, y + 1, bw, FH - 3); gc:drawRect(bx + 1, y + 2, bw - 2, FH - 5)
    else
      col(gc, C.box); gc:drawRect(bx, y + 1, bw, FH - 3)
    end
    local txt = VALS[cur][f[1]] or ""
    gc:setFont("sansserif", "r", 9)
    col(gc, C.black)
    local avail = bw - 8
    local start = 0
    if idx == fsel then
      while gc:getStringWidth(txt:sub(start + 1, cpos)) > avail do start = nextchar(txt, start) end
    end
    local shown = txt:sub(start + 1)
    while gc:getStringWidth(shown) > avail and #shown > 0 do shown = shown:sub(1, prevchar(shown, #shown)) end
    gc:drawString(shown, bx + 3, y + 2, "top")
    if idx == fsel then
      local cx = bx + 3 + gc:getStringWidth(txt:sub(start + 1, cpos))
      col(gc, C.focus); gc:fillRect(cx, y + 3, 1, FH - 7)
    end
  end
  -- barra de scroll
  if #fields > vis then
    local th = vis * FH
    local bh = math.max(8, th * vis / #fields)
    local by = top + (th - bh) * fscroll / (#fields - vis)
    col(gc, C.box); gc:fillRect(W() - 3, by, 2, bh)
  end
  -- ayuda del campo
  local y = H() - hintH
  col(gc, C.hint); gc:fillRect(0, y, W(), hintH)
  col(gc, C.box); gc:drawLine(0, y, W(), y)
  local f = fields[fsel]
  gc:setFont("sansserif", "r", 9)
  col(gc, C.black)
  if f then
    local ls = wraptext(gc, f[3] or "", W() - 6)
    for i = 1, math.min(2, #ls) do gc:drawString(ls[i], 3, y + 1 + (i - 1) * 14, "top") end
  end
end

local function paint_result(gc)
  local m = MODULES[cur]
  titlebar(gc, cur .. ". " .. m.title, rviews and ("tab: grafico (" .. #rviews .. ")") or "esc: volver")
  gc:setFont("sansserif", "r", 9)
  local lh = gc:getStringHeight("Ay")
  if lh < 11 then lh = 11 end
  local maxw = W() - 10
  if not rwrap or rwrapW ~= maxw then
    rwrap = {}
    for _, s in ipairs(rlines) do
      if s:sub(1, 3) == "== " then
        rwrap[#rwrap + 1] = {s:gsub("^== ", ""):gsub(" ==$", ""), "h"}
      else
        for _, w in ipairs(wraptext(gc, s, maxw)) do rwrap[#rwrap + 1] = {w, "n"} end
      end
    end
    rwrapW = maxw
  end
  local top = 19
  local vis = math.floor((H() - top - 2) / lh)
  local maxs = math.max(0, #rwrap - vis)
  if rscroll > maxs then rscroll = maxs end
  if rscroll < 0 then rscroll = 0 end
  for i = 1, vis do
    local e = rwrap[rscroll + i]
    if not e then break end
    local y = top + (i - 1) * lh
    if e[2] == "h" then
      col(gc, C.hdr); gc:fillRect(0, y, W() - 5, lh)
      col(gc, C.hdrtxt); gc:setFont("sansserif", "b", 9)
      gc:drawString(e[1], 3, y, "top")
      gc:setFont("sansserif", "r", 9)
    else
      col(gc, (rerror and rscroll + i <= 2) and C.err or C.black)
      gc:drawString(e[1], 3, y, "top")
    end
  end
  if #rwrap > vis then
    local th = H() - top - 2
    local bh = math.max(8, th * vis / #rwrap)
    local by = top + (th - bh) * rscroll / maxs
    col(gc, C.box); gc:fillRect(W() - 4, top, 3, th)
    col(gc, C.focus); gc:fillRect(W() - 4, by, 3, bh)
  end
end

--------------------------------------------------------------------------
-- grafico del plano s (localizacion de polos)
--------------------------------------------------------------------------
local GCOL = {
  pp = {0xD0, 0x10, 0x10}, pz = {0xD0, 0x10, 0x10},      -- planta: rojo
  pd = {0x10, 0x40, 0xE0},                               -- deseados: azul
  ci = {0x00, 0x96, 0x20}, cz = {0x00, 0x96, 0x20},      -- controlador: verde
  lc = {0x00, 0x96, 0x20}, txt = {0, 0, 0},
}

local function marker(gc, kind, x, y)
  col(gc, GCOL[kind] or C.black)
  if kind == "pp" or kind == "ci" or kind == "lc" then
    gc:drawLine(x - 4, y - 4, x + 4, y + 4); gc:drawLine(x - 4, y + 4, x + 4, y - 4)
    gc:drawLine(x - 3, y - 4, x + 5, y + 4); gc:drawLine(x - 3, y + 4, x + 5, y - 4)
  elseif kind == "pz" or kind == "cz" then
    gc:drawArc(x - 4, y - 4, 8, 8, 0, 360); gc:drawArc(x - 3, y - 3, 6, 6, 0, 360)
  elseif kind == "pd" then
    gc:fillArc(x - 4, y - 4, 8, 8, 0, 360)
  end
end

local function vtitle(gc, t)
  titlebar(gc, t, "tab " .. vidx .. "/" .. #rviews)
end

local function paint_splane(gc, P)
  vtitle(gc, "Plano s: angulos hacia pd")
  gc:setFont("sansserif", "r", 7)
  local lh = gc:getStringHeight("Ay")
  if lh < 9 then lh = 9 end
  local nleg = #P.legend + 1
  local maxleg = math.max(3, math.floor((H() - 18) * 0.45 / lh))
  local shown = math.min(nleg, maxleg)
  if lscroll > nleg - shown then lscroll = nleg - shown end
  if lscroll < 0 then lscroll = 0 end
  local legH = shown * lh + 3
  local top, bot = 19, H() - legH - 1
  local left, right = 4, W() - 4
  -- rango
  local sig, wd = P.sig, P.wd
  local R = math.max(sig, wd, 1e-6)
  local lim = 6 * R
  local xmin, xmax, ymin, ymax = -sig, 0, 0, wd
  for _, it in ipairs(P.items) do
    if math.abs(it.x) <= lim and math.abs(it.y) <= lim then
      if it.x < xmin then xmin = it.x end
      if it.x > xmax then xmax = it.x end
      if it.y < ymin then ymin = it.y end
      if it.y > ymax then ymax = it.y end
    end
  end
  local mx = 0.18 * R
  xmin, xmax = xmin - mx, xmax + 2.2 * mx
  ymax = ymax + 0.15 * R
  ymin = math.min(ymin - 0.1 * R, -0.25 * R)
  local sc = math.min((right - left) / (xmax - xmin), (bot - top) / (ymax - ymin))
  local cxw = (right - left) - sc * (xmax - xmin)
  local cyh = (bot - top) - sc * (ymax - ymin)
  local function X(x) return left + cxw / 2 + (x - xmin) * sc end
  local function Y(y) return top + cyh / 2 + (ymax - y) * sc end
  -- ejes
  col(gc, C.box)
  gc:drawLine(left, Y(0), right, Y(0))
  if X(0) >= left and X(0) <= right then gc:drawLine(X(0), top, X(0), bot) end
  col(gc, C.grey)
  gc:drawString(U.sigma, right - 8, Y(0) - lh, "top")
  gc:drawString("j" .. U.w, X(0) + 3, top, "top")
  local pdx, pdy = X(-sig), Y(wd)
  -- lineas y angulos
  local placed = {}
  for _, it in ipairs(P.items) do
    if it.ang then
      local x, y = X(it.x), Y(it.y)
      local clipped = false
      if x < left or x > right or y < top or y > bot then
        -- recortar sobre el borde, en la misma direccion
        local t = 1
        local function fit(a0, a1, lo, hi)
          if a1 < lo then t = math.min(t, (lo - a0) / (a1 - a0)) end
          if a1 > hi then t = math.min(t, (hi - a0) / (a1 - a0)) end
        end
        fit(pdx, x, left, right); fit(pdy, y, top, bot)
        x, y = pdx + (x - pdx) * t, pdy + (y - pdy) * t
        clipped = true
      end
      col(gc, GCOL[it.kind] or C.black)
      gc:drawLine(x, y, pdx, pdy)
      -- arcos escalonados si hay raices muy juntas
      local r = 13
      for _, q in ipairs(placed) do
        if math.abs(q[1] - x) < 16 and math.abs(q[2] - y) < 16 then r = r + 9 end
      end
      placed[#placed + 1] = {x, y}
      local a0 = 0
      local aa = it.ang
      gc:drawArc(x - r, y - r, 2 * r, 2 * r, a0, aa)
      local am = math.rad(aa / 2)
      gc:setFont("sansserif", "b", 9)
      local lx = x + (r + 3) * math.cos(am) - 3
      local ly2 = y - (r + 3) * math.sin(am) - 7
      lx = math.max(left, math.min(right - 10, lx))
      ly2 = math.max(top, math.min(bot - 14, ly2))
      gc:drawString(it.label, lx, ly2, "top")
      gc:setFont("sansserif", "r", 7)
      if clipped then
        gc:drawString(fmt(it.x, 3), x + 2, y + 2, "top")
      end
      if not clipped then marker(gc, it.kind, x, y) end
    end
  end
  -- polos de lazo cerrado adicionales
  for _, it in ipairs(P.items) do
    if not it.ang then
      local x, y = X(it.x), Y(it.y)
      if x >= left and x <= right and y >= top and y <= bot then marker(gc, it.kind, x, y) end
    end
  end
  -- polos deseados
  marker(gc, "pd", pdx, pdy)
  local cy = Y(-wd)
  if cy <= bot then marker(gc, "pd", pdx, cy) end
  col(gc, GCOL.pd)
  gc:setFont("sansserif", "r", 7)
  local pl = "pd=" .. fmtc(-sig, wd, 3)
  local pw = gc:getStringWidth(pl)
  local px = pdx + 6
  if px + pw > right then px = pdx - 6 - pw end
  if px < left then px = left end
  local py = pdy - lh - 2
  if py < top then py = pdy + 5 end
  gc:drawString(pl, px, py, "top")
  -- leyenda
  local ly = H() - legH
  col(gc, C.hint); gc:fillRect(0, ly, W(), legH)
  col(gc, C.box); gc:drawLine(0, ly, W(), ly)
  local lines = {{"rojo: planta   azul: pd deseado   verde: controlador / polos LC", "txt"}}
  for _, l in ipairs(P.legend) do lines[#lines + 1] = l end
  for i = 1, shown do
    local l = lines[lscroll + i]
    if l then
      col(gc, GCOL[l[2]] or C.black)
      gc:drawString(l[1], 3, ly + 2 + (i - 1) * lh, "top")
    end
  end
  if nleg > shown then
    col(gc, C.grey); gc:drawString(U.arrow, W() - 10, ly + 2, "top")
  end
end

-- utilidades de ejes ----------------------------------------------------
local function nicestep(range, n)
  local raw = range / n
  if raw <= 0 then return 1 end
  local e = 10 ^ math.floor(math.log10(raw))
  local f = raw / e
  local st = (f < 1.5) and 1 or ((f < 3.5) and 2 or ((f < 7.5) and 5 or 10))
  return st * e
end

local function frame(gc, x0, y0, x1, y1, xmin, xmax, ymin, ymax, logx, nox)
  local function X(x) return x0 + (x - xmin) / (xmax - xmin) * (x1 - x0) end
  local function Y(y) return y1 - (y - ymin) / (ymax - ymin) * (y1 - y0) end
  col(gc, C.box); gc:drawRect(x0, y0, x1 - x0, y1 - y0)
  gc:setFont("sansserif", "r", 7)
  -- marcas eje y
  local sy = nicestep(ymax - ymin, 4)
  local v = math.ceil(ymin / sy) * sy
  while v <= ymax + 1e-9 do
    local yy = Y(v)
    col(gc, {0xE4, 0xE4, 0xE4}); gc:drawLine(x0 + 1, yy, x1 - 1, yy)
    col(gc, C.grey); gc:drawString(fmt(v, 3), 1, yy - 6, "top")
    v = v + sy
  end
  -- marcas eje x
  if logx then
    for d = math.ceil(xmin), math.floor(xmax) do
      local xx = X(d)
      col(gc, {0xE4, 0xE4, 0xE4}); gc:drawLine(xx, y0 + 1, xx, y1 - 1)
      if not nox then col(gc, C.grey); gc:drawString("1e" .. d, xx - 8, y1 + 1, "top") end
    end
  else
    local sx = nicestep(xmax - xmin, 5)
    local u = math.ceil(xmin / sx) * sx
    while u <= xmax + 1e-9 do
      local xx = X(u)
      col(gc, {0xE4, 0xE4, 0xE4}); gc:drawLine(xx, y0 + 1, xx, y1 - 1)
      col(gc, C.grey); gc:drawString(fmt(u, 3), xx - 6, y1 + 1, "top")
      u = u + sx
    end
  end
  return X, Y
end

local function dashed(gc, on_)
  if gc.setPen then pcall(function() gc:setPen("thin", on_ and "dashed" or "smooth") end) end
end

local function legendbox(gc, lines, nmax)
  gc:setFont("sansserif", "r", 7)
  local lh = math.max(9, gc:getStringHeight("Ay"))
  local n = math.min(#lines, nmax or 4)
  local hgt = n * lh + 3
  local ly = H() - hgt
  col(gc, C.hint); gc:fillRect(0, ly, W(), hgt)
  col(gc, C.box); gc:drawLine(0, ly, W(), ly)
  local sc = math.max(0, math.min(lscroll, #lines - n))
  for i = 1, n do
    local l = lines[sc + i]
    if l then
      col(gc, l[2] or C.black)
      gc:drawString(l[1], 3, ly + 2 + (i - 1) * lh, "top")
    end
  end
  return ly
end

local CURVE = {{0x10, 0x40, 0xE0}, {0xD0, 0x10, 0x10}, {0x00, 0x96, 0x20}, {0x90, 0x30, 0xB0}}

-- respuesta temporal ----------------------------------------------------
local function paint_time(gc, V)
  vtitle(gc, V.title or "Respuesta y(t)")
  local leg = {}
  for i, c in ipairs(V.curves) do leg[#leg + 1] = {c.name, CURVE[i]} end
  for _, l in ipairs(V.legend or {}) do leg[#leg + 1] = {l, C.black} end
  local ly = legendbox(gc, leg, 4)
  local tmax, ymin, ymax = 0, 0, 0
  for _, c in ipairs(V.curves) do
    for i = 1, #c.t do
      if c.t[i] > tmax then tmax = c.t[i] end
      if c.y[i] > ymax then ymax = c.y[i] end
      if c.y[i] < ymin then ymin = c.y[i] end
    end
  end
  if V.ref then ymax = math.max(ymax, V.ref) end
  local pad = 0.08 * (ymax - ymin + 1e-9)
  ymax, ymin = ymax + pad, ymin - pad
  if tmax <= 0 then tmax = 1 end
  local X, Y = frame(gc, 26, 20, W() - 4, ly - 11, 0, tmax, ymin, ymax, false)
  -- referencia y banda de 2%
  if V.ref then
    col(gc, C.grey); dashed(gc, true)
    gc:drawLine(X(0), Y(V.ref), X(tmax), Y(V.ref))
    if V.band then
      col(gc, {0xB0, 0xB0, 0xB0})
      gc:drawLine(X(0), Y(V.ref * 1.02), X(tmax), Y(V.ref * 1.02))
      gc:drawLine(X(0), Y(V.ref * 0.98), X(tmax), Y(V.ref * 0.98))
    end
    dashed(gc, false)
  end
  for i, c in ipairs(V.curves) do
    col(gc, CURVE[i])
    for k = 2, #c.t do gc:drawLine(X(c.t[k - 1]), Y(c.y[k - 1]), X(c.t[k]), Y(c.y[k])) end
  end
  -- marcas de la primera curva
  local m = V.marks
  if m then
    gc:setFont("sansserif", "b", 7)
    if m.tp and m.ymax then
      col(gc, CURVE[1]); gc:fillArc(X(m.tp) - 3, Y(m.ymax) - 3, 6, 6, 0, 360)
      gc:drawString("Mp=" .. fmt(m.Mp, 3) .. "%", math.min(X(m.tp) + 4, W() - 60), Y(m.ymax) - 2, "top")
    end
    if m.ts then
      col(gc, CURVE[3]); dashed(gc, true)
      gc:drawLine(X(m.ts), 21, X(m.ts), ly - 12)
      dashed(gc, false)
      gc:drawString("ts", X(m.ts) + 2, 22, "top")
    end
  end
  col(gc, C.grey); gc:setFont("sansserif", "r", 7)
  gc:drawString("t [s]", W() - 28, ly - 22, "top")
end

-- Bode --------------------------------------------------------------------
local function paint_bode(gc, V)
  vtitle(gc, "Bode de G0(j" .. U.w .. ")")
  local leg = {}
  for _, l in ipairs(V.legend or {}) do leg[#leg + 1] = {l[1], l[2] and GCOL[l[2]] or C.black} end
  local ly = legendbox(gc, leg, 3)
  local top, bot = 19, ly - 10
  local mid = math.floor((top + bot) / 2)
  local lw0, lw1 = V.lw[1], V.lw[#V.lw]
  local mmin, mmax, pmin, pmax = 0, 0, -180, -90
  for i = 1, #V.lw do
    mmin = math.min(mmin, V.mag[i]); mmax = math.max(mmax, V.mag[i])
    pmin = math.min(pmin, V.ph[i]); pmax = math.max(pmax, V.ph[i])
  end
  mmin, mmax = math.max(mmin, -80) - 3, math.min(mmax, 80) + 3
  pmin, pmax = math.max(pmin, -540) - 10, math.min(pmax, 180) + 10
  local X, Ym = frame(gc, 26, top, W() - 4, mid - 2, lw0, lw1, mmin, mmax, true, true)
  local _, Yp = frame(gc, 26, mid + 2, W() - 4, bot, lw0, lw1, pmin, pmax, true)
  gc:setFont("sansserif", "r", 7); col(gc, C.grey)
  gc:drawString("dB", 28, top + 1, "top"); gc:drawString("fase " .. U.deg, 28, mid + 3, "top")
  -- 0 dB y -180
  col(gc, C.grey); dashed(gc, true)
  gc:drawLine(X(lw0), Ym(0), X(lw1), Ym(0))
  gc:drawLine(X(lw0), Yp(-180), X(lw1), Yp(-180))
  dashed(gc, false)
  col(gc, CURVE[1])
  for i = 2, #V.lw do
    local a, b = math.max(mmin, math.min(mmax, V.mag[i - 1])), math.max(mmin, math.min(mmax, V.mag[i]))
    gc:drawLine(X(V.lw[i - 1]), Ym(a), X(V.lw[i]), Ym(b))
    local c, d = math.max(pmin, math.min(pmax, V.ph[i - 1])), math.max(pmin, math.min(pmax, V.ph[i]))
    gc:drawLine(X(V.lw[i - 1]), Yp(c), X(V.lw[i]), Yp(d))
  end
  gc:setFont("sansserif", "b", 7)
  if V.wcg then
    local x = X(math.log10(V.wcg))
    col(gc, GCOL.cz); dashed(gc, true); gc:drawLine(x, top + 1, x, bot - 1); dashed(gc, false)
    gc:drawLine(x, Yp(V.phcg), x, Yp(-180))
    gc:fillArc(x - 2, Yp(V.phcg) - 2, 4, 4, 0, 360)
    gc:drawString("MF=" .. fmt(V.MF, 3) .. U.deg, math.min(x + 3, W() - 50), (Yp(V.phcg) + Yp(-180)) / 2 - 5, "top")
    gc:drawString(U.w .. "cg", math.min(x + 2, W() - 22), top + 1, "top")
  end
  if V.wcf then
    local x = X(math.log10(V.wcf))
    col(gc, GCOL.pp); dashed(gc, true); gc:drawLine(x, top + 1, x, bot - 1); dashed(gc, false)
    gc:drawLine(x, Ym(V.magcf), x, Ym(0))
    gc:fillArc(x - 2, Ym(V.magcf) - 2, 4, 4, 0, 360)
    gc:drawString("MG=" .. fmt(-V.magcf, 3) .. "dB", math.min(x + 3, W() - 55), (Ym(V.magcf) + Ym(0)) / 2 - 5, "top")
    gc:drawString(U.w .. "cf", math.min(x + 2, W() - 22), mid + 3, "top")
  end
end

-- lugar geometrico de las raices ------------------------------------------
local function paint_rlocus(gc, V)
  vtitle(gc, "Lugar de raices (k > 0)")
  local leg = {}
  for _, l in ipairs(V.legend or {}) do leg[#leg + 1] = {l[1], l[2] and GCOL[l[2]] or C.black} end
  local ly = legendbox(gc, leg, 3)
  local top, bot, left, right = 19, ly - 2, 4, W() - 4
  local R = V.R
  local xmin, xmax, ymin, ymax = V.xmin, V.xmax, -V.ymax, V.ymax
  local sc = math.min((right - left) / (xmax - xmin), (bot - top) / (ymax - ymin))
  local cxw = (right - left) - sc * (xmax - xmin)
  local cyh = (bot - top) - sc * (ymax - ymin)
  local function X(x) return left + cxw / 2 + (x - xmin) * sc end
  local function Y(y) return top + cyh / 2 + (ymax - y) * sc end
  col(gc, C.box)
  gc:drawLine(left, Y(0), right, Y(0))
  if X(0) >= left and X(0) <= right then gc:drawLine(X(0), top, X(0), bot) end
  col(gc, CURVE[1])
  for _, p in ipairs(V.pts) do
    local x, y = X(p[1]), Y(p[2])
    if x >= left and x <= right and y >= top and y <= bot then gc:fillRect(x, y, 2, 2) end
  end
  for _, p in ipairs(V.poles) do marker(gc, "pp", X(p[1]), Y(p[2])) end
  for _, z in ipairs(V.zeros) do marker(gc, "pz", X(z[1]), Y(z[2])) end
  gc:setFont("sansserif", "b", 7)
  for _, c in ipairs(V.cross) do
    col(gc, GCOL.cz)
    gc:fillArc(X(0) - 3, Y(c.w) - 3, 6, 6, 0, 360)
    gc:fillArc(X(0) - 3, Y(-c.w) - 3, 6, 6, 0, 360)
    gc:drawString("k=" .. fmt(c.k, 4), math.min(X(0) + 4, right - 50), Y(c.w) - 4, "top")
  end
  if V.kpoles then
    for _, p in ipairs(V.kpoles) do
      col(gc, GCOL.pd); gc:fillArc(X(p[1]) - 3, Y(p[2]) - 3, 6, 6, 0, 360)
    end
  end
  col(gc, C.grey); gc:setFont("sansserif", "r", 7)
  gc:drawString(U.sigma, right - 8, Y(0) - 11, "top")
  gc:drawString("j" .. U.w, X(0) + 3, top, "top")
end

-- mapa de polos con dominancia -------------------------------------------
local function paint_poles(gc, V)
  vtitle(gc, "Polos LC y dominancia")
  local leg = {}
  for _, l in ipairs(V.legend or {}) do leg[#leg + 1] = {l[1], l[2] and GCOL[l[2]] or C.black} end
  local ly = legendbox(gc, leg, 3)
  local top, bot, left, right = 19, ly - 2, 4, W() - 4
  local xmin, xmax, ymax = 0, 0, 0
  for _, p in ipairs(V.poles) do xmin = math.min(xmin, p[1]); xmax = math.max(xmax, p[1]); ymax = math.max(ymax, math.abs(p[2])) end
  for _, p in ipairs(V.zeros) do xmin = math.min(xmin, p[1]); xmax = math.max(xmax, p[1]); ymax = math.max(ymax, math.abs(p[2])) end
  local five = 5 * V.dom[1]
  if V.dom[1] < 0 then xmin = math.min(xmin, five) end
  local R = math.max(xmax - xmin, 2 * ymax, 1e-3)
  xmin, xmax = xmin - 0.08 * R, xmax + 0.12 * R
  ymax = math.max(ymax, 0.15 * R) * 1.2
  local sc = math.min((right - left) / (xmax - xmin), (bot - top) / (2 * ymax))
  local cxw = (right - left) - sc * (xmax - xmin)
  local cyh = (bot - top) - sc * 2 * ymax
  local function X(x) return left + cxw / 2 + (x - xmin) * sc end
  local function Y(y) return top + cyh / 2 + (ymax - y) * sc end
  col(gc, C.box)
  gc:drawLine(left, Y(0), right, Y(0))
  gc:drawLine(X(0), top, X(0), bot)
  if V.dom[1] < 0 then
    col(gc, GCOL.cz); dashed(gc, true)
    gc:drawLine(X(five), top, X(five), bot)
    dashed(gc, false)
    gc:setFont("sansserif", "r", 7)
    gc:drawString("5" .. "\194\183" .. "Re dom", X(five) + 2, top + 1, "top")
    col(gc, GCOL.pd); dashed(gc, true)
    gc:drawLine(X(V.dom[1]), top, X(V.dom[1]), bot)
    dashed(gc, false)
    gc:drawString("Re dom", X(V.dom[1]) + 2, bot - 12, "top")
  end
  for _, p in ipairs(V.zeros) do marker(gc, "pz", X(p[1]), Y(p[2])) end
  for i, p in ipairs(V.poles) do
    local isdom = (i <= V.nd)
    if isdom then
      col(gc, GCOL.pd)
      gc:fillArc(X(p[1]) - 4, Y(p[2]) - 4, 8, 8, 0, 360)
    else
      marker(gc, "pp", X(p[1]), Y(p[2]))
    end
    gc:setFont("sansserif", "r", 7)
    col(gc, isdom and GCOL.pd or GCOL.pp)
    if p[2] >= 0 then gc:drawString(fmtc(p[1], p[2], 3), math.min(X(p[1]) + 5, right - 60), Y(p[2]) - 12, "top") end
  end
  col(gc, C.grey); gc:setFont("sansserif", "r", 7)
  gc:drawString(U.sigma, right - 8, Y(0) - 11, "top")
  gc:drawString("j" .. U.w, X(0) + 3, top, "top")
end

local function paint_graph(gc)
  local V = rviews[vidx]
  if V.kind == "splane" then paint_splane(gc, V)
  elseif V.kind == "time" then paint_time(gc, V)
  elseif V.kind == "bode" then paint_bode(gc, V)
  elseif V.kind == "rlocus" then paint_rlocus(gc, V)
  elseif V.kind == "poles" then paint_poles(gc, V) end
end

function on.paint(gc)
  col(gc, C.white); gc:fillRect(0, 0, W(), H())
  if screen == "home" then paint_home(gc)
  elseif screen == "form" then paint_form(gc)
  elseif screen == "graph" and rviews then paint_graph(gc)
  else paint_result(gc) end
end

function on.resize() rwrap = nil; inval() end

--------------------------------------------------------------------------
-- teclado
--------------------------------------------------------------------------
local function pagesize()
  return math.max(1, math.floor((H() - 21) / 14) - 1)
end

function on.arrowUp()
  if screen == "home" then hsel = (hsel - 2) % #MODULES + 1
  elseif screen == "form" then
    local n = #MODULES[cur].fields
    fsel = (fsel - 2) % n + 1; cpos = #getval()
  elseif screen == "graph" then lscroll = lscroll - 1
  else rscroll = rscroll - 1 end
  inval()
end

function on.arrowDown()
  if screen == "home" then hsel = hsel % #MODULES + 1
  elseif screen == "form" then
    local n = #MODULES[cur].fields
    fsel = fsel % n + 1; cpos = #getval()
  elseif screen == "graph" then lscroll = lscroll + 1
  else rscroll = rscroll + 1 end
  inval()
end

function on.arrowLeft()
  if screen == "form" then cpos = prevchar(getval(), cpos)
  elseif screen == "result" then rscroll = rscroll - pagesize() end
  inval()
end

function on.arrowRight()
  if screen == "form" then cpos = nextchar(getval(), cpos)
  elseif screen == "result" then rscroll = rscroll + pagesize() end
  inval()
end

function on.tabKey()
  if screen == "form" then on.arrowDown()
  elseif screen == "result" and rviews then screen = "graph"; vidx = 1; lscroll = 0; inval()
  elseif screen == "graph" then
    if vidx < #rviews then vidx = vidx + 1; lscroll = 0 else screen = "result" end
    inval()
  end
end
function on.backtabKey() if screen == "form" then on.arrowUp() end end

function on.enterKey()
  if screen == "home" then openmodule(hsel)
  elseif screen == "form" then compute()
  elseif screen == "graph" then screen = "result"; inval()
  else
    if MODULES[cur].help then screen = "home" else screen = "form" end
    inval()
  end
end

function on.returnKey() on.enterKey() end

function on.escapeKey()
  if screen == "graph" then screen = "result"; inval(); return end
  if screen == "result" and not MODULES[cur].help then screen = "form"
  else screen = "home" end
  inval()
end

function on.charIn(ch)
  if screen == "home" then
    local d = tonumber(ch)
    if d == 0 then d = 10 end
    if d and d >= 1 and d <= #MODULES then hsel = d; openmodule(d) end
    return
  end
  -- en resultados: tecla c = guardar como caso A y volver a editar
  if screen == "result" and (ch == "c" or ch == "C") and rmet and not rerror then
    local copy = {}
    for k_, v_ in pairs(VALS[cur]) do copy[k_] = v_ end
    cmpbase = {mod = cur, vals = copy, met = rmet}
    screen = "form"; cpos = #getval()
    inval()
    return
  end
  -- en resultados: tecla numerica = enviar datos a otro modulo
  if screen == "result" and rchain and rchain[ch] then
    local c = rchain[ch]
    local t = c.target
    for k, val in pairs(c.vals) do VALS[t][k] = val end
    cur = t; fsel = 1; hsel = t
    compute()
    return
  end
  if screen ~= "form" then
    if ch == "+" or ch == "-" then rscroll = rscroll + (ch == "+" and pagesize() or -pagesize()); inval() end
    return
  end
  local add = normalize(ch)
  local s = getval()
  s = s:sub(1, cpos) .. add .. s:sub(cpos + 1)
  cpos = cpos + #add
  setval(s)
  inval()
end

function on.backspaceKey()
  if screen == "form" then
    local s = getval()
    if cpos > 0 then
      local p = prevchar(s, cpos)
      setval(s:sub(1, p) .. s:sub(cpos + 1))
      cpos = p
    end
  elseif screen == "result" then on.escapeKey(); return
  end
  inval()
end

function on.deleteKey()
  if screen == "form" then
    local s = getval()
    if cpos < #s then
      local q = nextchar(s, cpos)
      setval(s:sub(1, cpos) .. s:sub(q + 1))
    end
  end
  inval()
end

function on.clearKey()
  if screen == "form" then setval(""); cpos = 0 end
  inval()
end

function on.mouseDown(x, y)
  if screen == "home" then
    local i = HOME.scroll + math.floor((y - HOME.y0) / HOME.rh) + 1
    if y >= HOME.y0 and i >= 1 and i <= #MODULES and i <= HOME.scroll + HOME.vis then
      if i == hsel then openmodule(i) else hsel = i; inval() end
    end
  elseif screen == "form" then
    local top, hintH, vis = formlayout()
    local i = fscroll + math.floor((y - top) / FH) + 1
    if y >= top and i >= 1 and i <= #MODULES[cur].fields and i <= fscroll + vis then
      fsel = i; cpos = #getval(); inval()
    end
  end
end

--------------------------------------------------------------------------
-- guardar / restaurar datos con el documento
--------------------------------------------------------------------------
local STATE_VER = 3
function on.save()
  return {vals = VALS, hsel = hsel, ver = STATE_VER}
end

function on.restore(st)
  -- datos de una version anterior (otro orden de modulos): se ignoran
  if type(st) == "table" and st.ver == STATE_VER and type(st.vals) == "table" then
    for i, t in pairs(st.vals) do
      if type(t) == "table" and VALS[i] then
        for k, v in pairs(t) do if type(v) == "string" then VALS[i][k] = v end end
      end
    end
    hsel = tonumber(st.hsel) or 1
  end
  inval()
end

--------------------------------------------------------------------------
-- menu (tecla menu)
--------------------------------------------------------------------------
if toolpalette then
  local items = {"Modulos"}
  for i, m in ipairs(MODULES) do
    items[#items + 1] = {i .. " " .. m.title, function() openmodule(i) end}
  end
  local function clearall()
    if cur then
      for _, f in ipairs(MODULES[cur].fields) do VALS[cur][f[1]] = f[4] or "" end
      cpos = #getval(); screen = MODULES[cur].help and "home" or "form"; inval()
    end
  end
  toolpalette.register({
    items,
    {"Editar",
      {"Vaciar campo (clear)", function() if screen == "form" then setval(""); cpos = 0; inval() end end},
      {"Reiniciar modulo actual", clearall},
      {"Volver al inicio", function() screen = "home"; inval() end},
      {"Quitar comparacion (caso A)", function() cmpbase = nil; inval() end},
    },
  })
end
