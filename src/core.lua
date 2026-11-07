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
  arrow="\226\134\146", le="\226\137\164", ge="\226\137\165", eps="\206\181",
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
function proots(pin)
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
    for k = 1, n do
      local ang = 2 * pi * (k - 1) / n + 0.4
      zr[k] = R * math.cos(ang)
      zi[k] = R * math.sin(ang)
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
      if maxstep < 1e-15 then break end
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
        out[#out + 1] = fmt(r[1], d) .. " " .. U.pm .. " j" .. fmt(r[2], d)
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
