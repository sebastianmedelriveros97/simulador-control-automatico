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
        merged.views = res.views; merged.chain = res.chain; merged.met = res.met; merged.cansim = res.cansim
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
  titlebar(gc, "Control Automatico", "v14 - APN")
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
  for i, c in ipairs(V.curves) do if c.name then leg[#leg + 1] = {c.name, c.col or CURVE[i]} end end
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
    col(gc, c.col or CURVE[i])
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
    local m2 = V.marks2
    if m2 and m2.tp and m2.ymax then
      col(gc, {0x00, 0x96, 0x20}); gc:fillArc(X(m2.tp) - 3, Y(m2.ymax) - 3, 6, 6, 0, 360)
      gc:drawString("Mp=" .. fmt(m2.Mp, 3) .. "%", math.min(X(m2.tp) + 4, W() - 60), Y(m2.ymax) + 3, "top")
    end
    if m.ts then
      col(gc, V.marks2 and C.grey or CURVE[3]); dashed(gc, true)
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
  -- en resultados: tecla s = simular la respuesta real (curva verde)
  if screen == "result" and (ch == "s" or ch == "S") and rlines and rlines.cansim and not rerror then
    VALS[cur]._sim = "1"
    local keepcmp = cmpbase
    compute()
    VALS[cur]._sim = nil
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
