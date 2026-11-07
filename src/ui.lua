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
    rlines = res; rerror = false
  else
    local msg = tostring(res):gsub("^.-:%d+: ", "")
    rlines = {"ERROR", msg, "", "Revise los datos (esc para volver)."}
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

local function paint_home(gc)
  titlebar(gc, "Control Automatico", "APN - examen")
  local y0 = 19
  local n = #MODULES
  local rh = math.floor((H() - y0 - 12) / n)
  if rh > 24 then rh = 24 end
  for i, m in ipairs(MODULES) do
    local y = y0 + (i - 1) * rh
    if i == hsel then col(gc, C.sel); gc:fillRect(0, y, W(), rh) end
    col(gc, C.focus); gc:setFont("sansserif", "b", 10)
    gc:drawString(tostring(i), 6, y + 1, "top")
    col(gc, C.black)
    gc:drawString(m.title, 20, y + 1, "top")
    col(gc, C.grey); gc:setFont("sansserif", "r", 7)
    if rh >= 20 then gc:drawString(m.short, 20, y + rh - 11, "top") end
  end
  col(gc, C.grey); gc:setFont("sansserif", "r", 7)
  gc:drawString("1-9 o flechas + enter.  menu: modulos", 4, H() - 11, "top")
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
  titlebar(gc, cur .. ". " .. m.title, "enter: calcular")
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
  titlebar(gc, cur .. ". " .. m.title, "esc: volver")
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

function on.paint(gc)
  col(gc, C.white); gc:fillRect(0, 0, W(), H())
  if screen == "home" then paint_home(gc)
  elseif screen == "form" then paint_form(gc)
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
  else rscroll = rscroll - 1 end
  inval()
end

function on.arrowDown()
  if screen == "home" then hsel = hsel % #MODULES + 1
  elseif screen == "form" then
    local n = #MODULES[cur].fields
    fsel = fsel % n + 1; cpos = #getval()
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

function on.tabKey() if screen == "form" then on.arrowDown() end end
function on.backtabKey() if screen == "form" then on.arrowUp() end end

function on.enterKey()
  if screen == "home" then openmodule(hsel)
  elseif screen == "form" then compute()
  else
    if MODULES[cur].help then screen = "home" else screen = "form" end
    inval()
  end
end

function on.returnKey() on.enterKey() end

function on.escapeKey()
  if screen == "result" and not MODULES[cur].help then screen = "form"
  else screen = "home" end
  inval()
end

function on.charIn(ch)
  if screen == "home" then
    local d = tonumber(ch)
    if d and d >= 1 and d <= #MODULES then hsel = d; openmodule(d) end
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
    local y0 = 19
    local rh = math.floor((H() - y0 - 12) / #MODULES)
    if rh > 24 then rh = 24 end
    local i = math.floor((y - y0) / rh) + 1
    if i >= 1 and i <= #MODULES then
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
function on.save()
  return {vals = VALS, hsel = hsel}
end

function on.restore(st)
  if type(st) == "table" and type(st.vals) == "table" then
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
    },
  })
end
