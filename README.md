# simulador-control-automatico

Asistente de cálculo para el examen de **Control Automático** (APN) en la
**TI-Nspire CX II CAS**. Es un documento `.tns` con una aplicación Lua: usted
escribe los datos que entrega el enunciado y la aplicación calcula los valores
que faltan y muestra el desarrollo paso a paso. Así se evitan errores al
expandir polinomios, ingresar números o despejar fórmulas.

Archivo para la calculadora: **`dist/ControlAuto.tns`**

## Instalación

1. Descargue `dist/ControlAuto.tns`.
2. Con **TI-Nspire CX Student Software** o **TI Connect CE**, arrastre el archivo
   a la calculadora (por ejemplo, a *Mis documentos*).
3. En la calculadora: `home` → *Mis documentos* → `ControlAuto`.

Requiere un SO 3.2 o superior. Cualquier TI-Nspire CX II cumple este requisito.

## Uso

| Tecla | Acción |
|---|---|
| `1`–`8` o flechas + `enter` | abrir un módulo |
| flechas ▲▼ / `tab` | cambiar de campo |
| flechas ◄► | mover el cursor dentro del campo |
| `del` / `clear` | borrar un carácter / vaciar el campo |
| `enter` | **CALCULAR** |
| ▲▼ (◄► = página) | desplazarse por los resultados |
| `esc` | volver |
| `menu` | elegir módulo, reiniciar campos |

Los campos que no conoce se dejan **en blanco**. Los datos quedan guardados al
guardar el documento.

### Cómo escribir G(s)

Use `s` como variable. La multiplicación implícita está permitida:

```
4/((s+1)(s+2))      2/(s^2+2s+17)       (s+9)/s
k/(s(s+1)(s+10))    10s/(s^2+4s+8)      1.309(1+1/(5s))
```

- El **retardo** `e^(-Ls)` **no** se escribe en G(s): va en el campo `L` / `retardo L`.
- `k` es un **parámetro libre**. Si deja el campo `k` vacío, el módulo de Routh
  entrega el **rango de k estable** y el módulo de error **despeja k** para un
  error deseado.
- También acepta `pi`, `sqrt()`, `ln()`, `log()`, `exp()`, `sin()`, `cos()` y
  `tan()` en los campos numéricos (ej.: `MF = 60` o `pi/3` según el campo).

## Módulos

| # | Módulo | Calcula | Temario |
|---|---|---|---|
| 1 | Especificaciones 2º orden | ζ ↔ Mp, σ, ωn, ωd, ts (2 % y 5 %), tp, tr, polos y polinomio deseado | UT2 |
| 2 | Polos y respuesta (LA/LC) | G0, GLC = G0/(1+G0), ec. característica, polos y ceros, estabilidad, y(∞), polos dominantes (regla de 5×), modelo reducido con la misma ganancia estática, ts/Mp/tp/tr | UT2, UT3 |
| 3 | Error estado estacionario | tipo N, Kp/Kv/Ka, ep/ev/ea con amplitud A y verificación de estabilidad; **despeje de k** para un error objetivo | UT3 |
| 4 | Routh-Hurwitz | tabla de Routh (casos ε y fila nula con polinomio auxiliar), cambios de signo, raíces; **rango de k estable** y frecuencia de cruce jω | UT2, UT3 |
| 5 | Frecuencia y márgenes | ωcg, MF, ωcf, MG (lineal y dB), ganancia crítica, Ku/Pu, estabilidad; \|G(jω)\|, ∠G(jω) y salida senoidal permanente. Admite retardo exacto | UT4 |
| 6 | Localización de polos | polos deseados, ∠Gp(pd), deficiencia angular, **elección de estructura mínima** (P/PD/PI/PID) y su justificación, cero(s), Kc, formas ideal y paralela, verificación de polos de LC, error estacionario. Opción de PID algebraico con un tercer polo fijo | UT5 |
| 7 | Sintonía empírica PI/PID | ZN lazo abierto, ZN lazo cerrado (Ku dato, relé `4d/(πa)` o calculado del modelo), **Ho-Hang-Cao** (ωp, Kp, Ti, ωg, Td = τ2) y verificación de MF/MG con el **retardo exacto**; efecto de multiplicar Kp (ej. ×2) | UT7 |
| 8 | Ayuda | teclas y sintaxis | |

Las fórmulas son las del **formulario del examen**: `ts,2% = 4/(ζωn)`,
`ζ = −ln(Mp/100)/√(π²+ln²(Mp/100))`, Ho-Hang-Cao con MF en radianes y MG lineal,
ZN-PI `Kp = 0.9τ/(KL)`, `Ti = 3.33L`, `Kp = 0.45Ku`, `Ti = Pu/1.2` y `Ku ≈ 4d/(πa)`.

## Ejemplos (última evaluación)

**Pregunta 1**: `Gp = 4/((s+1)(s+2))`, ts(2 %) ≈ 2 s, Mp ≈ 16,3 %, sin seguimiento exacto.
Módulo 6 → `Gp = 4/((s+1)(s+2))`, `ts = 2`, `Mp = 16.3`, `error nulo = 0`:

- ζ = 0,5; ωn = 4; **pd = −2 ± j3,464**
- ∠Gp(pd) = −196,1° ⇒ φd = 16,1°. Como 0° < φd < 90°, la estructura mínima es **PD**
- **Gc = 0,25 (s + 14)** (Kp = 3,5; Td = 0,0714 s); LC: s² + 4s + 16 ✓
- Escalón unitario: Kp = 7 ⇒ **ess = 0,125 (12,5 %)**, y(∞) = 0,875

**Pregunta 2**: `Gp = 2e^(−1.5s)/(9s+1)`, MF = 60°, MG = 3.
Módulo 7 → `K = 2`, `τ = 9`, `L = 1.5`, `MG = 3`, `MF = 60`, `factor Kp = 2`:

- Ho-Hang-Cao: ωp = 1,047 rad/s; **Kp = 1,571; Ti = 9 s**; ωg = 0,349 rad/s
- Verificación con el retardo exacto: ωcg = 0,349, **MF = 60°**; ωcf = 1,047, **MG = 3**
- Con Kp × 2 (Ti fijo): MF = 30° y MG = 1,5. La respuesta es más rápida, pero la estabilidad relativa es menor

## Estructura del repositorio

```
src/core.lua      núcleo matemático: parser de G(s), polinomios, raíces,
                  Routh, márgenes, errores (Lua 5.1 puro)
src/modules.lua   módulos de cálculo (texto de resultados)
src/ui.lua        interfaz TI-Nspire (API Lua 2.0)
tools/build.py    une los fuentes y genera dist/ControlAuto.lua y .tns
tools/make_tns.py empaquetador .tns (port en Python del formato de Luna)
tests/            regresión contra las respuestas de las guías 1-7
dist/             ControlAuto.tns (para la calculadora) y ControlAuto.lua
```

Para regenerar el documento y ejecutar las pruebas:

```
pip install lupa cryptography
python3 tools/build.py
python3 tests/test_guias.py      # 46 casos de las guías
python3 tests/ui_smoke.py        # interfaz con API Nspire simulada
```

## Notas

- Las raíces y los márgenes se calculan numéricamente (Durand-Kerner y barrido
  de frecuencia con bisección), con 5 cifras significativas.
- La estabilidad de un lazo con retardo se evalúa con el criterio de márgenes
  (MF > 0 y MG > 1), como en el curso.
- Revise que el uso de la aplicación esté permitido en su evaluación.
