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
| `1`–`9`, `0` (módulo 10) o flechas + `enter` | abrir un módulo |
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
| 2 | Respuesta y(t) | Y(s) = G·R (escalón, rampa o impulso de magnitud A), **fracciones parciales** (incluye polos complejos y repetidos), **y(t)**, modos (converge/diverge, τ), modo más veloz y modo dominante, estabilidad, y(∞) | UT2 |
| 3 | Polos y respuesta (LA/LC) | G0, GLC = G0/(1+G0), ec. característica, polos y ceros, estabilidad, y(∞), polos dominantes (regla de 5×), modelo reducido con la misma ganancia estática, ts/Mp/tp/tr. Con `k` vacío: **GLC en función de K** y ωn², 2ζωn en función de K | UT2, UT3 |
| 4 | Error estado estacionario | tipo N, Kp/Kv/Ka, ep/ev/ea con amplitud A y verificación de estabilidad; **despeje de k** para un error objetivo | UT3 |
| 5 | Routh-Hurwitz | tabla de Routh (casos ε y fila nula con polinomio auxiliar), cambios de signo, raíces; **rango de k estable** y frecuencia de cruce jω | UT2, UT3 |
| 6 | Frecuencia y márgenes | ωcg, MF, ωcf, MG (lineal y dB), ganancia crítica, Ku/Pu, estabilidad; \|G(jω)\|, ∠G(jω) y salida senoidal permanente. Admite retardo exacto | UT4 |
| 7 | Localización de polos | polos deseados, ∠Gp(pd), deficiencia angular, **elección de estructura mínima** (P/PD/PI/PID) y su justificación, cero(s), Kc, formas ideal y paralela, verificación de polos de LC, error estacionario. Opción de PID algebraico con un tercer polo fijo | UT5 |
| 8 | Sintonía empírica PI/PID | ZN lazo abierto, ZN lazo cerrado (Ku dato, relé `4d/(πa)` o calculado del modelo, también con dos constantes de tiempo), **Ho-Hang-Cao** (ωp, Kp, Ti, ωg, Td = τ2) y verificación de MF/MG con el **retardo exacto**; efecto de multiplicar Kp (ej. ×2) | UT7 |
| 9 | Identificación (curva) | K, τ y L de un modelo de primer orden con retardo a partir de la curva de reacción (método de dos puntos 28,3 %/63,2 % o lectura con tangente); ZN inmediato y envío al módulo 8 | UT7 |
| 10 | Entrenamiento | ejercicios aleatorios del tipo de las pruebas (8 tipos), corrección de sus respuestas y solución paso a paso | todas |
| 11 | Ayuda | qué módulo usar según el enunciado, glosario, teclas y sintaxis | |

Las fórmulas son las del **formulario del examen**: `ts,2% = 4/(ζωn)`,
`ζ = −ln(Mp/100)/√(π²+ln²(Mp/100))`, Ho-Hang-Cao con MF en radianes y MG lineal,
ZN-PI `Kp = 0.9τ/(KL)`, `Ti = 3.33L`, `Kp = 0.45Ku`, `Ti = Pu/1.2` y `Ku ≈ 4d/(πa)`.

## Ejemplos (última evaluación)

**Pregunta 1**: `Gp = 4/((s+1)(s+2))`, ts(2 %) ≈ 2 s, Mp ≈ 16,3 %, sin seguimiento exacto.
Módulo 7 → `Gp = 4/((s+1)(s+2))`, `ts = 2`, `Mp = 16.3`, `error nulo = 0`:

- ζ = 0,5; ωn = 4; **pd = −2 ± j3,464**
- ∠Gp(pd) = −196,1° ⇒ φd = 16,1°. Como 0° < φd < 90°, la estructura mínima es **PD**
- **Gc = 0,25 (s + 14)** (Kp = 3,5; Td = 0,0714 s); LC: s² + 4s + 16 ✓
- Escalón unitario: Kp = 7 ⇒ **ess = 0,125 (12,5 %)**, y(∞) = 0,875

**Pregunta 2**: `Gp = 2e^(−1.5s)/(9s+1)`, MF = 60°, MG = 3.
Módulo 8 → `K = 2`, `τ = 9`, `L = 1.5`, `MG = 3`, `MF = 60`, `factor Kp = 2`:

- Ho-Hang-Cao: ωp = 1,047 rad/s; **Kp = 1,571; Ti = 9 s**; ωg = 0,349 rad/s
- Verificación con el retardo exacto: ωcg = 0,349, **MF = 60°**; ωcf = 1,047, **MG = 3**
- Con Kp × 2 (Ti fijo): MF = 30° y MG = 1,5. La respuesta es más rápida, pero la estabilidad relativa es menor

### Resultados: respuestas, gráficos y envío entre módulos

- Cada resultado empieza con un recuadro **RESPUESTAS** con lo que se escribe en la hoja; debajo va el **Desarrollo** paso a paso.
- **`tab`** recorre los gráficos del resultado:
  - **y(t)**: respuesta al escalón con Mp, tp y ts **reales** de la simulación exacta (módulos 1, 2, 3, 6 y 7). En los módulos 6 y 8, la simulación incluye el retardo exacto y compara los controladores HHC y ZN.
  - **Bode** (módulo 6), con ωcg, ωcf, MF y MG marcados.
  - **Lugar de raíces** (módulo 5 con `k`), con los cruces del eje jω y su k.
  - **Plano s** (módulo 7).
- **Teclas numéricas en el resultado**: envían el controlador o el modelo a otro módulo sin reescribirlo. Por ejemplo, en el módulo 7, `3` lleva a la respuesta LC y `6` a márgenes/Bode; en el módulo 9, `8` lleva a Sintonía.
- **Entrenamiento (10)**: elija el `tipo` y el `nro` del ejercicio, presione enter para ver el enunciado, resuelva a mano, escriba r1..r3 y presione enter para corregir. Con `?` en r1 se muestra la solución.

### Gráfico del plano s (módulo 7)

Después de calcular un diseño por localización de polos, presione **`tab`** para ver el plano s:

- **Rojo:** polos (×) y ceros (○) de la planta.
- **Azul:** polos deseados pd.
- **Verde:** lo que agrega el controlador (integrador ×, cero ○) y los polos adicionales de lazo cerrado.
- Cada raíz se une con pd y su ángulo se marca con una letra griega (α, β, γ…). La leyenda explica cada uno, la suma ∠G(pd), la deficiencia φd y el ángulo φ del cero del controlador.
- Las raíces muy lejanas se dibujan en el borde, con su valor.

Teclas: ▲▼ desplazan la leyenda; `tab`, `esc` o `enter` vuelven al texto.

### Pruebas de años anteriores

| Pregunta | Módulo y datos | Resultado de la app |
|---|---|---|
| y(t) de `20/((s+2)(s+1)(s-6))`, escalón de 5 | 2: `G`, `A = 5` | y(t) = 0,2976e^(6t) − 8,333 + 14,29e^(−t) − 6,25e^(−2t) ⇒ modo e^(6t) diverge: **inestable** |
| y(t) de `2/(s^3+4s^2+s+6)`, escalón de 5 | 2: `G`, `A = 5` | polos 0,056 ± j1,207 y −4,11 ⇒ **inestable** |
| K·Gp con `2/(s(s-1))`, K = 10 | 3: `Gc = 10` | polos 0,5 ± j4,44 ⇒ inestable |
| `Gc = 3`, `Gp = 3/(s(s+9)(s+7))` | 3, 4 y 5 | polo dominante −0,148 (relación 43 ≥ 5), ts = 26,96 s; error rampa; 0 < K < 336 |
| `Gc = 9/s`, `Gp = 1/(s^2+2s+5)` | 3 | −0,0585 ± j2,186 dominantes; ts = 68,4 s; Mp = 91,9 % |
| `K·4/(s(s+4))` en función de K; ev = 0,125 | 3 con `k` vacío; 4 con `e obj = 0.125`, `señal = 2` | GLC = 4K/(s²+4s+4K); **K = 8** |
| `K·10/(s(s^2+7s+20))` | 5 y 3 (`k = 5`) | 0 < K < 14; con K = 5: −1 ± j3, −5 (relación 5); Gr = 10/(s²+2s+10); ts5 % = 3 s, tr = 0,631 s, Mp = 35,1 % |
| PI, ts = 29 s, Mp = 9 %, `5/(s^2+10s+1)` | 7: `error nulo = 1` | Gc = 0,347(s + 0,288)/s (pauta: 0,35 y 0,29) |
| ZN-PI de `1.5e^(-10s)/(200s+1)` | 8: `K = 1.5`, `τ = 200`, `L = 10` | Kp = 12, Ti = 33,3 s; MF = 21,3°, MG = 1,52 (pauta: 21,35° y 1,52) |
| ZN-LC PI de `2e^(-18s)/((102s+1)(87s+1))` | 8: `τ = 102`, `τ2 = 87`, `L = 18` | Ku = 5,59, Pu = 185,6 s ⇒ Kp = 2,52, Ti = 154,7 s; MF y MG del lazo |

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
python3 tests/test_guias.py      # 60 casos (guías y pruebas)
python3 tests/ui_smoke.py        # interfaz con API Nspire simulada
```

## Notas

- Las raíces y los márgenes se calculan numéricamente (Durand-Kerner y barrido
  de frecuencia con bisección), con 5 cifras significativas.
- La estabilidad de un lazo con retardo se evalúa con el criterio de márgenes
  (MF > 0 y MG > 1), como en el curso.
- Revise que el uso de la aplicación esté permitido en su evaluación.
