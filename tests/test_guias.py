"""Regresion: compara salidas con respuestas de las guias del curso."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from run_cases import run

CASES = [
  # (modulo, campos, textos esperados en la salida)
  ("spec", dict(Mp="16.3", ts="2"), ["0.50004", "-2 ± j3.4637"]),
  ("spec", dict(Mp="25.38", wn="5"), ["ζ = 0.4"]),
  ("polos", dict(Gc="9", Gp="7/(s(s+10)(s+1))"), ["-0.19149 ± j2.4284", "-10.617", "20.889"]),
  ("polos", dict(Gp="40/(s^2+10s+16)"), ["-5 ± j5.5678", "Mp = 100e^(-ζπ/√(1-ζ²)) = 5.95", "0.41"]),
  ("polos", dict(Gp="10/(s^3+7s^2+13s+15)"), ["1.8324 / (s^2 + 1.5427s + 4.581)", "29.708"]),
  ("polos", dict(Gp="40/(s^3+12s^2+41s+42)"), ["4.363", "NO conviene"]),
  ("polos", dict(Gp="25/(s^2+4s+25)", lazo="0"), ["25.38", "tp = π/ωd = 0.6855"]),
  ("error", dict(Gp="63/(s(s+10)(s+1))"), ["ev = A/Kv = 0.15873"]),
  ("error", dict(Gp="7/(s^2(s^2+2s+1))"), ["ea = A/Ka = 0.14286"]),
  ("error", dict(Gp="10/(s^2+4s+8)"), ["ep = A/(1+Kp) = 0.44444"]),
  ("error", dict(Gp="10s/(s^2+4s+8)", Gc="1/s^2"), ["ev = A/Kv = 0.8"]),
  ("error", dict(Gp="2k/(s(s+1))", eobj="0.1", senal="2"), ["k = 5 "]),
  ("error", dict(Gp="10/((s+2)(s+k))", eobj="0.1"), ["k = 0.555556"]),
  ("routh", dict(P="s^3+4s^2+5s+k"), ["0 < k < 20", "±j2.2361"]),
  ("routh", dict(Gp="1/(s^3+4s^2+s+2)"), ["-2 < k < 2"]),
  ("routh", dict(Gp="(s+1)/(s^2+8s+1)"), ["-1 < k <"]),
  ("routh", dict(Gp="1/(s(s^2+5s+1))"), ["0 < k < 5"]),
  ("routh", dict(Gp="7/(s(s+1)(s+10))"), ["0 < k < 15.7143", "±j3.1623"]),
  ("routh", dict(P="s^4+2s^3+3s^2+6s+4"), ["Cambios de signo: 2"]),
  ("routh", dict(P="s^4+5s^2+4"), ["CRITICAMENTE"]),
  ("freq", dict(Gp="1/(s(s+1))"), ["0.78615", "MF = 51.827"]),
  ("freq", dict(Gp="20/((s+1)^2(s+6))"), ["MF = 53.574", "MG = 4.9 "]),
  ("freq", dict(Gp="10/(s^2+5s+1)"), ["1.925", "74.29"]),
  ("freq", dict(Gp="20/(s+1)"), ["19.975", "92.866"]),
  ("freq", dict(Gc="(s+9)/s", Gp="1/(s^2+6s+5)"), ["MF = 35.149", "MG = 10 "]),
  ("freq", dict(Gp="2/(s+2)", w="2", A="3"), ["yss = 2.1213 sin(2t - 45°)"]),
  ("place", dict(Gp="4/((s+1)(s+2))", ts="2", Mp="16.3"), ["PD", "0.25(s + 13.997)", "ep = A/(1+Kp) = 0.125"]),
  ("place", dict(Gp="1/(s^2-2s+2)", ts="3", Mp="28", crit="5"), ["4(s + 1.2727)"]),
  ("place", dict(Gp="3/(s^2+8s+2)", ts="5", Mp="15", enulo="1"), ["PI", "3.545(s + 1.4413)/s", "-6.4"]),
  ("place", dict(Gp="2/(s^2+7s+1)", ts="19", Mp="11", enulo="1"), ["0.95", "(s + 0.463"]),
  ("place", dict(Gp="9/(s^2+10s+8)", ts="3", Mp="7", enulo="1"), ["1.75", "(s + 1.97"]),
  ("place", dict(Gp="0.5/(s-1)", ts="4", Mp="20.8", enulo="1"), ["6(s + 1.6677)/s"]),
  ("place", dict(Gp="36/(s^2+2s+2)", ts="2", Mp="16", enulo="1"), ["PID", "0.10675(s + 2.7488)", "Ti = 2/a = 0.72759"]),
  ("place", dict(Gp="25/(s^2+2s+14)", ts="1", Mp="28", enulo="1"), ["0.391", "(s + 6.62"]),
  ("place", dict(Gp="34/(s^2+9s+23)", ts="1", Mp="6", enulo="1"), ["0.112", "(s + 6.73"]),
  ("place", dict(Gp="4/(s^2+9s+3)", ts="1", Mp="18", enulo="1"), ["1.91", "(s + 8.87"]),
  ("place", dict(Gp="1/(s^2+1)", ts="4", Mp="16.3", enulo="1", p3="5"), ["Kd = (7-0)/1 = 7", "-5"]),
  ("emp", dict(K="2", tau="10", L="2"), ["Kp = 0.9τ/(KL) = 2.25", "Ti = 3.33L = 6.66", "Td = 0.5L = 1"]),
  ("emp", dict(Ku="13", Pu="6.7"), ["5.85", "5.5833", "0.8375"]),
  ("emp", dict(K="15.1", tau="16.9", L="2.42"), ["0.41623", "MF = 27.358", "MG = 1.603"]),
  ("emp", dict(K="1", tau="0", L="0"), ["ERROR"]),
  ("emp", dict(d="2", a="0.5", Pu="6"), ["5.093", "2.2918", "Ti = Pu/1.2 = 5"]),
  ("emp", dict(K="1.5", tau="8", L="2", MG="3", MF="60"), ["0.7854", "1.3963", "= 8 s", "0.2618"]),
  ("emp", dict(K="3", tau="12", L="1.5", MG="4", MF="55"), ["1.0084", "6.32"]),
  ("emp", dict(K="3", tau="8", L="1.5", tau2="2", MG="3", MF="45"), ["0.98175", "0.87266", "4.03", "Td = τ2 = 2"]),
  ("emp", dict(K="2", tau="5", L="1", MG="3", MF="60", fKp="2"), ["MF: 60° → 30°", "MG: 3 → 1.5"]),
]

def main():
    bad = 0
    for mod, f, exp in CASES:
        out = "\n".join(run(mod, **f))
        miss = [e for e in exp if e not in out]
        if miss:
            bad += 1
            print("FALLA", mod, f, "->", miss)
            print("   " + out.replace("\n", "\n   "))
    print("%d/%d casos OK" % (len(CASES) - bad, len(CASES)))
    sys.exit(1 if bad else 0)

if __name__ == "__main__":
    main()
