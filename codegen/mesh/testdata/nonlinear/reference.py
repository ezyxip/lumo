"""Independent reduced calculation for a binary three-stage Rating column.

Uses scalar bisection and overall/reboiler balances, not the full MESH vector,
its Jacobian, or generated C. Standard Python only. Run this file to regenerate
reference.json next to it. This is a numerical regression reference for the
chosen ideal model, not experimental validation or an external simulator.
"""
import hashlib
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
raw = (ROOT / "case.json").read_bytes()
case = json.loads(raw)
cfg, op = case["compiler_input"], case["operating"]
props = [c["properties"] for c in cfg["components"]]
pressure = cfg["column"]["pressure"]["value_Pa"]

assert cfg["column"]["stages"] == 3
assert cfg["column"]["feeds"] == [{"id": "main_feed", "stage": 2}]
assert cfg["column"]["condenser"]["type"] == "total"
assert cfg["column"]["reboiler"]["type"] == "equilibrium"
assert cfg["column"]["pressure"]["type"] == "constant"
assert cfg["thermo"] == {"phase_equilibrium": "antoine_raoult", "enthalpy": "ideal"}
assert len(props) == 2 and op["feed_state"] == "saturated_liquid"
assert op["feed_mole_fractions"] == [0.5, 0.5]
assert op["distillate_mol_s"] == op["feed_flow_mol_s"]/2
for p in props:
    assert p["source"] == "inline" and p["enthalpy"]["type"] == "constant_cp"
    assert p["antoine"]["pressure_unit"] == "mmHg" and p["antoine"]["temperature_unit"] == "C"
    assert p["temperature_range_K"] == {"min": 350, "max": 390}


def bisect(function, lo, hi):
    flo, fhi = function(lo), function(hi)
    if flo * fhi > 0:
        raise ValueError("root is not bracketed")
    for _ in range(100):
        mid = (lo + hi) / 2
        fm = function(mid)
        if fm == 0 or mid == lo or mid == hi:
            return mid
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
    raise ValueError("bisection did not converge")


def phase(x):
    def k(i, temperature):
        a = props[i]["antoine"]
        return 133.3224 * 10 ** (a["A"] - a["B"] / (a["C"] + temperature - 273.15)) / pressure
    temperature = bisect(lambda t: x * k(0, t) + (1-x) * k(1, t) - 1, 350, 390)
    return temperature, x * k(0, temperature)


def enthalpy(x, temperature, vapor=False):
    terms = []
    for p in props:
        h = p["enthalpy"]
        dt = temperature - h["reference_temperature_K"]
        terms.append(h["vaporization_enthalpy_at_reference_J_mol"] + h["vapor_cp_J_mol_K"] * dt
                     if vapor else h["liquid_cp_J_mol_K"] * dt)
    return x * terms[0] + (1-x) * terms[1]


F, z = op["feed_flow_mol_s"], op["feed_mole_fractions"][0]
D, R = op["distillate_mol_s"], op["reflux_ratio"]
B, reflux, vapor_top = F-D, R*D, (R+1)*D
feed_temperature, _ = phase(z)
feed_enthalpy = enthalpy(z, feed_temperature)


def reduced(top_x):
    # Condenser fixes the vapor composition leaving the feed tray.
    feed_x = bisect(lambda x: phase(x)[1] - top_x, 0, 1)
    bottom_x = (F*z - D*top_x)/B  # overall light-component balance
    top_t, formal_top_y = phase(top_x)
    feed_t, feed_y = phase(feed_x)
    bottom_t, bottom_y = phase(bottom_x)
    # Reboiler total + light-component balances determine its boilup.
    denominator = bottom_y - feed_x
    if denominator <= 0:
        return None
    vapor_bottom = B*(feed_x-bottom_x)/denominator
    if vapor_bottom <= 0:
        return None
    liquid_feed = B + vapor_bottom
    hl = [enthalpy(x, t) for x, t in zip((top_x, feed_x, bottom_x), (top_t, feed_t, bottom_t))]
    hv = [enthalpy(y, t, True) for y, t in zip((formal_top_y, feed_y, bottom_y), (top_t, feed_t, bottom_t))]
    balance = reflux*hl[0] + vapor_bottom*hv[2] + F*feed_enthalpy - liquid_feed*hl[1] - vapor_top*hv[1]
    qd = vapor_top*hv[1] - (reflux+D)*hl[0]
    qw = B*hl[2] + vapor_bottom*hv[2] - liquid_feed*hl[1]
    state = [top_x, 1-top_x, formal_top_y, 1-formal_top_y, top_t, reflux, qd,
             feed_x, 1-feed_x, feed_y, 1-feed_y, feed_t, liquid_feed, vapor_top,
             bottom_x, 1-bottom_x, bottom_y, 1-bottom_y, bottom_t, B, vapor_bottom,
             qw]
    return balance, state


previous = None
bracket = None
for n in range(1, 1000):
    top_x = 0.5 + 0.499*n/1000
    result = reduced(top_x)
    if result is None:
        previous = None
        continue
    if previous is not None and previous[1]*result[0] <= 0:
        bracket = previous[0], top_x
        break
    previous = top_x, result[0]
if bracket is None:
    raise ValueError("no physical reduced root found")
root = bisect(lambda x: reduced(x)[0], *bracket)
balance, state = reduced(root)
assert abs(balance) < 1e-5, balance
assert all(math.isfinite(v) for v in state)
assert all(0 < state[j*7] < 1 for j in range(3))
reference = {
    "method": "Reduced three-stage balances and scalar bisection, reference.py",
    "input_sha256": hashlib.sha256(raw).hexdigest(),
    "state": state,
    "operating": {
        "feed_flow_mol_s": F,
        "feed_mole_fractions": op["feed_mole_fractions"],
        "feed_enthalpy_J_mol": feed_enthalpy,
        "distillate_mol_s": D,
        "reflux_ratio": R,
    },
}
(ROOT / "reference.json").write_text(json.dumps(reference, indent=2) + "\n")
print(f"Reduced solution xD={state[0]:.12f}, xB={state[14]:.12f}, feed energy residual={balance:.3g} W")
