"""Independent binary stage shooting reference, Python standard library only.

Equations and sources are documented in README.md. No Lumo imports, generated
C, full MESH vector, Jacobian or Newton iteration are used.
"""
import argparse
import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def bisect(fn, lo, hi):
    flo, fhi = fn(lo), fn(hi)
    if not math.isfinite(flo + fhi) or flo * fhi > 0:
        raise ValueError("root not bracketed")
    for _ in range(100):
        mid = (lo + hi) / 2
        fm = fn(mid)
        if fm == 0 or mid == lo or mid == hi:
            return mid
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
    raise ValueError("bisection did not converge")


def calculate():
    cfg = json.loads((ROOT / 'compiler_input.json').read_text())
    model = json.loads((ROOT / 'model_input.json').read_text())
    assert cfg['column']['stages'] == 5
    assert cfg['column']['feeds'] == [{'id': 'main_feed', 'stage': 3}]
    assert cfg['column']['condenser']['type'] == 'total'
    assert cfg['column']['reboiler']['type'] == 'equilibrium'
    assert cfg['column']['pressure']['type'] == 'constant'
    assert cfg['thermo'] == {'phase_equilibrium': 'antoine_raoult', 'enthalpy': 'ideal'}
    props = [c['properties'] for c in cfg['components']]
    assert [c['id'] for c in cfg['components']] == ['benzene', 'toluene']
    for p in props:
        assert p['source'] == 'inline'
        assert p['antoine']['pressure_unit'] == 'mmHg'
        assert p['antoine']['temperature_unit'] == 'C'
        assert p['enthalpy']['type'] == 'constant_cp'
    feed = model['feeds']['main_feed']
    specs = model['specifications']
    assert model['schema_version'] == cfg['schema_version'] == 1
    assert set(model['feeds']) == {'main_feed'}
    assert feed['state'] == {'type': 'saturated_liquid'}
    assert feed['flow']['unit'] == specs['distillate_flow']['unit'] == 'mol_s'
    assert feed['composition']['basis'] == specs['reflux_ratio']['basis'] == 'mole'
    assert specs['type'] == 'rating'
    fractions = feed['composition']['fractions']
    assert set(fractions) == {'benzene', 'toluene'}
    assert abs(sum(fractions.values()) - 1) < 1e-14
    F, z = feed['flow']['value'], fractions['benzene']
    D, R = specs['distillate_flow']['value'], specs['reflux_ratio']['value']
    B = F - D
    assert 0 < D < F and R > 0 and 0 < z < 1
    pressure = cfg['column']['pressure']['value_Pa']
    tmin = max(p['temperature_range_K']['min'] for p in props)
    tmax = min(p['temperature_range_K']['max'] for p in props)

    def k(i, t):
        a = props[i]['antoine']
        return 133.3224 * 10**(a['A'] - a['B'] / (a['C'] + t - 273.15)) / pressure

    def h(fraction, t, vapor=False):
        values = []
        for p in props:
            e = p['enthalpy']
            dt = t - e['reference_temperature_K']
            values.append(e['vaporization_enthalpy_at_reference_J_mol'] + e['vapor_cp_J_mol_K']*dt
                          if vapor else e['liquid_cp_J_mol_K']*dt)
        return fraction*values[0] + (1-fraction)*values[1]

    @lru_cache(maxsize=8192)
    def liquid(x):
        t = bisect(lambda t: x*k(0,t)+(1-x)*k(1,t)-1, tmin, tmax)
        y = x*k(0,t)
        return {'T_K': t, 'x': [x,1-x], 'y': [y,1-y],
                'hL_J_mol': h(x,t), 'hV_J_mol': h(y,t,True)}

    def vapor(y):
        t = bisect(lambda t: y/k(0,t)+(1-y)/k(1,t)-1, tmin, tmax)
        x = y/k(0,t)
        return {'T_K': t, 'x': [x,1-x], 'y': [y,1-y],
                'hL_J_mol': h(x,t), 'hV_J_mol': h(y,t,True)}

    hf = liquid(z)['hL_J_mol']

    def shoot(xd):
        top = liquid(xd) | {'L_mol_s': R*D, 'V_mol_s': None}
        current = vapor(xd) | {'V_mol_s': (R+1)*D}
        stages = [top]
        for j in range(2,5):
            prev = stages[-1]
            fj = F if j == 3 else 0
            v = current['V_mol_s']
            x, y = current['x'][0], current['y'][0]
            c = prev['L_mol_s'] + fj - v
            a = c*x + v*y - prev['L_mol_s']*prev['x'][0] - fj*z
            b = c*current['hL_J_mol'] + v*current['hV_J_mol'] - prev['L_mol_s']*prev['hL_J_mol'] - fj*hf

            def next_stage(yn):
                vn = a/(yn-x)
                return vapor(yn) | {'V_mol_s': vn}

            def energy(yn):
                nxt = next_stage(yn)
                return a*(nxt['hV_J_mol'] - current['hL_J_mol']) - b*(yn-x)

            # Positive boilup determines which side of x can contain the root.
            lo, hi = (x+1e-12, 1.0) if a > 0 else (0.0, x-1e-12)
            yn = bisect(energy,lo,hi)
            nxt = next_stage(yn)
            current['L_mol_s'] = nxt['V_mol_s'] + c
            if min(nxt['V_mol_s'], current['L_mol_s']) <= 0:
                raise ValueError('nonphysical flow')
            stages.append(current)
            current = nxt
        current['L_mol_s'] = B
        stages.append(current)
        closure = F*z-D*xd-B*current['x'][0]
        return closure, stages

    bracket = None
    previous = None
    for i in range(1,1000):
        xd = z+(1-z)*i/1000
        try:
            value, _ = shoot(xd)
        except ValueError:
            previous = None
            continue
        if previous and previous[1]*value <= 0:
            bracket = previous[0], xd
            break
        previous = xd,value
    if bracket is None:
        raise ValueError('no physical shooting root found')
    xd = bisect(lambda x: shoot(x)[0], *bracket)
    closure, stages = shoot(xd)
    qd = stages[1]['V_mol_s']*stages[1]['hV_J_mol']-(R+1)*D*stages[0]['hL_J_mol']
    qw = B*stages[-1]['hL_J_mol']+stages[-1]['V_mol_s']*stages[-1]['hV_J_mol']-stages[-2]['L_mol_s']*stages[-2]['hL_J_mol']
    # Direct local balance checks, separate from the shooting closure.
    mass, component, energy = [], [], []
    for j in range(1,4):
        p,s,n = stages[j-1:j+2]
        fj = F if j == 2 else 0
        mass.append(p['L_mol_s']+n['V_mol_s']+fj-s['L_mol_s']-s['V_mol_s'])
        component.append(p['L_mol_s']*p['x'][0]+n['V_mol_s']*n['y'][0]+fj*z-s['L_mol_s']*s['x'][0]-s['V_mol_s']*s['y'][0])
        energy.append(p['L_mol_s']*p['hL_J_mol']+n['V_mol_s']*n['hV_J_mol']+fj*hf-s['L_mol_s']*s['hL_J_mol']-s['V_mol_s']*s['hV_J_mol'])
    p,s = stages[-2:]
    mass.append(p['L_mol_s']-B-s['V_mol_s'])
    component.append(p['L_mol_s']*p['x'][0]-B*s['x'][0]-s['V_mol_s']*s['y'][0])
    global_energy = F*hf+qw-qd-D*stages[0]['hL_J_mol']-B*s['hL_J_mol']
    mass.append(stages[1]['V_mol_s']-(R+1)*D)
    component.append(stages[1]['V_mol_s']*stages[1]['y'][0]-(R+1)*D*xd)
    equilibrium = [s['y'][i]-k(i,s['T_K'])*s['x'][i] for s in stages for i in range(2)]
    checks = {'max_equilibrium_residual': max(map(abs,equilibrium)),
              'max_total_balance_mol_s': max(map(abs,mass)),
              'max_component_balance_mol_s': max(map(abs,component)),
              'max_tray_energy_balance_W': max(map(abs,energy)),
              'overall_component_balance_mol_s': closure,
              'overall_energy_balance_W': global_energy}
    assert checks['max_equilibrium_residual'] < 1e-12
    assert checks['max_total_balance_mol_s'] < 1e-9
    assert checks['max_component_balance_mol_s'] < 1e-9
    assert checks['max_tray_energy_balance_W'] < 1e-5
    assert abs(closure) < 1e-9 and abs(global_energy) < 1e-5
    for i,s in enumerate(stages,1):
        s['stage'] = i
        assert all(math.isfinite(v) for key,v in s.items() if isinstance(v,(int,float)))
        assert all(0 <= v <= 1 for v in s['x']+s['y'])
        assert abs(sum(s['y'])-1) < 1e-14
    return {'method': 'Binary stage shooting with scalar bisection, calculate.py',
            'input_sha256': {n:hashlib.sha256((ROOT/n).read_bytes()).hexdigest() for n in ['compiler_input.json','model_input.json']},
            'component_ids': ['benzene','toluene'], 'feed_enthalpy_J_mol': hf,
            'distillate_mol_s': D, 'bottoms_mol_s': B,
            'condenser_heat_removed_W': qd, 'reboiler_heat_added_W': qw,
            'stages': stages, 'checks': checks}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', action='store_true', help='verify the saved reference without overwriting it')
    args = parser.parse_args()
    result = calculate()
    path = ROOT/'reference.json'
    if args.check:
        saved = json.loads(path.read_text())
        def compare(actual, expected):
            if isinstance(actual, dict):
                assert actual.keys() == expected.keys()
                for key in actual:
                    compare(actual[key], expected[key])
            elif isinstance(actual, list):
                assert len(actual) == len(expected)
                for a,b in zip(actual, expected):
                    compare(a,b)
            elif isinstance(actual, (int,float)):
                assert math.isclose(actual,expected,rel_tol=1e-12,abs_tol=1e-7), (actual,expected)
            else:
                assert actual == expected, (actual,expected)
        compare(result,saved)
    else:
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result['checks'],indent=2))
