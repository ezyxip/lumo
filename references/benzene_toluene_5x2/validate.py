"""Validate the real Lumo CLI and generated C99 against the saved 5x2 reference.

Python 3, Go and a C99 compiler are required. No numerical Python dependencies.
This script never regenerates or updates reference.json.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
TOLERANCES = {'T_K': 1e-6, 'composition': 1e-8, 'flow_mol_s': 1e-6, 'heat_W': 0.1}


def run(args, **kwargs):
    result = subprocess.run([str(a) for a in args], cwd=REPO, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, **kwargs)
    if result.returncode:
        raise RuntimeError(f"{shlex.join(map(str,args))}\n{result.stdout}{result.stderr}")
    return result.stdout


def require(condition, message):
    if not condition:
        raise ValueError(message)


def input_header(cfg, model):
    # Deliberately restricted fixture adapter, not a general model JSON runtime.
    require(cfg['column']['stages'] == 5 and len(cfg['components']) == 2,
            'expected the 5x2 fixture')
    require(cfg['column']['feeds'] == [{'id': 'main_feed', 'stage': 3}], 'unexpected feeds')
    feed = model['feeds']['main_feed']
    spec = model['specifications']
    require(set(model['feeds']) == {'main_feed'}, 'unexpected runtime feeds')
    require(feed['flow']['unit'] == spec['distillate_flow']['unit'] == 'mol_s', 'expected mol_s')
    require(feed['composition']['basis'] == spec['reflux_ratio']['basis'] == 'mole', 'expected mole basis')
    require(feed['state'] == {'type': 'saturated_liquid'} and spec['type'] == 'rating', 'unsupported fixture state')
    ids = [c['id'] for c in cfg['components']]
    require(set(feed['composition']['fractions']) == set(ids), 'component IDs differ')
    values = [feed['flow']['value'], *[feed['composition']['fractions'][i] for i in ids],
              spec['distillate_flow']['value'], spec['reflux_ratio']['value']]
    require(all(math.isfinite(v) for v in values), 'nonfinite input')
    f,z0,z1,d,r = (format(v,'.17g') for v in values)
    return f'''/* Adapted from model_input.json by the validation harness. */
static const struct lumo_prepare_input input = {{
    .feeds = {{{{.flow={f}, .flow_unit=LUMO_FLOW_MOL_S,
                .composition_basis=LUMO_COMPOSITION_MOLE, .fractions={{{z0},{z1}}},
                .state=LUMO_FEED_SATURATED_LIQUID, .state_value=0}}}},
    .distillate_flow={d}, .distillate_unit=LUMO_FLOW_MOL_S, .reflux_ratio={r}
}};
'''


def compare(actual, reference, cfg, model):
    require(actual['status'] == 0 and actual['residual_norm'] <= 1e-10, 'solver did not converge')
    require(len(actual['state']) == len(actual['residual']) == 36, 'wrong MESH dimension')
    require(all(math.isfinite(v) for v in actual['state']+actual['residual']), 'nonfinite output')
    require(actual['events'] == (actual['iterations']+1 if actual['diagnostics'] else 0), 'wrong event count')
    rows = []
    maxima = dict.fromkeys(TOLERANCES, 0.0)

    def check(label, value, expected, kind):
        error = abs(value-expected)
        require(math.isfinite(value) and error <= TOLERANCES[kind],
                f'{label}: actual={value}, reference={expected}, error={error}, tolerance={TOLERANCES[kind]}')
        maxima[kind] = max(maxima[kind],error)
        rows.append({'quantity':label,'actual':value,'reference':expected,'absolute_error':error})

    state = actual['state']
    for j,stage in enumerate(reference['stages']):
        b = j*7
        check(f'T[{j+1}]',state[b+4],stage['T_K'],'T_K')
        for i in range(2):
            check(f'x[{j+1},{i}]',state[b+i],stage['x'][i],'composition')
            check(f'y[{j+1},{i}]',state[b+2+i],stage['y'][i],'composition')
        check(f'L[{j+1}]',state[b+5],stage['L_mol_s'],'flow_mol_s')
        if j:
            check(f'V[{j+1}]',state[b+6],stage['V_mol_s'],'flow_mol_s')
    check('QD',state[6],reference['condenser_heat_removed_W'],'heat_W')
    check('QW',state[35],reference['reboiler_heat_added_W'],'heat_W')
    # Independently recompute overall balances in physical units.
    feed = model['feeds']['main_feed']
    F = feed['flow']['value']
    D = model['specifications']['distillate_flow']['value']
    B = state[33]
    properties = [c['properties']['enthalpy'] for c in cfg['components']]

    def hl(offset):
        return sum(state[offset+i]*p['liquid_cp_J_mol_K']*(state[offset+4]-p['reference_temperature_K'])
                   for i,p in enumerate(properties))

    mass = abs(F-D-B)
    components = [abs(F*feed['composition']['fractions'][c['id']]-D*state[i]-B*state[28+i])
                  for i,c in enumerate(cfg['components'])]
    energy = abs(F*actual['feed_enthalpy_J_mol']+state[35]-state[6]-D*hl(0)-B*hl(28))
    # Raw residual row groups, deliberately not a mixture of units.
    residual = actual['residual']
    raw_mass = max(abs(residual[j*7+i]) for j in range(5) for i in range(2))
    raw_equilibrium = max(abs(residual[j*7+i]) for j in range(5) for i in range(2,6))
    raw_energy = max(abs(residual[j*7+6]) for j in range(5))
    checks = {'overall_total_mol_s':mass, 'overall_component_mol_s':max(components),
              'overall_energy_W':energy, 'max_M_mol_s':raw_mass,
              'max_E_S_dimensionless':raw_equilibrium, 'max_H_W':raw_energy,
              'reflux_closure_mol_s':abs(residual[-1]),
              'feed_enthalpy_error_J_mol':abs(actual['feed_enthalpy_J_mol']-reference['feed_enthalpy_J_mol'])}
    require(max(mass,*components,raw_mass,abs(residual[-1])) <= 1e-6, f'material balance: {checks}')
    require(max(energy,raw_energy) <= 0.1, f'energy balance: {checks}')
    require(raw_equilibrium <= 1e-8, f'equilibrium or normalization: {checks}')
    require(checks['feed_enthalpy_error_J_mol'] <= 1e-6, 'feed preparation differs')
    return {'iterations':actual['iterations'], 'backtracks':actual['backtracks'],
            'residual_norm':actual['residual_norm'], 'events':actual['events'],
            'max_absolute_errors':maxima, 'checks':checks, 'comparison':rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cc', default=os.environ.get('CC','cc'), help='C compiler executable')
    parser.add_argument('--output', type=Path, default=HERE/'validation', help='report directory')
    args = parser.parse_args()
    cfg = json.loads((HERE/'compiler_input.json').read_text())
    model = json.loads((HERE/'model_input.json').read_text())
    reference = json.loads((HERE/'reference.json').read_text())
    for name,digest in reference['input_sha256'].items():
        require(hashlib.sha256((HERE/name).read_bytes()).hexdigest() == digest, f'stale reference input: {name}')
    flags = ['-std=c99','-O2','-Wall','-Wextra','-Werror','-pedantic']
    report = {'input_sha256':reference['input_sha256'],
              'reference_sha256':hashlib.sha256((HERE/'reference.json').read_bytes()).hexdigest(),
              'revision':run(['git','rev-parse','HEAD']).strip(),
              'working_tree':run(['git','status','--porcelain']).splitlines(),
              'go_version':run(['go','version']).strip(),
              'cc_version':run([args.cc,'--version']).splitlines()[0],
              'c_flags':flags, 'tolerances':TOLERANCES, 'variants':{}}
    outputs = {}
    with tempfile.TemporaryDirectory(prefix='lumo-validation-') as temp:
        temp = Path(temp)
        lumo = temp/'lumo'
        run(['go','build','-o',lumo,'./cmd/lumo'])
        for diagnostics in (False,True):
            name = 'diagnostics' if diagnostics else 'off'
            generated = temp/name
            command = [lumo,'-o',generated]
            if diagnostics:
                command.append('-diagnostics')
            run([*command,HERE/'compiler_input.json'])
            (generated/'input.h').write_text(input_header(cfg,model))
            binary = temp/f'run-{name}'
            run([args.cc,*flags,'-I',generated,HERE/'driver.c',*sorted(generated.glob('*.c')),'-lm','-o',binary])
            actual = json.loads(run([binary]))
            require(actual['diagnostics'] == int(diagnostics),'wrong generated diagnostics mode')
            report['variants'][name] = compare(actual,reference,cfg,model)
            report['variants'][name]['generated_sha256'] = {
                path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                for path in sorted(generated.iterdir())
            }
            outputs[name] = actual
    for key in ('state','residual','status','iterations','backtracks','residual_norm'):
        require(outputs['off'][key] == outputs['diagnostics'][key],f'diagnostics changes {key}')
    report['diagnostics_results_equal'] = True
    args.output.mkdir(parents=True,exist_ok=True)
    (args.output/'report.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    lines = ['# Сравнение Lumo с эталоном 5×2', '',
             'Обе сборки сошлись, состояния и невязки совпали точно в этом прогоне', '',
             '| Сборка | Итерации | Отклонённые пробы | Норма невязки | События |',
             '|---|---|---|---|---|']
    for name,result in report['variants'].items():
        lines.append(f"| {name} | {result['iterations']} | {result['backtracks']} | {result['residual_norm']:.6g} | {result['events']} |")
    lines += ['', '| Величина | Максимальное отклонение OFF | Допуск |','|---|---|---|']
    for kind,error in report['variants']['off']['max_absolute_errors'].items():
        lines.append(f'| {kind} | {error:.6g} | {TOLERANCES[kind]:.6g} |')
    lines += ['', 'Полное сравнение, проверки балансов и сведения о сборке находятся в `report.json`',
              'Это проверка численной реализации идеальной модели, точность реальных свойств смеси здесь не оценивается','']
    (args.output/'report.md').write_text('\n'.join(lines))
    print('\n'.join(lines))


if __name__ == '__main__':
    main()
