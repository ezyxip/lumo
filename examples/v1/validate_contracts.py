"""Validate v1 schemas and contract examples.

Run: uv run --with jsonschema python examples/v1/validate_contracts.py
Synthetic values below exercise schema branches; they are not property data
or calculated results.
"""

import copy
import json
from pathlib import Path

from jsonschema import Draft202012Validator


ROOT = Path(__file__).resolve().parents[2]


def load(path):
    return json.loads((ROOT / path).read_text())


def main():
    validators = {}
    for name in ("compiler_input", "model_input", "model_output"):
        schema = load(f"examples/schema/v1/{name}.schema.json")
        Draft202012Validator.check_schema(schema)
        validators[name] = Draft202012Validator(schema)

    count = 0

    def check(name, data, valid=True):
        nonlocal count
        errors = list(validators[name].iter_errors(data))
        if valid and errors:
            raise AssertionError(f"{name}: {errors[0].message}")
        if not valid and not errors:
            raise AssertionError(f"{name}: invalid case was accepted: {data}")
        count += 1

    compiler = load("examples/v1/benzene_toluene/compiler_input.json")
    model = load("examples/v1/benzene_toluene/model_input.json")
    check("compiler_input", compiler)
    check("model_input", model)

    # Multiple components/feeds have the same shape as the binary example.
    expanded = copy.deepcopy(compiler)
    expanded["components"].append({**copy.deepcopy(compiler["components"][0]), "id": "component_c"})
    expanded["column"]["feeds"].append({"id": "second_feed", "stage": 6})
    check("compiler_input", expanded)
    expanded_model = copy.deepcopy(model)
    feed = expanded_model["feeds"]["main_feed"]
    feed["composition"]["fractions"] = {"benzene": 0.4, "toluene": 0.4, "component_c": 0.2}
    expanded_model["feeds"]["second_feed"] = copy.deepcopy(feed)
    check("model_input", expanded_model)

    # Synthetic inline property data, only for structural checks.
    inline = copy.deepcopy(compiler)
    inline["components"][0]["properties"] = {
        "source": "inline",
        "antoine": {"A": 1, "B": 1, "C": 1, "pressure_unit": "mmHg", "temperature_unit": "C"},
        "enthalpy": {
            "type": "constant_cp", "reference_temperature_K": 300,
            "liquid_cp_J_mol_K": 100, "vapor_cp_J_mol_K": 50,
            "vaporization_enthalpy_at_reference_J_mol": 30000,
        },
        "temperature_range_K": {"min": 290, "max": 400},
    }
    check("compiler_input", inline)
    inline["components"][0]["properties"]["key"] = "mixed_branch"
    check("compiler_input", inline, valid=False)

    for state in (
        {"type": "temperature", "temperature_K": 350},
        {"type": "vapor_fraction", "vapor_fraction": 0},
        {"type": "vapor_fraction", "vapor_fraction": 1},
    ):
        candidate = copy.deepcopy(model)
        candidate["feeds"]["main_feed"]["state"] = state
        check("model_input", candidate)

    for state in (
        {},
        {"type": "temperature"},
        {"type": "saturated_liquid", "temperature_K": 350},
        {"type": "vapor_fraction", "vapor_fraction": 1.1},
        {"type": "temperature", "temperature_K": 350, "vapor_fraction": 0.5},
    ):
        candidate = copy.deepcopy(model)
        candidate["feeds"]["main_feed"]["state"] = state
        check("model_input", candidate, valid=False)

    for path, value in (
        (("schema_version",), 2),
        (("specifications",), {"type": "rating"}),
        (("specifications", "reflux_ratio", "basis"), "mass"),
        (("feeds", "main_feed", "composition", "fractions", "benzene"), -0.1),
        (("feeds", "main_feed", "flow", "unit"), "unknown"),
    ):
        candidate = copy.deepcopy(model)
        target = candidate
        for key in path[:-1]:
            target = target[key]
        target[path[-1]] = value
        check("model_input", candidate, valid=False)
    for field in ("schema_version", "name", "components", "thermo", "column"):
        candidate = copy.deepcopy(compiler)
        del candidate[field]
        check("compiler_input", candidate, valid=False)
    candidate = copy.deepcopy(compiler)
    candidate["column"]["pressure"]["type"] = "linear_drop"
    check("compiler_input", candidate, valid=False)
    candidate = copy.deepcopy(compiler)
    candidate["column"]["feeds"] = None
    check("compiler_input", candidate, valid=False)

    # A failed solve may omit results; a successful one must supply them.
    output = {
        "schema_version": 1,
        "provenance": {"solver": "schema_test", "compiler_input": "test", "model_input": "test"},
        "diagnostics": {"converged": False, "iterations": 0, "residual_norm": 1},
    }
    check("model_output", output)
    output["diagnostics"]["converged"] = True
    check("model_output", output, valid=False)
    stream = {"flow_mol_s": 1, "mole_fractions": {"benzene": 0.5, "toluene": 0.5}}
    output["products"] = {"distillate": stream, "bottoms": stream, "condenser_duty_W": 1, "reboiler_duty_W": 1}
    output["profiles"] = [
        {"stage": j, "temperature_K": 350, "liquid": stream, "vapor": None if j == 1 else stream}
        for j in range(1, 4)
    ]
    check("model_output", output, valid=False)  # missing balance diagnostic
    output["diagnostics"]["mass_balance_error"] = 0
    check("model_output", output)
    output["diagnostics"]["mass_balance_error"] = None
    check("model_output", output, valid=False)

    # Cross-file semantic checks for the published example pair.
    ids = {component["id"] for component in compiler["components"]}
    feed_ids = {feed["id"] for feed in compiler["column"]["feeds"]}
    assert len(ids) == len(compiler["components"])
    assert len(feed_ids) == len(compiler["column"]["feeds"])
    assert feed_ids == set(model["feeds"])
    for feed in compiler["column"]["feeds"]:
        assert 2 <= feed["stage"] < compiler["column"]["stages"]
    for feed in model["feeds"].values():
        fractions = feed["composition"]["fractions"]
        assert set(fractions) == ids
        assert abs(sum(fractions.values()) - 1) <= 1e-10
    print(f"Validated 3 schemas, {count} contract cases, and example cross-file references.")


if __name__ == "__main__":
    main()
