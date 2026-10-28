#!/usr/bin/env python3
"""Expand and build the periodic driver without Flash-X modules or headers."""
import argparse
from pathlib import Path
import shutil
import subprocess
import macro_expand as macros

ROOT = Path(__file__).resolve().parent

def build(output, layout="variable-first", compiler="gfortran"):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    definitions = [ROOT / name for name in [
        "hydro_layout.ini", "hydro_helpers.ini", "update_helpers.ini",
        "diagnostic_helpers.ini", "driver_helpers.ini"]]
    overrides = {"spatial-first": "hydro_layout_spatial_first.ini", "records": "hydro_layout_records.ini"}
    if layout in overrides:
        definitions.append(ROOT / overrides[layout])
    expander = macros.MacroExpander(macros.load_macros(definitions))
    sources = []
    for name in [
        "hy_getFaceFlux", "hy_updateSolution", "hy_shockDetect",
        "hy_computeDt", "hydro_grid_contract", "hy_prepareAdvance",
        "periodic_grid_example", "periodic_shock_tube"]:
        template = ROOT / (name + ".F90-mc")
        source = output / (name + ".f90")
        if template.exists():
            macros.expand_file(template, source, expander)
        else:
            source.write_text(macros.wrap_fortran((ROOT / (name + ".F90")).read_text()))
        sources.append(source.name)
    subprocess.run([compiler, "-std=f2003", "-O2", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
                    *sources, "-o", "periodic_shock_tube"], cwd=output, check=True)
    return output / "periodic_shock_tube"

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "build" / "standalone")
    parser.add_argument("--layout", choices=["variable-first", "spatial-first", "records"],
                        default="variable-first")
    parser.add_argument("--compiler", default="gfortran")
    args = parser.parse_args()
    if not shutil.which(args.compiler):
        parser.error("Fortran compiler is unavailable")
    print(build(args.output, args.layout, args.compiler))
