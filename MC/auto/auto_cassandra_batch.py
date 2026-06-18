#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path


# ---------------------------------------------------------------------------
# User settings
# ---------------------------------------------------------------------------

TABLE_FILE = "experimental_points.csv"
OUTPUT_DIR = "runs"

# One row creates one experiment.  The liquid box starts from x, and the vapor
# box starts from the same x but with fewer particles.
LIQUID_TOTAL_MOLECULES = 2400
VAPOR_TOTAL_MOLECULES = 600

# Initial box lengths in Angstrom.  GEMC-NPT volume moves will then adjust them.
LIQUID_BOX_LENGTH = 71.4
VAPOR_BOX_LENGTH = 300.0

TEMPERATURE_K = 329.55
PRESSURE_BAR = 1.01325

# If the table has mole fractions, use "mole_fraction".
# Other supported modes: "mole_percent", "mass_percent".
COMPOSITION_MODE = "mole_fraction"

# Edit this list to define n stages.  Stage 1 starts from make_config; later
# stages start from the previous stage checkpoint.
STAGES = [
    {
        "name": "eq_10k",
        "steps": 10_000,
        "run_type": "equilibration 1000 100",
        "prop_freq": 100,
        "coord_freq": 1000,
        "seed": (31415, 92653),
        "widom": False,
    },
]

SMOKE_TEST_SETTINGS = {
    "liquid_total": 600,
    "vapor_total": 150,
    "liquid_box_length": 45.0,
    "vapor_box_length": 190.0,
    "stages": [
        {
            "name": "smoke_1k",
            "steps": 1_000,
            "run_type": "equilibration 100 20",
            "prop_freq": 10,
            "coord_freq": 100,
            "seed": (27182, 81828),
            "widom": False,
        },
    ],
}

# Move settings copied from the azeotrope GEMC model.
MOVE_PROB_TRANSLATION = 0.35
MOVE_TRANSLATION_WIDTHS = ((1.0, 1.0, 1.0), (20.0, 20.0, 20.0))

MOVE_PROB_ROTATION = 0.25
MOVE_ROTATION_WIDTHS = ((30.0, 30.0, 30.0), (180.0, 180.0, 180.0))

MOVE_PROB_REGROWTH = 0.15
MOVE_REGROWTH_SPECIES = (0.0, 0.3, 0.7)

MOVE_PROB_SWAP = 0.20
MOVE_PROB_VOLUME = 0.05
MOVE_VOLUME_WIDTHS = (5000.0, 500000.0)

VDW_CUTOFFS = (12.0, 50.0)
CHARGE_CUTOFFS = (12.0, 50.0)
EWALD_ACCURACY = "1E-5"
RCUTOFF_LOW = 0.8
KAPPA_INS = 12
KAPPA_DIH = 10

# Widom settings are only used when a stage has "widom": True.
WIDOM_INSERTIONS = [
    ("cbmc 5000 5000 50", "cbmc 2000 5000 20"),
    ("cbmc 10000 5000 50", "cbmc 3000 5000 20"),
    ("cbmc 15000 5000 50", "cbmc 5000 5000 20"),
]


# ---------------------------------------------------------------------------
# Internal definitions
# ---------------------------------------------------------------------------

SPECIES = ("water", "ethanol", "hexane")
MOLAR_MASSES = {
    "water": 18.01528,
    "ethanol": 46.06844,
    "hexane": 86.17536,
}

X_COLUMN_ALIASES = {
    "water": (
        "x_water",
        "x_h2o",
        "x1",
        "water",
        "h2o",
        "x_вода",
        "x_вода,%",
        "вода",
        "вода,%",
    ),
    "ethanol": (
        "x_ethanol",
        "x_etoh",
        "x2",
        "ethanol",
        "etoh",
        "x_спирт",
        "x_спирт,%",
        "спирт",
        "спирт,%",
        "этанол",
        "x_этанол",
    ),
    "hexane": (
        "x_hexane",
        "x_nhexane",
        "x3",
        "hexane",
        "nhexane",
        "x_гексан",
        "x_гексан,%",
        "гексан",
        "гексан,%",
    ),
}

ID_COLUMN_ALIASES = ("experiment", "experiment_id", "exp", "id", "n", "номер", "точка")


@dataclass(frozen=True)
class ExperimentPoint:
    exp_id: str
    fractions: tuple[float, float, float]
    source_row: dict[str, str]


def project_root() -> Path:
    return Path(__file__).resolve().parents[2]


def auto_dir() -> Path:
    return Path(__file__).resolve().parent


def mc_data_dir() -> Path:
    return project_root() / "MC" / "data"


def cassandra_exe() -> Path:
    return project_root() / "Cassandra" / "bin" / "cassandra.exe"


def cassandra_lib_dir() -> Path:
    return project_root() / "Cassandra" / "Libraries" / "locals" / "gfortran" / "local" / "lib"


def normalize_column(name: str) -> str:
    return name.strip().lower().replace(" ", "_").replace("-", "_")


def parse_float(value: str) -> float:
    return float(value.strip().replace(",", "."))


def find_column(fieldnames: list[str], aliases: tuple[str, ...]) -> str | None:
    normalized = {normalize_column(name): name for name in fieldnames}
    for alias in aliases:
        key = normalize_column(alias)
        if key in normalized:
            return normalized[key]
    return None


def sniff_dialect(path: Path) -> csv.Dialect:
    sample = path.read_text(encoding="utf-8-sig")[:4096]
    try:
        return csv.Sniffer().sniff(sample, delimiters=",;\t")
    except csv.Error:
        return csv.excel


def read_table(path: Path) -> list[ExperimentPoint]:
    if not path.exists():
        raise FileNotFoundError(f"Experimental table not found: {path}")

    dialect = sniff_dialect(path)
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, dialect=dialect)
        if not reader.fieldnames:
            raise ValueError(f"Table has no header row: {path}")
        fieldnames = list(reader.fieldnames)

        composition_columns = []
        for species in SPECIES:
            column = find_column(fieldnames, X_COLUMN_ALIASES[species])
            if not column:
                raise ValueError(
                    f"Cannot find liquid composition column for {species}. "
                    f"Known columns: {', '.join(fieldnames)}"
                )
            composition_columns.append(column)

        id_column = find_column(fieldnames, ID_COLUMN_ALIASES)

        points = []
        exp_counter = 0
        for row in reader:
            if not any((value or "").strip() for value in row.values()):
                continue
            exp_counter += 1
            raw_values = [parse_float(row[column]) for column in composition_columns]
            fractions = composition_to_mole_fractions(raw_values)
            exp_id = row[id_column].strip() if id_column and row.get(id_column) else f"{exp_counter:04d}"
            points.append(ExperimentPoint(sanitize_id(exp_id), fractions, dict(row)))

    return points


def sanitize_id(value: str) -> str:
    clean = "".join(char if char.isalnum() or char in ("_", "-") else "_" for char in value)
    return clean or "unnamed"


def composition_to_mole_fractions(values: list[float]) -> tuple[float, float, float]:
    if any(value < 0 for value in values):
        raise ValueError(f"Composition cannot contain negative values: {values}")

    if COMPOSITION_MODE == "mole_fraction":
        total = sum(values)
        if total > 1.5:
            values = [value / 100.0 for value in values]
    elif COMPOSITION_MODE == "mole_percent":
        values = [value / 100.0 for value in values]
    elif COMPOSITION_MODE == "mass_percent":
        mass_fractions = values
        if sum(mass_fractions) <= 1.5:
            mass_fractions = [value * 100.0 for value in mass_fractions]
        moles = [
            mass_fractions[0] / MOLAR_MASSES["water"],
            mass_fractions[1] / MOLAR_MASSES["ethanol"],
            mass_fractions[2] / MOLAR_MASSES["hexane"],
        ]
        values = moles
    else:
        raise ValueError(f"Unknown COMPOSITION_MODE: {COMPOSITION_MODE}")

    total = sum(values)
    if total <= 0:
        raise ValueError(f"Composition sums to zero: {values}")

    return tuple(value / total for value in values)  # type: ignore[return-value]


def integer_counts(total: int, fractions: tuple[float, float, float]) -> tuple[int, int, int]:
    exact = [total * fraction for fraction in fractions]
    counts = [math.floor(value) for value in exact]
    remainder = total - sum(counts)
    order = sorted(range(len(exact)), key=lambda i: exact[i] - counts[i], reverse=True)
    for index in order[:remainder]:
        counts[index] += 1

    for index, fraction in enumerate(fractions):
        if fraction > 0 and counts[index] == 0 and total >= len(fractions):
            donor = max(range(len(counts)), key=lambda i: counts[i])
            counts[donor] -= 1
            counts[index] = 1

    if sum(counts) != total:
        raise AssertionError("Integer count rounding failed")
    return tuple(counts)  # type: ignore[return-value]


def rel(path: Path, start: Path) -> str:
    return os.path.relpath(path, start)


def molecule_file(species: str, stage_dir: Path) -> str:
    path = auto_dir() / f"{species}.mcf"
    if not path.exists():
        path = mc_data_dir() / f"{species}.mcf"
    return rel(path, stage_dir)


def fragment_lines(stage_dir: Path) -> list[str]:
    data = mc_data_dir()
    fragments = [
        ("species1/frag1/frag1.dat", 1),
        ("species2/frag1/frag1.dat", 2),
        ("species2/frag2/frag2.dat", 3),
        ("species2/frag3/frag3.dat", 4),
        ("species3/frag1/frag1.dat", 5),
        ("species3/frag2/frag2.dat", 6),
        ("species3/frag3/frag3.dat", 7),
        ("species3/frag4/frag4.dat", 8),
        ("species3/frag5/frag5.dat", 9),
        ("species3/frag6/frag6.dat", 10),
    ]
    return [f"{rel(data / path, stage_dir)}  {index}" for path, index in fragments]


def counts_line(counts: tuple[int, int, int]) -> str:
    return f"{counts[0]} {counts[1]} {counts[2]}"


def triple_line(values: tuple[float, float, float]) -> str:
    return f"{values[0]} {values[1]} {values[2]}"


def widom_section(stage: dict) -> str:
    if not stage.get("widom", False):
        return ""
    lines = ["", "# Widom_Insertion", "true"]
    for box1_entry, box2_entry in WIDOM_INSERTIONS:
        lines.append(f"{box1_entry} {box2_entry}")
    lines.extend(["", "# Cell_List_Overlap", "true"])
    return "\n".join(lines)


def build_input(
    stage_dir: Path,
    run_name: str,
    stage: dict,
    liquid_counts: tuple[int, int, int],
    vapor_counts: tuple[int, int, int],
    previous_checkpoint: Path | None,
) -> str:
    total_counts = tuple(liquid_counts[i] + vapor_counts[i] for i in range(3))
    seed = stage.get("seed", (12345, 67890))
    start_type = (
        f"checkpoint {rel(previous_checkpoint, stage_dir)}"
        if previous_checkpoint
        else f"make_config {counts_line(liquid_counts)}\nmake_config {counts_line(vapor_counts)}"
    )

    return f"""# Run_Name
{run_name}

# Sim_Type
gemc_npt

# Nbr_Species
3

# Molecule_Files
{molecule_file("water", stage_dir)} {total_counts[0]}
{molecule_file("ethanol", stage_dir)} {total_counts[1]}
{molecule_file("hexane", stage_dir)} {total_counts[2]}

# Fragment_Files
{chr(10).join(fragment_lines(stage_dir))}

# Box_Info
2
cubic
{LIQUID_BOX_LENGTH}

cubic
{VAPOR_BOX_LENGTH}

# VDW_Style
lj cut_tail {VDW_CUTOFFS[0]}
lj cut_tail {VDW_CUTOFFS[1]}

# Charge_Style
coul ewald {CHARGE_CUTOFFS[0]} {EWALD_ACCURACY}
coul ewald {CHARGE_CUTOFFS[1]} {EWALD_ACCURACY}

# Mixing_Rule
lb

# Temperature_Info
{TEMPERATURE_K}
{TEMPERATURE_K}

# Pressure_Info
{PRESSURE_BAR}
{PRESSURE_BAR}

# Move_Probability_Info

# Prob_Translation
{MOVE_PROB_TRANSLATION}
{triple_line(MOVE_TRANSLATION_WIDTHS[0])}
{triple_line(MOVE_TRANSLATION_WIDTHS[1])}

# Prob_Rotation
{MOVE_PROB_ROTATION}
{triple_line(MOVE_ROTATION_WIDTHS[0])}
{triple_line(MOVE_ROTATION_WIDTHS[1])}

# Prob_Regrowth
{MOVE_PROB_REGROWTH}
{triple_line(MOVE_REGROWTH_SPECIES)}

# Prob_Swap
{MOVE_PROB_SWAP}
cbmc cbmc cbmc

# Prob_Volume
{MOVE_PROB_VOLUME}
{MOVE_VOLUME_WIDTHS[0]}
{MOVE_VOLUME_WIDTHS[1]}

# Done_Probability_Info

# Start_Type
{start_type}

# Run_Type
{stage["run_type"]}

# Simulation_Length_Info
units steps
prop_freq {stage["prop_freq"]}
coord_freq {stage["coord_freq"]}
run {stage["steps"]}

# Property_Info 1
energy_total
pressure
volume
mass_density
density
nmols

# Property_Info 2
energy_total
pressure
volume
mass_density
density
nmols

# CBMC_Info
kappa_ins {KAPPA_INS}
kappa_dih {KAPPA_DIH}
rcut_cbmc {VDW_CUTOFFS[0]} {VDW_CUTOFFS[1]}

# Rcutoff_Low
{RCUTOFF_LOW}

# Seed_Info
{seed[0]} {seed[1]}

# Average_Info
1

# Pair_Energy
true
{widom_section(stage)}

END
"""


def read_last_property_row(path: Path) -> tuple[int, list[float], list[float]] | None:
    if not path.exists():
        return None
    for line in reversed(path.read_text(errors="ignore").splitlines()):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split()
        if len(parts) < 11:
            continue
        try:
            step = int(float(parts[0]))
            nmols = [float(value) for value in parts[8:11]]
        except ValueError:
            continue
        total = sum(nmols)
        if total <= 0:
            return None
        return step, nmols, [value / total for value in nmols]
    return None


def read_last_log_line(path: Path) -> str:
    if not path.exists():
        return "log is not created yet"
    for line in reversed(path.read_text(errors="ignore").splitlines()):
        stripped = line.strip()
        if stripped:
            return stripped[:180]
    return "log is empty"


def log_has_fatal_error(path: Path) -> bool:
    if not path.exists():
        return False
    text = path.read_text(errors="ignore")
    return "Fatal Error" in text or "************ERROR************" in text


def fmt(values: list[float], digits: int = 4) -> str:
    return "/".join(f"{value:.{digits}f}" for value in values)


def monitor_process(process: subprocess.Popen, run_prefix: str, stage_dir: Path) -> int:
    box1_path = stage_dir / f"{run_prefix}.box1.prp"
    box2_path = stage_dir / f"{run_prefix}.box2.prp"
    log_path = stage_dir / f"{run_prefix}.log"
    last_step = None
    last_wait = 0.0

    while process.poll() is None:
        box1 = read_last_property_row(box1_path)
        box2 = read_last_property_row(box2_path)
        if box1 and box2:
            step = min(box1[0], box2[0])
            if step != last_step:
                last_step = step
                k_values = [
                    box2[2][i] / box1[2][i] if box1[2][i] > 0 else float("nan")
                    for i in range(3)
                ]
                print(
                    f"  MC_STEP {step:>8} | "
                    f"x_liq={fmt(box1[2])} y_vap={fmt(box2[2])} K={fmt(k_values)}",
                    flush=True,
                )
        else:
            now = time.monotonic()
            if now - last_wait >= 5:
                last_wait = now
                print(
                    f"  waiting for property output... log: {read_last_log_line(log_path)}",
                    flush=True,
                )
        time.sleep(2)

    status = process.wait()
    if log_has_fatal_error(log_path):
        print(f"  Cassandra reported a fatal error in {log_path.name}", flush=True)
        return status if status != 0 else 1
    return status


def write_metadata(path: Path, point: ExperimentPoint, liquid_counts, vapor_counts) -> None:
    rows = [
        "experiment_id,component,mole_fraction,liquid_count,vapor_count,total_count",
    ]
    for i, species in enumerate(SPECIES):
        rows.append(
            f"{point.exp_id},{species},{point.fractions[i]:.10f},"
            f"{liquid_counts[i]},{vapor_counts[i]},{liquid_counts[i] + vapor_counts[i]}"
        )
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")


def prepare_experiment(point: ExperimentPoint, output_root: Path) -> list[tuple[Path, str]]:
    liquid_counts = integer_counts(LIQUID_TOTAL_MOLECULES, point.fractions)
    vapor_counts = integer_counts(VAPOR_TOTAL_MOLECULES, point.fractions)

    exp_dir = output_root / f"exp_{point.exp_id}"
    exp_dir.mkdir(parents=True, exist_ok=True)
    write_metadata(exp_dir / "metadata.csv", point, liquid_counts, vapor_counts)

    generated = []
    previous_checkpoint = None
    for index, stage in enumerate(STAGES, start=1):
        stage_name = sanitize_id(stage["name"])
        stage_dir = exp_dir / f"stage_{index:02d}_{stage_name}"
        stage_dir.mkdir(parents=True, exist_ok=True)
        run_name = f"exp_{point.exp_id}_s{index:02d}_{stage_name}"
        input_path = stage_dir / f"{run_name}.inp"
        input_path.write_text(
            build_input(
                stage_dir=stage_dir,
                run_name=run_name,
                stage=stage,
                liquid_counts=liquid_counts,
                vapor_counts=vapor_counts,
                previous_checkpoint=previous_checkpoint,
            ),
            encoding="utf-8",
        )
        generated.append((input_path, run_name))
        previous_checkpoint = stage_dir / f"{run_name}.chk"

    return generated


def run_stage(input_path: Path, run_name: str) -> int:
    env = os.environ.copy()
    lib_dir = str(cassandra_lib_dir())
    env["DYLD_LIBRARY_PATH"] = (
        lib_dir if not env.get("DYLD_LIBRARY_PATH") else lib_dir + ":" + env["DYLD_LIBRARY_PATH"]
    )

    print(f"Running {input_path}", flush=True)
    process = subprocess.Popen(
        [str(cassandra_exe()), input_path.name],
        cwd=str(input_path.parent),
        env=env,
        stdout=None,
        stderr=None,
    )

    return monitor_process(process, run_name, input_path.parent)


def parse_experiment_filter(value: str | None) -> set[str] | None:
    if not value:
        return None
    return {sanitize_id(item.strip()) for item in value.split(",") if item.strip()}


def main() -> int:
    global LIQUID_TOTAL_MOLECULES
    global VAPOR_TOTAL_MOLECULES
    global LIQUID_BOX_LENGTH
    global VAPOR_BOX_LENGTH
    global STAGES

    parser = argparse.ArgumentParser(
        description="Generate and optionally run Cassandra GEMC jobs for experimental points."
    )
    parser.add_argument("--table", default=TABLE_FILE, help="CSV/TSV table with x composition columns.")
    parser.add_argument("--output", default=OUTPUT_DIR, help="Output directory inside MC/auto.")
    parser.add_argument("--prepare-only", action="store_true", help="Generate folders and .inp files only.")
    parser.add_argument("--experiments", help="Comma-separated experiment ids to prepare/run.")
    parser.add_argument(
        "--smoke-test",
        action="store_true",
        help="Use a small fast system to verify the pipeline before the full 3000-molecule run.",
    )
    parser.add_argument(
        "--steps",
        type=int,
        help="Override the number of MC steps in every configured stage.",
    )
    args = parser.parse_args()

    if args.smoke_test:
        LIQUID_TOTAL_MOLECULES = SMOKE_TEST_SETTINGS["liquid_total"]
        VAPOR_TOTAL_MOLECULES = SMOKE_TEST_SETTINGS["vapor_total"]
        LIQUID_BOX_LENGTH = SMOKE_TEST_SETTINGS["liquid_box_length"]
        VAPOR_BOX_LENGTH = SMOKE_TEST_SETTINGS["vapor_box_length"]
        STAGES = list(SMOKE_TEST_SETTINGS["stages"])

    if args.steps is not None:
        if args.steps <= 0:
            print("--steps must be positive.", file=sys.stderr)
            return 1
        STAGES = [dict(stage, steps=args.steps) for stage in STAGES]

    table_path = Path(args.table)
    if not table_path.is_absolute():
        table_path = auto_dir() / table_path

    output_root = Path(args.output)
    if not output_root.is_absolute():
        output_root = auto_dir() / output_root
    output_root.mkdir(parents=True, exist_ok=True)

    selected = parse_experiment_filter(args.experiments)
    points = read_table(table_path)
    if selected:
        points = [point for point in points if point.exp_id in selected]
    if not points:
        print("No experiment points selected.", file=sys.stderr)
        return 1

    if not cassandra_exe().exists():
        print(f"Cassandra executable not found: {cassandra_exe()}", file=sys.stderr)
        return 1

    print(f"Loaded {len(points)} experiment point(s) from {table_path}")
    print(f"Output root: {output_root}")
    print(f"Stages: {', '.join(stage['name'] for stage in STAGES)}")

    all_inputs = []
    for point in points:
        generated = prepare_experiment(point, output_root)
        all_inputs.extend(generated)
        print(
            f"Prepared exp_{point.exp_id}: "
            f"x={fmt(list(point.fractions), 6)} "
            f"stages={len(generated)}"
        )

    if args.prepare_only:
        print("Preparation complete. Cassandra was not run.")
        return 0

    for input_path, run_name in all_inputs:
        status = run_stage(input_path, run_name)
        if status != 0:
            print(f"Stage failed with status {status}: {input_path}", file=sys.stderr)
            return status

    print("All Cassandra stages finished.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
