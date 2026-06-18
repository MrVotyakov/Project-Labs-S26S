# Automatic Cassandra Batch Runs

Put the experimental table in this directory as `experimental_points.csv`, or
pass another file with `--table`.

Required liquid-composition columns:

```csv
experiment,x_water,x_ethanol,x_hexane
0001,0.105,0.236667,0.658333
```

Accepted aliases include `x1/x2/x3`, `water/ethanol/hexane`, and Russian
headers such as `x_вода`, `x_спирт`, `x_гексан`.

Run:

```bash
./run_auto_batch.sh --table experimental_points.csv
```

Generate inputs without running Cassandra:

```bash
./run_auto_batch.sh --table experimental_points.csv --prepare-only
```

Run only selected experiments:

```bash
./run_auto_batch.sh --table experimental_points.csv --experiments 0001,0003
```

Fast smoke test on a smaller system:

```bash
./run_auto_batch.sh --table exp21withGamma.csv --output runs_exp21_smoke --experiments 0001,0002 --smoke-test
```

The smoke test uses 600 liquid molecules, 150 vapor molecules, 1000 MC steps,
and `prop_freq 10`, so property output should appear quickly. It is meant to
test the automation pipeline, not to produce final VLE statistics.

Edit `auto_cassandra_batch.py` to change:

- `LIQUID_TOTAL_MOLECULES`
- `VAPOR_TOTAL_MOLECULES`
- `STAGES`
- `COMPOSITION_MODE`
- box sizes and move settings

Each experiment is written to:

```text
runs/exp_<id>/stage_XX_<name>/
```

Each experiment folder also receives `metadata.csv` with mole fractions and
the integer molecule counts used for liquid, vapor, and total system.
