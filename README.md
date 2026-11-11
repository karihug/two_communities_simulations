# Two-community SIR simulations

Generate contact networks, simulate epidemics with travel between two communities, and plot baseline, travel-ban, and social-distancing trajectories.

## Files

| Task | Required files |
| --- | --- |
| Graph generation | `graphgen.cpp` |
| Epidemic simulation | `main.cpp`, `SIRSimulation.cpp`, `SIRSimulation.h` |
| Visualization | `plot_sims.py` |

Keep these five files in the same directory and run the commands below from that directory. No notebooks are required. Graphs, simulation CSVs, figures, and compiled executables are generated separately.

## Compilation

Use a C++17 compiler. Compile on the machine where the simulations will run.

```sh
g++ -O3 -std=c++17 graphgen.cpp -o graphgen
g++ -O3 -std=c++17 main.cpp SIRSimulation.cpp -o simulation_final
```

## Graph generation

The generator takes the **total population size**: `20000000` creates two communities of 10 million individuals each. The parameters `mu` and `sigma` describe the underlying normal distribution of a lognormal degree distribution.

### Lognormal graphs

Each individual receives independent home-community and travel-community degree draws from the same distribution. The draws are floored to integers. Home edges connect individuals within a community; travel edges connect individuals in opposite communities. Stub balancing and removal of self-loops and duplicate edges can reduce the realized degrees, so individual home and travel degrees are not required to be equal.

Hill graph:

```sh
mkdir -p final_submission_sims/hill_sims
./graphgen ln 20000000 1.4 1.27 \
  > final_submission_sims/hill_sims/graph_ln_20m.txt
```

Copenhagen graph, using the fitted lognormal parameters:

```sh
mkdir -p final_submission_sims/copenhagen_sims
./graphgen ln 20000000 3.1814878275740055 0.9337859022096622 \
  > final_submission_sims/copenhagen_sims/copenhagen_graph_20m.txt
```

The `ln` mode uses a time-based random seed. Keep the generated graph to reuse the same network across simulation runs.

### Degree-assortative and degree-disassortative graphs

These modes draw one total degree per individual and use an 80% home / 20% travel preference for each stub. They use the degree classes and mixing matrices encoded from the supplementary information, with a balance correction for the sampled class stub totals. Unmatched stubs, self-loops, and duplicate edges can reduce realized degrees. Generation diagnostics are printed to the terminal.

The final argument is the graph-generation seed, separate from the epidemic seed.

```sh
mkdir -p final_submission_sims/assortative_sims
./graphgen assortative 20000000 1.4 1.27 20261111 \
  > final_submission_sims/assortative_sims/assortative_graph_20m.txt

mkdir -p final_submission_sims/disassortative_sims
./graphgen disassortative 20000000 1.4 1.27 20261111 \
  > final_submission_sims/disassortative_sims/disassortative_graph_20m.txt
```

All graph files contain the total vertex count on the first line, followed by one undirected edge `u v` per line. Vertex IDs are zero-based; the first half belongs to community 1 and the second half to community 2.

## Epidemic simulation

The simulator takes the **population per community**: `-n 10000000` requires a graph containing 20 million vertices. Each run starts with one infected individual in community 1.

The following examples use the Copenhagen graph. Change `GRAPH` and `OUT` to run another network.

```sh
GRAPH=final_submission_sims/copenhagen_sims/copenhagen_graph_20m.txt
OUT=final_submission_sims/copenhagen_sims

./simulation_final --graph "$GRAPH" --output-dir "$OUT" \
  -n 10000000 -b 0.1 -g 0.2041 -T 0.0003 -H 1 -c 5 \
  --max-time 100
```

| Argument | Meaning |
| --- | --- |
| `--graph` | Input graph file |
| `--output-dir` | Directory for simulation CSVs |
| `-n` | Number of individuals in each community |
| `-b` | Transmission-attempt rate, beta |
| `-g` | Gamma recovery-time **scale**, with fixed shape 3.5 |
| `-T` | Rate of departing the home community |
| `-H` | Rate of returning home |
| `-c` | Nominal degree used in printed theoretical diagnostics; does not change the loaded graph |
| `--max-time` | Simulation time limit |
| `--seed` | Epidemic random seed; generated and printed when omitted |
| `--travel-ban` | Enable a travel ban |
| `--social-distancing 5` | Enable social distancing with strength 5 |

Recovery times follow a Gamma distribution with shape 3.5 and scale `-g`; their mean is `3.5 * g`. Thus `-g 0.2041` gives mean recovery time 0.71435. The printed theoretical diagnostics use exponential-recovery formulas and do not describe this Gamma-recovery model exactly.

### TB and SD5 with the same seed

If a baseline outbreak dies out, rerun the baseline command without `--seed` to sample another seed. Once a baseline outbreak takes off, record its printed seed. Set `SEED` below to that value; `12345` is only an example.

```sh
SEED=12345

./simulation_final --graph "$GRAPH" --output-dir "$OUT" \
  -n 10000000 -b 0.1 -g 0.2041 -T 0.0003 -H 1 -c 5 \
  --max-time 100 --seed "$SEED" --travel-ban

./simulation_final --graph "$GRAPH" --output-dir "$OUT" \
  -n 10000000 -b 0.1 -g 0.2041 -T 0.0003 -H 1 -c 5 \
  --max-time 100 --seed "$SEED" --social-distancing 5
```

Both interventions activate when cumulative infections (`I1 + R1`) reach 1% of community 1. TB stops departures and returns. SD5 reduces the transmission rate to `beta / 5` for contacts physically located in community 2. Use the same graph, seed, and other parameters for all three scenarios.

Output filenames encode the parameters and seed, with `_TB` or `_SD5` identifying the intervention. CSV columns are `time,S1,I1,R1,S2,I2,R2`; community counts refer to individuals' home communities.

## Visualization

Install the Python dependencies and a working LaTeX installation accessible on `PATH`:

```sh
python3 -m pip install numpy pandas matplotlib
python3 plot_sims.py --data-dir final_submission_sims/copenhagen_sims
```

The script creates a linear-percentage PDF and a log-count PDF for each available scenario in `DATA_DIR/figures`. Use `--output-dir PATH` to choose another destination. Linear plots omit the red intervention-time line; log-count intervention plots retain it.

The script selects the newest CSV independently for the baseline, TB, and SD scenarios. Keep the intended matching runs as the newest files in the directory, or copy just those three CSVs into a separate plotting directory. Preserve their generated filenames.

Optional animated GIFs require Pillow:

```sh
python3 -m pip install pillow
python3 plot_sims.py --data-dir final_submission_sims/copenhagen_sims --animations
```

For server runs, copy the resulting CSVs back to the local machine to plot them; visualization does not require the graph file.
