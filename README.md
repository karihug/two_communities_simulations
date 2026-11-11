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

Lognormal graphs are drawn according to the description given in Section 6A of https://arxiv.org/pdf/2510.08895. The Copenhagen graph uses data from https://figshare.com/articles/dataset/The_Copenhagen_Networks_Study_interaction_data/7267433/1.

```sh
./graphgen ln total_population_size mu sigma \
  > output_dir/ln_graph_20m.txt
```

### Degree-assortative and degree-disassortative graphs

Degree-(dis)assortative graphs are drawn according to the descriptions given in Sections 6B and 6B of https://arxiv.org/pdf/2510.08895.

```sh
./graphgen assortative total_population_size mu sigma \
  > output_dir/assortative_graph_20m.txt

./graphgen disassortative 20000000 mu population \
  > output_dir/disassortative_graph_20m.txt
```
Graph generation is stochastic and outputs from this code may slightly differ from the networks generated in the paper.

## Epidemic simulation

The simulator takes the **population per community**: `-n 10000000` requires a graph containing 20 million vertices. Each run starts with one infected individual in community 1.

Change `GRAPH` and `OUT` to run specific networks.

```sh
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
| `--max-time` | Simulation time limit |
| `--seed` | Epidemic random seed; generated and printed when omitted |
| `--travel-ban` | Enable a travel ban |
| `--social-distancing 5` | Enable social distancing with strength 5 |


### TB and SD

```sh
./simulation_final --graph "$GRAPH" --output-dir "$OUT" \
  -n 10000000 -b 0.1 -g 0.2041 -T 0.0003 -H 1 -c 5 \
  --max-time 100 --travel-ban

./simulation_final --graph "$GRAPH" --output-dir "$OUT" \
  -n 10000000 -b 0.1 -g 0.2041 -T 0.0003 -H 1 -c 5 \
  --max-time 100 --social-distancing 5
```

Both interventions activate when cumulative infections (`I1 + R1`) reach 1% of community 1. TB stops departures and returns. SD5 reduces the transmission rate to `beta / 5` for contacts physically located in community 2. 

Output filenames encode the parameters and seed, with `_TB` or `_SD5` identifying the intervention. CSV columns are `time,S1,I1,R1,S2,I2,R2`; community counts refer to individuals' home communities.

## Visualization

Install the Python dependencies and a working LaTeX installation accessible on `PATH`:

```sh
python3 -m pip install numpy pandas matplotlib
python3 plot_sims.py --data-dir final_submission_sims/copenhagen_sims
```
