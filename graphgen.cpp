#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <array>
#include <cstdint>
#include <iomanip>
#include <string>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <stdexcept>
using namespace std;
using ll = long long;

#define BETA 1.8
#define GAMMA 1.5

static mt19937_64 rng(chrono::steady_clock::now().time_since_epoch().count());
static std::exponential_distribution<> infect_clock(BETA);
static std::exponential_distribution<> recover_clock(GAMMA);

// split s by delim into tokens
vector<string> split(const string &s, char delim) {
    vector<string> out;
    string t;
    istringstream in(s);
    while (getline(in, t, delim)) if (!t.empty()) out.push_back(t);
    return out;
}

void gen_recovery_clocks(int n) {

    for (int i = 0; i < n; i++) {
        double t = recover_clock(rng);
        cout << t << "\n";
    }
}

void gen_ER(int n, double p) {
    cout << n << "\n";

    gen_recovery_clocks(n);

    for (int i = 0; i < n - 1; i++) {
        int R = (n - 1) - i;
        // Set binomial to Binomial(R, p) and sample d
        binomial_distribution<int> binom(R, p);
        int d = binom(rng);
        if (d == 0) continue;

        // Sample d unique neighbors from [i+1 .. n-1]
        unordered_set<int> chosen;
        chosen.reserve(d * 2);
        uniform_int_distribution<int> pick(i+1, n-1);
        while ((int)chosen.size() < d) {
            chosen.insert(pick(rng));
        }
        for (int j : chosen) {
            double itoj = infect_clock(rng);
            double jtoi = infect_clock(rng);
            cout << i << " " << j << " " << itoj << " " << jtoi << "\n";
        }

    }
}

// SBM: sizes[K], mat[K][K]
// assume mat[a][b] is probability for any i in block a, j in block b (i<j)
void gen_SBM(int n, const vector<int> &sizes,
             const vector<vector<double>> &mat) {
    int K = sizes.size();
    vector<int> start(K);
    start[0] = 0;
    for (int a = 1; a < K; a++) start[a] = start[a-1] + sizes[a-1];
    uniform_real_distribution<double> D(0,1);
    cout << n << "\n";

    gen_recovery_clocks(n);

    for (int a = 0; a < K; a++) {
        for (int b = a; b < K; b++) {
            double p = mat[a][b];
            int i0 = start[a], i1 = start[a] + sizes[a];
            int j0 = start[b], j1 = start[b] + sizes[b];
            for (int i = i0; i < i1; i++) {
                int jbeg = (a==b ? i+1 : j0);
                for (int j = jbeg; j < j1; j++) {
                    if (D(rng) < p) {
                        double itoj = infect_clock(rng);
                        double jtoi = infect_clock(rng);
                        cout << i << " " << j << " " << itoj << " " << jtoi << "\n";
                    }
                }
            }
        }
    }
}

// PA: Barabási–Albert, start with m0 = m nodes fully connected, then
// each new node attaches to m existing nodes by degree-proportional sampling
void gen_PA(int n, int m) {
    if (n <= m) {
        cerr << "Need n > m\n";
        return;
    }
    cout << n << "\n";

    gen_recovery_clocks(n);

    vector<pair<int,int>> edges;
    vector<int> deg(n,0), pool;
    auto add_edge = [&](int u,int v){
        edges.emplace_back(u,v);
        deg[u]++; deg[v]++;
        pool.push_back(u);
        pool.push_back(v);
    };
    // initial clique of nodes [0..m-1]
    for (int i = 0; i < m; i++)
        for (int j = i+1; j < m; j++)
            add_edge(i,j);
    // attach nodes m..n-1
    uniform_int_distribution<ll> D;
    for (int u = m; u < n; u++) {
        unordered_set<int> chosen;
        int need = m;
        while (chosen.size() < (size_t)m) {
            int idx = D(rng) % pool.size();
            int v = pool[idx];
            if (v == u || chosen.count(v)) continue;
            chosen.insert(v);
        }
        for (int v : chosen) add_edge(u,v);
    }
    for (auto [u,v] : edges) {
        double vtou = infect_clock(rng);
        double utov = infect_clock(rng);
        cout << u << " " << v << " " << utov << " " << vtou << "\n";
    }
}

// NB: Negative binomial configuration model
void gen_NB(int n, double mean_k, double disp_r) {
    cout << n << "\n";

    // Gamma–Poisson mixture: NB with mean = mean_k, dispersion = disp_r
    double shape = disp_r;               // r
    double scale = mean_k / disp_r;      // so E[lambda] = mean_k

    std::gamma_distribution<double> gamma(shape, scale);

    vector<int> degrees(n);
    for (int i = 0; i < n; ++i) {
        double lambda = std::max(0.0, gamma(rng));
        std::poisson_distribution<int> pois(lambda);
        degrees[i] = std::max(0, pois(rng));
    }

    // Configuration model stubs
    vector<int> stubs;
    stubs.reserve(std::accumulate(degrees.begin(), degrees.end(), 0LL));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < degrees[i]; ++j) {
            stubs.push_back(i);
        }
    }
    if (stubs.size() % 2 == 1) stubs.pop_back();

    std::shuffle(stubs.begin(), stubs.end(), rng);

    // Remove loops and duplicate edges.
    auto encode_edge = [](int u, int v) -> ll {
        if (u > v) std::swap(u, v);
        return ( (ll)u << 32 ) | (unsigned int)v;
    };

    std::unordered_set<ll> edge_set;
    edge_set.reserve(stubs.size());

    for (size_t i = 0; i + 1 < stubs.size(); i += 2) {
        int u = stubs[i];
        int v = stubs[i + 1];

        if (u == v) continue; // remove loops

        ll code = encode_edge(u, v);
        if (edge_set.count(code)) continue; // remove multiedges

        edge_set.insert(code);
        cout << u << " " << v << "\n";
    }
}

// Write a simple undirected edge list: remove loops and repeated pairs.
void write_simple_edges(vector<pair<int, int>>& edges) {
    for (auto& edge : edges) {
        if (edge.first > edge.second) swap(edge.first, edge.second);
    }
    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());

    for (const auto& edge : edges) {
        if (edge.first != edge.second) {
            cout << edge.first << " " << edge.second << "\n";
        }
    }
}

// Two equally sized communities. n is the TOTAL number of vertices.
// Each vertex gets independent home and travel degrees from the same law.
void gen_LN(int n, double mu, double sigma) {
    const int community_size = n / 2;
    cout << n << "\n";
    lognormal_distribution<double> degree_distribution(mu, sigma);

    vector<int> travel_stubs[2];

    // Build the home edges in each community separately.
    // Community 0 has IDs [0, community_size); community 1 has the remaining IDs.
    for (int community = 0; community < 2; ++community) {
        vector<int> home_stubs;
        const int first_vertex = community * community_size;

        for (int vertex = first_vertex; vertex < first_vertex + community_size; ++vertex) {
            // Floor both draws. Bound them by the number of possible neighbors.
            int home_degree = static_cast<int>(min(
                floor(degree_distribution(rng)), double(community_size - 1)));
            int travel_degree = static_cast<int>(min(
                floor(degree_distribution(rng)), double(community_size)));

            home_stubs.insert(home_stubs.end(), home_degree, vertex);
            travel_stubs[community].insert(travel_stubs[community].end(), travel_degree, vertex);
        }

        shuffle(home_stubs.begin(), home_stubs.end(), rng);
        // Pairing requires an even total; drop a uniformly chosen stub if odd.
        if (home_stubs.size() % 2 == 1) home_stubs.pop_back();

        vector<pair<int, int>> home_edges;
        home_edges.reserve(home_stubs.size() / 2);
        for (size_t i = 0; i < home_stubs.size(); i += 2) {
            home_edges.emplace_back(home_stubs[i], home_stubs[i + 1]);
        }
        write_simple_edges(home_edges);
    }

    // Travel edges always connect opposite communities. Independent draws
    // need not give equal stub totals, so randomly discard the excess stubs.
    shuffle(travel_stubs[0].begin(), travel_stubs[0].end(), rng);
    shuffle(travel_stubs[1].begin(), travel_stubs[1].end(), rng);
    const size_t travel_edge_count = min(travel_stubs[0].size(), travel_stubs[1].size());
    cerr << "Unmatched travel stubs discarded: community 1 = "
         << travel_stubs[0].size() - travel_edge_count
         << ", community 2 = " << travel_stubs[1].size() - travel_edge_count << "\n";

    vector<pair<int, int>> travel_edges;
    travel_edges.reserve(travel_edge_count);
    for (size_t i = 0; i < travel_edge_count; ++i) {
        travel_edges.emplace_back(travel_stubs[0][i], travel_stubs[1][i]);
    }
    write_simple_edges(travel_edges);
}

// Rows and columns are high, medium, low degree classes.
using MixingMatrix = array<array<double, 3>, 3>;

// A shuffled ranking gives reproducible random tie-breaking at the cutoffs.
vector<unsigned char> classify_degrees(const vector<int>& degrees, bool assortative) {
    vector<int> ranking(degrees.size());
    iota(ranking.begin(), ranking.end(), 0);
    shuffle(ranking.begin(), ranking.end(), rng);
    stable_sort(ranking.begin(), ranking.end(), [&](int a, int b) {
        return degrees[a] > degrees[b];
    });
    size_t high_end = degrees.size() * (assortative ? 15 : 5) / 100;
    size_t medium_end = degrees.size() * (assortative ? 50 : 20) / 100;
    vector<unsigned char> classes(degrees.size(), 2);
    for (size_t i = 0; i < high_end; ++i) classes[ranking[i]] = 0;
    for (size_t i = high_end; i < medium_end; ++i) classes[ranking[i]] = 1;
    return classes;
}

// Metropolis-Hastings correction: S[a] P[a][b] = S[b] P[b][a].
// S counts stubs, not vertices. Rejected probability stays in the same class.
MixingMatrix balance_mixing_matrix(MixingMatrix proposal, const array<uint64_t, 3>& stubs) {
    for (auto& row : proposal) {
        double sum = accumulate(row.begin(), row.end(), 0.0);
        for (double& probability : row) probability /= sum;
    }
    MixingMatrix balanced{};
    for (int a = 0; a < 3; ++a) {
        balanced[a][a] = 1.0;
        if (stubs[a] == 0) continue;
        for (int b = 0; b < 3; ++b) {
            if (a == b) continue;
            balanced[a][b] = min(proposal[a][b], double(stubs[b]) * proposal[b][a] / stubs[a]);
            balanced[a][a] -= balanced[a][b];
        }
    }
    cerr << "Adjusted class matrix (H, M, L):\n";
    for (const auto& row : balanced) cerr << row[0] << ' ' << row[1] << ' ' << row[2] << '\n';
    return balanced;
}

// A compact checksum lets an independent edge-file scan check all vertex counts.
uint64_t degree_checksum(const vector<uint32_t>& degrees, int begin, int end) {
    uint64_t checksum = 14695981039346656037ULL;
    for (int vertex = begin; vertex < end; ++vertex) {
        checksum ^= degrees[vertex];
        checksum *= 1099511628211ULL;
    }
    return checksum;
}

// SI assortative/disassortative construction: ONE total degree per vertex,
// with independent stub preferences of 80% home and 20% travel.
void gen_mixed_LN(int n, double mu, double sigma, bool assortative) {
    const int community_size = n / 2;
    lognormal_distribution<double> sample_degree(mu, sigma);
    vector<int> degrees(n);
    for (int& degree : degrees) degree = static_cast<int>(min(floor(sample_degree(rng)), double(n - 1)));
    vector<unsigned char> classes = classify_degrees(degrees, assortative);

    array<uint64_t, 3> class_stubs{};
    array<uint64_t, 3> class_vertices{};
    for (int vertex = 0; vertex < n; ++vertex) {
        class_stubs[classes[vertex]] += degrees[vertex];
        ++class_vertices[classes[vertex]];
    }
    MixingMatrix proposal = assortative
        ? MixingMatrix{{{0.939, 0.0516, 0.009}, {0.1, 0.865, 0.035}, {0.05, 0.1, 0.85}}}
        : MixingMatrix{{{0.056, 0.08, 0.864}, {0.091, 0.822, 0.087}, {0.9, 0.08, 0.02}}};
    cerr << setprecision(12);
    cerr << "Class vertex counts: " << class_vertices[0] << ' ' << class_vertices[1] << ' ' << class_vertices[2] << '\n';
    cerr << "Class stub counts: " << class_stubs[0] << ' ' << class_stubs[1] << ' ' << class_stubs[2] << '\n';
    MixingMatrix mixing = balance_mixing_matrix(proposal, class_stubs);

    // Six groups: community 1 H/M/L, then community 2 H/M/L.
    // buckets[source group][requested destination group] holds vertex IDs.
    vector<int> buckets[6][6];
    vector<uint32_t> assigned_home(n, 0), assigned_travel(n, 0);
    bernoulli_distribution choose_home(0.8);
    uniform_real_distribution<double> uniform(0.0, 1.0);
    for (int vertex = 0; vertex < n; ++vertex) {
        int community = vertex / community_size;
        int degree_class = classes[vertex];
        for (int stub = 0; stub < degrees[vertex]; ++stub) {
            bool home = choose_home(rng);
            int destination_community = home ? community : 1 - community;
            double draw = uniform(rng);
            int destination_class = 0;
            double cumulative = mixing[degree_class][0];
            while (destination_class < 2 && draw >= cumulative) {
                cumulative += mixing[degree_class][++destination_class];
            }
            buckets[3 * community + degree_class][3 * destination_community + destination_class].push_back(vertex);
            if (home) ++assigned_home[vertex];
            else ++assigned_travel[vertex];
        }
    }
    cerr << "Sampled and assigned all stubs.\n";

    uint64_t total_stubs = accumulate(class_stubs.begin(), class_stubs.end(), uint64_t(0));
    vector<pair<int, int>> edges;
    edges.reserve(total_stubs / 2);
    vector<int> remaining_home[2], remaining_travel[2];

    // Pair reciprocal requests. Leftovers retain their community preference.
    for (int source = 0; source < 6; ++source) {
        for (int destination = source; destination < 6; ++destination) {
            auto& forward = buckets[source][destination];
            shuffle(forward.begin(), forward.end(), rng);
            if (source == destination) {
                size_t paired = forward.size() - forward.size() % 2;
                for (size_t i = 0; i < paired; i += 2) edges.emplace_back(forward[i], forward[i + 1]);
                if (paired < forward.size()) remaining_home[source / 3].push_back(forward.back());
            } else {
                auto& reverse = buckets[destination][source];
                shuffle(reverse.begin(), reverse.end(), rng);
                size_t paired = min(forward.size(), reverse.size());
                for (size_t i = 0; i < paired; ++i) edges.emplace_back(forward[i], reverse[i]);
                auto& leftovers = (source / 3 == destination / 3) ? remaining_home : remaining_travel;
                leftovers[source / 3].insert(leftovers[source / 3].end(), forward.begin() + paired, forward.end());
                leftovers[destination / 3].insert(leftovers[destination / 3].end(), reverse.begin() + paired, reverse.end());
                vector<int>().swap(reverse);
            }
            vector<int>().swap(forward);
        }
    }

    // Uniformly pair remaining home stubs within each community.
    for (int community = 0; community < 2; ++community) {
        auto& stubs = remaining_home[community];
        cerr << "Community " << community + 1 << " fallback home stubs: " << stubs.size() << '\n';
        shuffle(stubs.begin(), stubs.end(), rng);
        for (size_t i = 0; i + 1 < stubs.size(); i += 2) edges.emplace_back(stubs[i], stubs[i + 1]);
        cerr << "Community " << community + 1 << " unmatched home stubs: " << stubs.size() % 2 << '\n';
    }
    // Uniformly pair remaining travel stubs across communities; discard excess.
    for (auto& stubs : remaining_travel) shuffle(stubs.begin(), stubs.end(), rng);
    size_t paired = min(remaining_travel[0].size(), remaining_travel[1].size());
    for (size_t i = 0; i < paired; ++i) edges.emplace_back(remaining_travel[0][i], remaining_travel[1][i]);
    for (int community = 0; community < 2; ++community) {
        cerr << "Community " << community + 1 << " fallback travel stubs: " << remaining_travel[community].size() << '\n';
        cerr << "Community " << community + 1 << " unmatched travel stubs: " << remaining_travel[community].size() - paired << '\n';
    }

    // Simplify the complete graph, including edges from fallback pairing.
    for (auto& edge : edges) if (edge.first > edge.second) swap(edge.first, edge.second);
    size_t before = edges.size();
    edges.erase(remove_if(edges.begin(), edges.end(), [](const auto& edge) {
        return edge.first == edge.second;
    }), edges.end());
    cerr << "Self-loops removed: " << before - edges.size() << '\n';
    before = edges.size();
    sort(edges.begin(), edges.end());
    edges.erase(unique(edges.begin(), edges.end()), edges.end());
    cerr << "Duplicate edges removed: " << before - edges.size() << '\n';

    // Count realized degrees from the FINAL edges, after all removals.
    vector<uint32_t> home(n, 0), travel(n, 0);
    for (const auto& edge : edges) {
        auto& counts = (edge.first / community_size == edge.second / community_size) ? home : travel;
        ++counts[edge.first];
        ++counts[edge.second];
    }
    for (int vertex = 0; vertex < n; ++vertex) {
        if (home[vertex] > assigned_home[vertex] || travel[vertex] > assigned_travel[vertex]
            || assigned_home[vertex] + assigned_travel[vertex] != uint32_t(degrees[vertex])) {
            throw runtime_error("Vertex degree accounting failed");
        }
    }
    cerr << "Per-vertex degree accounting passed for " << n << " vertices.\n";
    for (int community = 0; community < 2; ++community) {
        int begin = community * community_size, end = begin + community_size;
        uint64_t home_sum = accumulate(home.begin() + begin, home.begin() + end, uint64_t(0));
        uint64_t travel_sum = accumulate(travel.begin() + begin, travel.begin() + end, uint64_t(0));
        cerr << "Community " << community + 1 << " means: home=" << double(home_sum) / community_size
             << " travel=" << double(travel_sum) / community_size
             << " total=" << double(home_sum + travel_sum) / community_size << '\n';
        cerr << "Community " << community + 1 << " checksums: home=" << degree_checksum(home, begin, end)
             << " travel=" << degree_checksum(travel, begin, end) << '\n';
    }
    long double endpoint_sum = 0, endpoint_square_sum = 0, product_sum = 0;
    for (const auto& edge : edges) {
        long double a = home[edge.first] + travel[edge.first];
        long double b = home[edge.second] + travel[edge.second];
        endpoint_sum += a + b;
        endpoint_square_sum += a * a + b * b;
        product_sum += a * b;
    }
    if (!edges.empty()) {
        long double mean = endpoint_sum / (2 * edges.size());
        long double variance = endpoint_square_sum / (2 * edges.size()) - mean * mean;
        cerr << "Degree assortativity: " << (product_sum / edges.size() - mean * mean) / variance << '\n';
    }
    cerr << "Writing " << edges.size() << " edges.\n";
    cout << n << '\n';
    for (const auto& edge : edges) cout << edge.first << ' ' << edge.second << '\n';
}

// PA with user-specified degree distribution
void gen_PA_custom(int n, const vector<double>& degree_dist) {
    if (n <= 1) {
        cerr << "Need n > 1\n";
        return;
    }

    // Validate degree distribution
    double sum_prob = accumulate(degree_dist.begin(), degree_dist.end(), 0.0);
    if (abs(sum_prob - 1.0) > 1e-6) {
        cerr << "Degree distribution probabilities must sum to 1\n";
        return;
    }

    // Precompute cumulative distribution for sampling
    vector<double> cumulative_dist(degree_dist.size());
    partial_sum(degree_dist.begin(), degree_dist.end(), cumulative_dist.begin());

    auto sample_degree = [&]() {
        double r = uniform_real_distribution<double>(0, 1)(rng);
        return lower_bound(cumulative_dist.begin(), cumulative_dist.end(), r) - cumulative_dist.begin();
    };

    cout << n << "\n";

    gen_recovery_clocks(n);

    vector<pair<int, int>> edges;
    vector<int> deg(n, 0), pool;

    auto add_edge = [&](int u, int v) {
        edges.emplace_back(u, v);
        deg[u]++;
        deg[v]++;
        pool.push_back(u);
        pool.push_back(v);
    };

    // Use the length of degree_dist as m0 (maximum degree)
    int m0 = degree_dist.size();

    // Initial clique of nodes [0..m0-1]
    for (int i = 0; i < m0; i++) {
        for (int j = i + 1; j < m0; j++) {
            add_edge(i, j);
        }
    }

    // Attach nodes m0..n-1
    uniform_int_distribution<ll> D;
    for (int u = m0; u < n; u++) {
        unordered_set<int> chosen;
        int m = sample_degree(); // Sample degree for the new node
        while (chosen.size() < (size_t)m) {
            int idx = D(rng) % pool.size();
            int v = pool[idx];
            if (v == u || chosen.count(v)) continue;
            chosen.insert(v);
        }
        for (int v : chosen) add_edge(u, v);
    }

    for (auto [u, v] : edges) {
        double vtou = infect_clock(rng);
        double utov = infect_clock(rng);
        cout << u << " " << v << " " << utov << " " << vtou << "\n";
    }
}

int main(int argc, char **argv) {
    ios::sync_with_stdio(false);
    if (argc < 2) {
        cerr << "Usage:\n"
             << "  " << argv[0] << " erdosrenyi n p\n"
             << "  " << argv[0] << " sbm n sizes mat\n"
             << "  " << argv[0] << " pa n m\n"
             << "  " << argv[0] << " nb n mu r_param\n"
             << "  " << argv[0] << " ln total_nodes mu sigma (two equal communities)\n"
             << "  " << argv[0] << " assortative|disassortative total_nodes mu sigma seed\n";
        return 1;
    }
    string mode = argv[1];
    if (mode == "erdosrenyi") {
        if (argc!=4) { cerr<<"bad args\n"; return 1; }
        int n = stoi(argv[2]);
        double c = stod(argv[3]);
        double p = c / (n-1);
        gen_ER(n,p);
    }
    else if (mode == "sbm") {
        if (argc!=5) { cerr<<"bad args\n"; return 1; }
        int n = stoi(argv[2]);
        // e.g. sizes="10,20,30"
        auto vs = split(argv[3],',');
        vector<int> sizes;
        for (auto &t:vs) sizes.push_back(stoi(t));
        int K = sizes.size();
        // e.g. mat="0.1,0.05;0.05,0.2"
        vector<string> rows = split(argv[4],';');
        if ((int)rows.size()!=K) { cerr<<"matrix size mismatch\n"; return 1; }
        vector<vector<double>> mat(K, vector<double>(K));
        for (int i = 0; i < K; i++) {
            auto cols = split(rows[i],',');
            if ((int)cols.size()!=K) { cerr<<"matrix row "<<i<<" bad\n"; return 1; }
            for (int j = 0; j < K; j++)
                mat[i][j] = stod(cols[j]);
        }
        gen_SBM(n,sizes,mat);
    }
    else if (mode == "pa") {
        if (argc != 4) { cerr << "bad args\n"; return 1; }
        string arg3 = argv[3];
        if (all_of(arg3.begin(), arg3.end(), [](char c) { return isdigit(c) || c == '-'; })) {
            // If argv[3] is an integer
            int n = stoi(argv[2]);
            int m = stoi(argv[3]);
            gen_PA(n, m);
        } else {
            // If argv[3] is a comma-separated list of doubles
            int n = stoi(argv[2]);
            auto vs = split(argv[3], ',');
            vector<double> degree_dist;
            for (auto &t : vs) degree_dist.push_back(stod(t));
            gen_PA_custom(n, degree_dist);
        }
    }
    else if (mode == "nb") {
        if (argc != 5) { cerr << "bad args\n"; return 1; }
        int n = stoi(argv[2]);
        double mu = stod(argv[3]);
        double r_param = stod(argv[4]);
        gen_NB(n, mu, r_param);
    }
    else if (mode == "ln") {
        if (argc != 5) { cerr << "bad args\n"; return 1; }
        int n = stoi(argv[2]);
        double mu = stod(argv[3]);
        double sigma = stod(argv[4]);
        if (n < 2 || n % 2 != 0 || !isfinite(mu) || !isfinite(sigma) || sigma <= 0) {
            cerr << "ln requires an even total_nodes >= 2, finite mu, and finite sigma > 0\n";
            return 1;
        }
        gen_LN(n, mu, sigma);
    }
    else if (mode == "assortative" || mode == "disassortative") {
        if (argc != 6) { cerr << "Expected total_nodes mu sigma seed\n"; return 1; }
        int n = stoi(argv[2]);
        double mu = stod(argv[3]), sigma = stod(argv[4]);
        if (n < 2 || n % 2 != 0 || !isfinite(mu) || !isfinite(sigma) || sigma <= 0) {
            cerr << "Require even total_nodes >= 2, finite mu, and finite sigma > 0\n";
            return 1;
        }
        uint64_t seed = stoull(argv[5]);
        rng.seed(seed);
        cerr << "Mode=" << mode << " total_nodes=" << n << " mu=" << mu << " sigma=" << sigma << " seed=" << seed << '\n';
        gen_mixed_LN(n, mu, sigma, mode == "assortative");
        if (!cout) { cerr << "Graph output failed\n"; return 1; }
    }
    else {
        cerr<<"unknown mode\n";
        return 1;
    }
    return 0;
}
