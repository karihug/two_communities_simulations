#include "SIRSimulation.h"
#include <cmath>
#include <iostream>
SIRSimulation::SIRSimulation(int n_, double beta_, double gamma_,
                             double rho_T_, double rho_H_,
                             double c_, double record_interval_,
                             bool enable_travel_ban_,
                             double social_distancing_strength_,
                             unsigned seed,
                             const char* graph_file)
  : n(n_), beta(beta_), gamma(gamma_),
    rho_T(rho_T_), rho_H(rho_H_),
    prob_travel(rho_T_/(rho_T_+rho_H_)),
    c(c_), record_interval(record_interval_),
    time(0.0), last_record_time(0.0),
    enable_travel_ban(enable_travel_ban_),
    travel_ban_active(false),
    travel_ban_time(std::numeric_limits<double>::infinity()),
    social_distancing_strength(social_distancing_strength_),
    social_distancing_active(false),
    social_distancing_time(std::numeric_limits<double>::infinity()),
    intervention_threshold(static_cast<int>(std::ceil(0.01 * n_)))
{
    generator.seed(seed);

    int total = 2 * n;
    states.assign(total, 0);
    is_traveling.assign(total, false);
    infection_times.assign(total, std::numeric_limits<double>::infinity());
    if (graph_file != nullptr) {
        load_graph_from_file(graph_file);
    } else {
        generate_contact_network();
    }
    init_counts();

    for (int i = 0; i < total; i++)
        is_traveling[i] = (random_uniform() < prob_travel);

    lambda_val = c*beta - beta - gamma;
    R0         = c*beta/(beta+gamma);
    tau_scale  = (lambda_val>0) ? std::log(n)/lambda_val
                                : std::numeric_limits<double>::infinity();
    if (R0 > 1) {
        double r = 0.5;
        for (int it=0; it<1000; it++) {
            double nr = 1 - std::exp(-R0*r);
            if (std::abs(nr-r)<1e-6) { r=nr; break; }
            r = nr;
        }
        r_infty = r;
    } else r_infty = 0.0;

    std::cout << "Params: n="<<n<<" beta="<<beta
              <<" gamma="<<gamma<<" rho_T="<<rho_T
              <<" rho_H="<<rho_H<<" c="<<c<<"\n";
    std::cout << "R0="<<R0<<" tau="<<tau_scale
              <<" r_inf="<<r_infty<<"\n";

    // Count initially traveling individuals in each community
    int traveling_infectious_community1 = 0;
    int traveling_infectious_community2 = 0;
    for (int i = 0; i < n; ++i) {
        if (is_traveling[i]) {
            traveling_infectious_community1++;
        }
    }
    for (int i = n; i < 2 * n; ++i) {
        if (is_traveling[i]) {
            traveling_infectious_community2++;
        }
    }
    std::cout << "Traveling infectious individuals in Community 1: "
              << traveling_infectious_community1 << "\n";
    std::cout << "Traveling infectious individuals in Community 2: "
              << traveling_infectious_community2 << "\n";
}

void SIRSimulation::init_counts() {
    S1_count = n;  I1_count = 0;  R1_count = 0;
    S2_count = n;  I2_count = 0;  R2_count = 0;
}

void SIRSimulation::on_infection(int individual) {
    if (individual < n) {
        --S1_count; ++I1_count;
    } else {
        --S2_count; ++I2_count;
    }
}

void SIRSimulation::on_recovery(int individual) {
    if (individual < n) {
        --I1_count; ++R1_count;
    } else {
        --I2_count; ++R2_count;
    }
}

bool SIRSimulation::is_in_physical_community_one(int individual) {
    return (individual < n && !is_traveling[individual]) ||
           (individual >= n && is_traveling[individual]);
}

bool SIRSimulation::is_in_physical_community_two(int individual) {
    return (individual < n && is_traveling[individual]) ||
           (individual >= n && !is_traveling[individual]);
}

double SIRSimulation::random_uniform() {
    std::uniform_real_distribution<double> dist(0.0,1.0);
    return dist(generator);
}

double SIRSimulation::random_exponential(double rate) {
    std::exponential_distribution<double> dist(rate);
    return dist(generator);
}

double SIRSimulation::random_gamma(double scale) {
    std::gamma_distribution<double> dist(3.5, scale);
    return dist(generator);
}

void SIRSimulation::generate_contact_network() {
    int total = 2*n;
    neighbor_list.assign(total, {});
    
    // Negative binomial parameters
    double mu = 12.0;
    double r_param = 2;
    double p_param = static_cast<double>(r_param) / (r_param + mu);
    
    std::negative_binomial_distribution<int> deg_dist(r_param, p_param);
    
    std::vector<int> degrees(total);
    for (int i = 0; i < total; i++) {
        degrees[i] = deg_dist(generator);
    }
    
    std::vector<int> stubs;
    for (int i = 0; i < total; i++) {
        for (int j = 0; j < degrees[i]; j++) {
            stubs.push_back(i);
        }
    }

    if (stubs.size() % 2 == 1) {
        stubs.pop_back();  // Remove one stub
    }

    std::shuffle(stubs.begin(), stubs.end(), generator);

    for (size_t i = 0; i < stubs.size(); i += 2) {
        int u = stubs[i];
        int v = stubs[i + 1];
        if (u == v && i + 2 < stubs.size()) {
            std::swap(stubs[i + 1], stubs[i + 2]);
            v = stubs[i + 1];
        }
        if (u != v) {
            neighbor_list[u].push_back(v);
            neighbor_list[v].push_back(u);
        }
    }
}

void SIRSimulation::load_graph_from_file(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open graph file: " << filename << std::endl;
        exit(1);
    }
    
    int file_n;
    file >> file_n;
    
    if (file_n != 2 * n) {
        std::cerr << "Error: Graph has " << file_n << " nodes but simulation expects " << (2*n) << std::endl;
        exit(1);
    }
    
    neighbor_list.assign(2 * n, {});
    
    int u, v;
    while (file >> u >> v) {
        neighbor_list[u].push_back(v);
        neighbor_list[v].push_back(u);
    }
    
    file.close();
}

void SIRSimulation::schedule_event(double t,
                                   const std::string &type,
                                   int ind, int from_ind) {
    event_queue.push(Event{t,type,ind,from_ind});
}

void SIRSimulation::schedule_travel_events() {
    int total=2*n;
    for(int i=0;i<total;i++){
        if(is_traveling[i])
            schedule_event(time+random_exponential(rho_H),
                           "return_home",i,-1);
        else
            schedule_event(time+random_exponential(rho_T),
                           "travel",i,-1);
    }
}

void SIRSimulation::schedule_infections(int src) {
    for (int nb: neighbor_list[src]) {
        if (states[nb]==0) {
            double delay = random_exponential(beta);
            // Propose at rate beta; apply location and intervention rules at event time.
            schedule_event(time + delay, "infection", nb, src);
        }
    }
}

void SIRSimulation::initialize_epidemic() {
    std::uniform_int_distribution<int> pick(0,n-1);
    
    // Keep picking until we find a node with at least one neighbor
    int seed;
    do {
        seed = pick(generator);
    } while (neighbor_list[seed].empty());
    
    std::cout << "Initial infected node: " << seed 
              << " (degree: " << neighbor_list[seed].size() << ")\n";
    
    states[seed]=1;
    infection_times[seed]=time;
    schedule_event(time+random_gamma(gamma), // Gamma recovery time; gamma is the scale.
                   "recovery",seed,-1);
    schedule_travel_events();
    schedule_infections(seed);
    on_infection(seed);
}

void SIRSimulation::get_counts(int &S1,int &I1,int &R1,
                               int &S2,int &I2,int &R2)
{
    S1=I1=R1=0; S2=I2=R2=0;
    for(int i=0;i<n;i++){
        if(states[i]==0) S1++;
        else if(states[i]==1) I1++;
        else R1++;
    }
    for(int i=n;i<2*n;i++){
        if(states[i]==0) S2++;
        else if(states[i]==1) I2++;
        else R2++;
    }
}

void SIRSimulation::run(double max_time) {
    while(!event_queue.empty() && time<max_time) {
        auto ev = event_queue.top();
        event_queue.pop();

        time = ev.time;

        if (ev.type=="infection" && states[ev.individual]==0) {
            if ((ev.from_individual>=0 && states[ev.from_individual]==1)
                && states[ev.individual]==0)
            {
                int i=ev.individual, j=ev.from_individual;
                bool same_comm = (i<n&&j<n)||(i>=n&&j>=n);
                if ((same_comm && is_traveling[i]==is_traveling[j])||
                    (!same_comm && is_traveling[i]!=is_traveling[j]))
                {
                    // Thin rate-beta attempts using the intervention in effect now,
                    // including attempts queued before activation or a location change.
                    bool is_infecting_community_two_under_sd =
                        (is_in_physical_community_two(i) &&
                         is_in_physical_community_two(j) &&
                         social_distancing_active);
                    if (!is_infecting_community_two_under_sd ||
                        random_uniform() < 1.0 / social_distancing_strength)
                    {
                        states[i]=1;
                        infection_times[i]=time;
                        on_infection(i);
                        schedule_event(time+random_gamma(gamma), // Gamma recovery time; gamma is the scale.
                                       "recovery",i,-1);
                        schedule_infections(i);
                    }
                    if (enable_travel_ban &&
                        !travel_ban_active &&
                        I1_count + R1_count >= intervention_threshold)
                    {
                        travel_ban_active = true;
                        travel_ban_time   = time;
                        std::cout << "*** Travel ban ON at t=" << time << " ***\n";

                        // Count how many infectious individuals are traveling from each community
                        int traveling_infectious_community1 = 0;
                        int traveling_infectious_community2 = 0;
                        for (int k = 0; k < n; ++k) {
                            if (states[k] == 1 && is_traveling[k]) {
                                traveling_infectious_community1++;
                            }
                        }
                        for (int k = n; k < 2 * n; ++k) {
                            if (states[k] == 1 && is_traveling[k]) {
                                traveling_infectious_community2++;
                            }
                        }
                        std::cout << "Traveling infectious individuals in Community 1: "
                                << traveling_infectious_community1
                                << " out of " << I1_count << "\n";
                        std::cout << "Traveling infectious individuals in Community 2: "
                                << traveling_infectious_community2
                                << " out of " << I2_count << "\n";
                    }
                    if (!social_distancing_active && social_distancing_strength > 1 &&
                        I1_count + R1_count >= intervention_threshold)
                    {
                        social_distancing_active = true;
                        social_distancing_time   = time;
                        std::cout << "*** Social distancing ON at t=" << time << " ***\n";
                    }
                }
                // A failed attempt does not remove future transmission opportunities.
                // Stop retrying once the target is infected or the source has recovered.
                if (states[i]==0) {
                    schedule_event(time + random_exponential(beta), "infection", i, j);
                }
            }
        }
        else if (ev.type=="recovery") {
            states[ev.individual]=2;
            on_recovery(ev.individual);
        }
        else if (ev.type=="travel" && !(enable_travel_ban &&
                                        travel_ban_active)) {
            is_traveling[ev.individual]=true;
            if(states[ev.individual]!=2)
                schedule_event(time+random_exponential(rho_H),
                               "return_home",ev.individual,-1);
        }
        else if (ev.type=="return_home" && !(enable_travel_ban &&
                                             travel_ban_active)) {
            is_traveling[ev.individual]=false;
            if(states[ev.individual]!=2)
                schedule_event(time+random_exponential(rho_T),
                               "travel",ev.individual,-1);
        }

        if (time - last_record_time >= record_interval) {
            history.push_back(
              HistoryEntry{time, 
                           S1_count, I1_count, R1_count,
                           S2_count, I2_count, R2_count});
            last_record_time=time;
        }
    }
}

void SIRSimulation::output_history(const std::string &fname) {
    std::ofstream out(fname);
    if (!out.is_open()) {
        std::cerr<<"Error opening "<<fname<<"\n";
        return;
    }
    out<<"time,S1,I1,R1,S2,I2,R2\n";
    for (auto &h:history) {
        out<<h.time<<","
           <<h.S1<<","<<h.I1<<","<<h.R1<<","
           <<h.S2<<","<<h.I2<<","<<h.R2<<"\n";
    }
}

void SIRSimulation::print_final_counts() {
    int S1,I1,R1,S2,I2,R2;
    get_counts(S1,I1,R1,S2,I2,R2);
    std::cout<<"Com1: S="<<S1<<" I="<<I1<<" R="<<R1<<"\n"
             <<"Com2: S="<<S2<<" I="<<I2<<" R="<<R2<<"\n";
}
