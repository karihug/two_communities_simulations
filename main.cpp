#include "SIRSimulation.h"
#include <iostream>
#include <sstream>
#include <chrono>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <random>
#include <filesystem>

int find_arg_idx(int argc, char** argv, const char* option) {
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], option) == 0) {
            return i;
        }
    }
    return -1;
}

int find_int_arg(int argc, char** argv, const char* option, int default_value) {
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc - 1) {
        return std::stoi(argv[iplace + 1]);
    }

    return default_value;
}

double find_double_arg(int argc, char** argv, const char* option, double default_value) {
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc - 1) {
        return std::stod(argv[iplace + 1]);
    }

    return default_value;
}

char* find_string_option(int argc, char** argv, const char* option, char* default_value) {
    int iplace = find_arg_idx(argc, argv, option);

    if (iplace >= 0 && iplace < argc - 1) {
        return argv[iplace + 1];
    }

    return default_value;
}

int main(int argc, char** argv) {
    if (find_arg_idx(argc, argv, "-h") >= 0) {
        std::cout << "Options:" << std::endl;
        std::cout << "-h: see this help" << std::endl;
        std::cout << "-n <int>: size of population" << std::endl;
        std::cout << "-s <int>: set particle initialization seed" << std::endl;
        std::cout << "-b <float>: set infection rate beta" << std::endl;
        std::cout << "-g <float>: set Gamma recovery-time scale (shape 3.5)" << std::endl;
        std::cout << "--output-dir <directory>: CSV output directory (default: lognormal_sims)" << std::endl;
        std::cout << "-T <float>: set rate of traveling rho_T" << std::endl;
        std::cout << "-H <float>: set rate of returning home rho_H" << std::endl;
        std::cout << "-c <float>: set average number of neighbors c" << std::endl;
        std::cout << "--travel-ban: enable travel ban at 1% cumulatively infected in community 1" << std::endl;
        std::cout << "--social-distancing <float>: enable social distancing at 1% cumulatively infected in community 1 with given strength" << std::endl;
        std::cout << "--seed <int>: set RNG seed (default: random_device)" << std::endl;
        std::cout << "--graph <filename>: load contact network from file" << std::endl;
        return 0;
    }

    int n = find_int_arg(argc, argv, "-n", 1000000);          // Community size (each community has n individuals, total 2*n)
    double beta = find_double_arg(argc, argv, "-b", 1.5);     // Infection rate.
    double gamma = find_double_arg(argc, argv, "-g", 3.0);    // Gamma recovery-time scale (shape 3.5).
    double rho_T = find_double_arg(argc, argv, "-T", 0.0002);// Rate of traveling.
    double rho_H = find_double_arg(argc, argv, "-H", 1.0);   // Rate of returning home.
    double c = find_double_arg(argc, argv, "-c", 5.0);       // Average number of neighbors.
    double record_interval = 0.05; // Interval for recording simulation history.
    double max_time = find_double_arg(argc, argv, "--max-time", 20.0);        // Maximum simulation time.

    bool enable_travel_ban        = (find_arg_idx(argc, argv, "--travel-ban")        >= 0);
    bool enable_social_distancing = (find_arg_idx(argc, argv, "--social-distancing") >= 0);
    double social_distancing_factor = find_double_arg(argc, argv, "--social-distancing", 1.0); // Default to 1.0 if not provided
     unsigned seed;
     int si = find_arg_idx(argc, argv, "--seed");
     if (si >= 0 && si < argc-1) {
         seed = static_cast<unsigned>(std::stoul(argv[si+1]));
     } else {
         seed = std::random_device{}();
         std::cout << "Generated random seed: " << seed << std::endl;
     }

    char* graph_file = find_string_option(argc, argv, "--graph", nullptr);
    char* output_dir_arg = find_string_option(argc, argv, "--output-dir", nullptr);
    const std::filesystem::path output_dir = output_dir_arg ? output_dir_arg : "lognormal_sims";
    std::filesystem::create_directories(output_dir);

    auto start = std::chrono::high_resolution_clock::now();
    std::cout << "Starting simulation..." << std::endl;

    SIRSimulation sim(n, beta, gamma, rho_T, rho_H, c,
                      record_interval,
                      enable_travel_ban,
                      social_distancing_factor,
                      seed,
                      graph_file);

    auto declared = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_declared = declared - start;
    std::cout << "Declared simulation. Elapsed time: " << elapsed_declared.count() << " seconds" << std::endl;

    sim.initialize_epidemic();

    auto initialized = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_initialized = initialized - declared;
    std::cout << "Initialized epidemic. Elapsed time: " << elapsed_initialized.count() << " seconds" << std::endl;

    sim.run(max_time);

    auto sim_run = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_run = sim_run - initialized;
    std::cout << "Ran simulation. Elapsed time: " << elapsed_run.count() << " seconds" << std::endl;

    sim.print_final_counts();

    std::stringstream ss;
    ss << "simulation_history_n" << n
       << "_beta" << beta
       << "_gamma" << gamma
       << "_rhoT" << rho_T
       << "_rhoH" << rho_H
       << "_c"    << c
       << "_maxTime" << max_time
       << "_seed" << seed;

    if (enable_travel_ban) {
        ss << "_TB";
    }
    if (social_distancing_factor != 1.0) {
        ss << "_SD" << social_distancing_factor;
    }

    ss << ".csv";
    std::string filename = (output_dir / ss.str()).string();

    sim.output_history(filename);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Simulation complete. History saved to " << filename << std::endl;
    std::cout << "Elapsed time: " << elapsed.count() << " seconds" << std::endl;
    return 0;
}
