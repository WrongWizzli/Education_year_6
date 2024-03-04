#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <string>
#include <vector>
#include <fstream>

#include "omp.h"

typedef std::vector<std::vector<std::vector<double>>> vector3d;
typedef std::vector<std::vector<double>> vector2d;

struct Point {
    double x;
    double y;
    double z;
};

struct Uanalytical {
    const double at = 2 * M_PI;
    const double cx = 3 * M_PI;
    const double cy = 2 * M_PI;
    const double cz = 2 * M_PI;
    Point L;

    explicit Uanalytical(const Point& L): L(L) {}

    double Value(double x, double y, double z, double t) const {
        return std::sin(cx * x / L.x) * std::sin(cy * y / L.y) * std::sin(cz * z / L.z) * std::cos(at * t);
    }

    double Dx(double x, double y, double z, double t) const {
        return cx / L.x * std::cos(cx * x / L.x) * std::sin(cy * y / L.y) * std::sin(cz * z / L.z) * std::cos(at * t);
    }

    double Dy(double x, double y, double z, double t) const {
        return cy / L.y * std::sin(cx * x / L.x) * std::cos(cy * y / L.y) * std::sin(cz * z / L.z) * std::cos(at * t);
    }

    double Dz(double x, double y, double z, double t) const {
        return cz / L.z * std::sin(cx * x / L.x) * std::sin(cy * y / L.y) * std::cos(cz * z / L.z) * std::cos(at * t);
    }
};

struct Uapproximated {
    Point h;
    double tau;
    int N, K;
    vector3d grid_n0;
    vector3d grid_n1;
    vector3d grid_n2;

    explicit Uapproximated(int N, int K, const Point& L, double T): N(N), K(K) {
        h = Point{
            .x = L.x / N,
            .y = L.y / N,
            .z = L.z / N
        };
        tau = T / K;
        grid_n0 = vector3d(N + 1, vector2d(N + 1, std::vector<double>(N + 1, 0)));
        grid_n1 = vector3d(N + 1, vector2d(N + 1, std::vector<double>(N + 1, 0)));
        grid_n2 = vector3d(N + 1, vector2d(N + 1, std::vector<double>(N + 1, 0)));
    }

    void FillZerothLayer(const Uanalytical& u) {
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < N + 1; ++i) {
            for (int j = 0; j < N + 1; ++j) {
                for (int k = 0; k < N + 1; ++k) {
                    grid_n0[i][j][k] = u.Value(h.x * i, h.y * j, h.z * k, 0);
                }
            }
        }
    }

    int SubstitutePeriodicIndex(int index) const {
        if (index < 0) return N - 1;
        if (index >= N) return 0; 
        return index;
    }

    double GetLaplasOperator(const vector3d& grid, int i, int j, int k) const {
        double x_component = -2 * grid[i][j][k];
        if (i > 0) {
            x_component += grid[i - 1][j][k];
        }
        if (i < N - 1) {
            x_component += grid[i + 1][j][k];
        }
        x_component /= (h.x * h.x); 

        double y_component = -2 * grid[i][j][k];
        y_component += grid[i][SubstitutePeriodicIndex(j - 1)][k];
        y_component += grid[i][SubstitutePeriodicIndex(j + 1)][k];
        y_component /= (h.y * h.y);

        double z_component = -2 * grid[i][j][k];
        z_component += grid[i][j][SubstitutePeriodicIndex(k - 1)];
        z_component += grid[i][j][SubstitutePeriodicIndex(k + 1)];
        z_component /= (h.z * h.z);
        return x_component + y_component + z_component;
    }

    void FillFirstLayer(const double a2, const Uanalytical& u) {
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < N + 1; ++i) {
            for (int j = 0; j < N + 1; ++j) {
                for (int k = 0; k < N + 1; ++k) {
                    double delta_uh = GetLaplasOperator(grid_n0, i, j, k);
                    grid_n1[i][j][k] = grid_n0[i][j][k] + 0.5 * a2 * tau * tau * delta_uh;
                }
            }
        }
    }

    void FillNextLayer(const Uanalytical& u, const double a2) {
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < N + 1; ++i) {
            for (int j = 0; j < N + 1; ++j) {
                for (int k = 0; k < N + 1; ++k) {
                    double delta_uh = GetLaplasOperator(grid_n1, i, j, k);
                    grid_n2[i][j][k] = a2 * tau * tau * delta_uh - grid_n0[i][j][k] + 2 * grid_n1[i][j][k];
                }
            }
        }
    }

    double CountMaxError(const Uanalytical& u, int step) const {
        double max_error = 0;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                for (int k = 0; k < N; ++k) {
                    double true_val = u.Value(i * h.x, j * h.y, k * h.z, tau * step);
                    double error = std::fabs(true_val - grid_n2[i][j][k]);
                    if (error > max_error) {
                        max_error = error;
                    }
                }
            }
        }
        return max_error;
    }

    void DumpResult(const std::string& folder, int step) const {
        std::ofstream out(folder + "-solution_step_" + std::to_string(step) + ".json");
        out << std::setprecision(15);
        out << "{\n    \"values\": [\n";
        for (int i = 0; i < N + 1; ++i) {
            out << "        [\n";
            for (int j = 0; j < N + 1; ++j) {
                out << "            [";
                for (int k = 0; k < N + 1; ++k) {
                    out << grid_n2[i][j][k];
                    if (k != N) {
                        out << ", ";
                    }
                }
                out << "]";
                if (j != N) {
                    out << ",";
                }
                out << "\n";
            }
            out << "       ]";
            if (i != N) {
                out << ",";
            }
            out << "\n";
        }
        out << "   ]\n}";
        out.close();
    }

    void MoveGridLayer() {
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < N + 1; ++i) {
            for (int j = 0; j < N + 1; ++j) {
                for (int k = 0; k < N + 1; ++k) {
                    grid_n0[i][j][k] = grid_n1[i][j][k];
                    grid_n1[i][j][k] = grid_n2[i][j][k];
                }
            }
        }
    }

    bool operator==(const Uapproximated& c) const {
        for (int i = 0; i < N + 1; ++i) {
            for (int j = 0; j < N + 1; ++j) {
                for (int k = 0; k < N + 1; ++k) {
                    if (std::fabs(grid_n0[i][j][k] - c.grid_n0[i][j][k]) >= 1e-8) {
                        return false;
                    }
                    if (std::fabs(grid_n1[i][j][k] - c.grid_n1[i][j][k]) >= 1e-8) {
                        return false;
                    }
                    if (std::fabs(grid_n2[i][j][k] - c.grid_n2[i][j][k]) >= 1e-8) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
};

Point GetBoarderValues(int argc, char** argv) {
    Point p = Point{
        .x = std::stod(argv[1]),
        .y = std::stod(argv[2]),
        .z = std::stod(argv[3])
    };
    if (p.x <= 0 || p.y <= 0 || p.z <= 0) {
        std::logic_error("Non-positive Lx, Ly or Lz!");
    }
    return p;
}

double GetTimeDiff(const std::chrono::system_clock::time_point& t) {
    return (std::chrono::high_resolution_clock::now() - t).count();
}

void DumpMetaInfo(const std::string& path_prefix, const std::vector<double>& errors, double memory_time, double diffur_time, int Nthreads) {
    std::ofstream out(path_prefix + "_meta.json");
    out << std::setprecision(15);
    out << "{\n";
    out << "    \"errors\": [";
    for (int i = 0; i < errors.size(); ++i) {
        out << errors[i];
        if (i != errors.size() - 1) {
            out << ", ";
        }
    }
    out << "],\n";
    out << "    \"memory_time_nanosec\": " << memory_time << ",\n";
    out << "    \"diffur_time_nanosec\": " << diffur_time << ",\n";
    out << "    \"nthreads\": " << Nthreads << "\n";
    out << "}";
    out.close();
}

int main(int argc, char** argv) {
    // Parse console arguments
    if (argc < 8) {
        std::cout << "Bad initial arguements! The list is following:\n";
        std::cout << "    1) Lx (double) - Ox axis limit;\n";
        std::cout << "    2) Ly (double) - Oy axis limit;\n";
        std::cout << "    3) Lz (double) - Oz axis limit;\n";
        std::cout << "    4) T (double) - timeline axis limit;\n";
        std::cout << "    5) N (int) - grid size for space axis;\n";
        std::cout << "    6) K (int) - grid size for time axis;\n";
        std::cout << "    7) Nthreads (int) - number of threads.\n";
        std::cout << "Example: `./a.out 1 1 1 1 10 10 10`.\n";
        throw std::logic_error("Not enough initial values");
    }
    Point L = GetBoarderValues(argc, argv);
    const double T = std::stod(argv[4]);
    const int N = std::stoi(argv[5]);
    const int K = std::stoi(argv[6]);
    const int Nthreads = std::stoi(argv[7]);

    const std::string kDumpPrefixFolder{
        std::string("results/") + argv[1] + "_" + argv[2] + "_" + argv[3] + "_" + argv[4] + "_" + argv[5] + "_" + argv[6] + "_" + argv[7]
    };

    // Setup number of threads
    std::cout << Nthreads << std::endl;
    omp_set_num_threads(Nthreads);

    // Setup initial conditions
    const double a2 = 4 / (9 / (L.x * L.x) + 4 / (L.y * L.y) + 4 / (L.z * L.z));
    Uanalytical u(L);

    // Init timers and approximate solution
    double diffur_time = 0;
    double memory_time = 0;
    std::chrono::system_clock::time_point t_diffur;
    std::chrono::system_clock::time_point t_memory;
    std::vector<double> errors{};
    t_memory = std::chrono::high_resolution_clock::now();
    Uapproximated solution(N, K, L, T);
    memory_time += GetTimeDiff(t_memory);
    // Count first two steps of calculation
    t_diffur = std::chrono::high_resolution_clock::now();
    solution.FillZerothLayer(u);
    solution.FillFirstLayer(a2, u);
    diffur_time += GetTimeDiff(t_diffur);
    // Count the rest steps
    for (int t = 2; t < K; ++t) {
        t_diffur = std::chrono::high_resolution_clock::now();
        solution.FillNextLayer(u, a2);
        diffur_time += GetTimeDiff(t_diffur);
        double max_error = solution.CountMaxError(u, t);
        errors.push_back(max_error);
        // if (solution.N < 200 && t % 9 == 0) // Not to dump too many big files!
        //     solution.DumpResult(kDumpPrefixFolder, t);
        t_memory = std::chrono::high_resolution_clock::now();
        solution.MoveGridLayer();
        memory_time += GetTimeDiff(t_memory);
    }
    std::cout << "Worst: " << errors[errors.size() - 1] << std::endl;
    DumpMetaInfo(kDumpPrefixFolder, errors, memory_time, diffur_time, omp_get_num_threads());
    return 0;
}