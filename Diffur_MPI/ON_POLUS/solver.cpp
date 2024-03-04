#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <string>
#include <stdexcept>
#include <vector>
#include <fstream>

#include "omp.h"
#include "mpi.h"

const int kBufTagXy = 1;
const int kBufTagXz = 2;
const int kBufTagYz = 3;

struct Point {
    double x;
    double y;
    double z;
};

struct GridDivision{
    int by_x;
    int by_y;
    int by_z;
};

std::ostream &operator<<(std::ostream &os, const GridDivision &m) { 
    return os << "{" << m.by_x << ", " << m.by_y << ", " << m.by_z << "}";
}

std::vector<int> GetDivisors(int x) {
    std::vector<int> divisors;
    while (!(x % 2)) {
        x >>= 1;
        divisors.push_back(2);
    }
    int upper_bound = std::sqrt(x);
    for (int i = 3; i <= upper_bound; i += 2) {
        while (!(x % i)) {
            x /= i;
            divisors.push_back(i);
        }
    }
    divisors.push_back(x);
    return divisors;
}

GridDivision GetDivision(int x) {
    std::vector<int> divisors = GetDivisors(x);
    GridDivision g {
        .by_x = 1,
        .by_y = 1,
        .by_z = 1
    };
    bool reverse = false;
    for (int i = 0; i < divisors.size(); i += 3) {
        if (!reverse) {
            g.by_x *= divisors[i];
            if (i + 2 < divisors.size()) {
                g.by_z *= divisors[i + 2];
            }
        } else {
            g.by_z *= divisors[i];
            if (i + 2 < divisors.size()) {
                g.by_x *= divisors[i + 2];
            }
        }
        if (i + 1 < divisors.size()) {
            g.by_y *= divisors[i + 1];
        }
    }
    return g;
}

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
    GridDivision g, proc_pos, sizes;

    
    double tau;
    int ntasks, task_id;
    int N, K;
    double*** grid_n0;
    double*** grid_n1;
    double*** grid_n2;
    double *xy_boarder_u = nullptr, *xy_boarder_d = nullptr;
    double *xz_boarder_u = nullptr, *xz_boarder_d = nullptr;
    double *yz_boarder_u = nullptr, *yz_boarder_d = nullptr;
    double *buf_xy = nullptr, *buf_xz = nullptr, *buf_yz = nullptr;

    explicit Uapproximated(int N, int K, const Point& L, double T, int ntasks, int task_id): N(N), K(K), ntasks(ntasks), task_id(task_id) {
        h = Point{
            .x = L.x / N,
            .y = L.y / N,
            .z = L.z / N
        };
        g = GetDivision(ntasks);
        if (!task_id) {
            std::cout << "Grid division: " << g << std::endl;
        }
        proc_pos.by_z = task_id % g.by_z;
        proc_pos.by_y = (task_id / g.by_z) % g.by_y;
        proc_pos.by_x = (task_id / g.by_y / g.by_z) % g.by_x;
        sizes.by_x = N / g.by_x + (N % g.by_x > proc_pos.by_x);
        sizes.by_y = N / g.by_y + (N % g.by_y > proc_pos.by_y);
        sizes.by_z = N / g.by_z + (N % g.by_z > proc_pos.by_z);
        tau = T / K;
        grid_n0 = new double**[sizes.by_x];
        grid_n1 = new double**[sizes.by_x];
        grid_n2 = new double**[sizes.by_x];
        for (int i = 0; i < sizes.by_x; ++i) {
            grid_n0[i] = new double*[sizes.by_y];
            grid_n1[i] = new double*[sizes.by_y];
            grid_n2[i] = new double*[sizes.by_y];
            for (int j = 0; j < sizes.by_y; ++j) {
                grid_n0[i][j] = new double[sizes.by_z];
                grid_n1[i][j] = new double[sizes.by_z];
                grid_n2[i][j] = new double[sizes.by_z];
            }
        }
        buf_xy = new double[sizes.by_x * sizes.by_y];
        buf_xz = new double[sizes.by_x * sizes.by_z];
        buf_yz = new double[sizes.by_y * sizes.by_z];
        xy_boarder_u = new double[sizes.by_x * sizes.by_y];
        xy_boarder_d = new double[sizes.by_x * sizes.by_y];
        xz_boarder_u = new double[sizes.by_x * sizes.by_z];
        xz_boarder_d = new double[sizes.by_x * sizes.by_z];
        yz_boarder_u = new double[sizes.by_y * sizes.by_z];
        yz_boarder_d = new double[sizes.by_y * sizes.by_z];
    }

    ~Uapproximated() {
        for (int i = 0; i < sizes.by_x; ++i) {
            for (int j = 0; j < sizes.by_y; ++j) {
                delete[] grid_n0[i][j]; 
                delete[] grid_n1[i][j]; 
                delete[] grid_n2[i][j];
            }
            delete[] grid_n0[i]; 
            delete[] grid_n1[i]; 
            delete[] grid_n2[i];
        }
        delete[] grid_n0; 
        delete[] grid_n1; 
        delete[] grid_n2;
        delete[] xz_boarder_d;
        delete[] xz_boarder_u;
        delete[] xy_boarder_d;
        delete[] xy_boarder_u;
        delete[] buf_xy;
        delete[] buf_xz;
        delete[] buf_yz;
        delete[] yz_boarder_d;
        delete[] yz_boarder_u;
    }

    void SendRecvLowerBuffers(double*** grid) {
        // Copy buffers
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < sizes.by_x; ++i) {
            for (int j = 0; j < sizes.by_y; ++j) {
                buf_xy[i * sizes.by_y + j] = grid[i][j][0];
            }
        }
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < sizes.by_x; ++i) {
            for (int j = 0; j < sizes.by_z; ++j) {
                buf_xz[i * sizes.by_z + j] = grid[i][0][j];
            }
        }
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < sizes.by_y; ++i) {
            for (int j = 0; j < sizes.by_z; ++j) {
                buf_yz[i * sizes.by_z + j] = grid[0][i][j];
            }
        }

        // ISend buffers
        MPI_Request send_request_z, send_request_y, send_request_x;
        int dest;
        int receiving_task_id_z = proc_pos.by_z - 1;
        if (receiving_task_id_z < 0) receiving_task_id_z = g.by_z - 1;
        dest = receiving_task_id_z + proc_pos.by_y * g.by_z + proc_pos.by_x * g.by_y * g.by_z;
        MPI_Isend(buf_xy, sizes.by_x * sizes.by_y, MPI_DOUBLE, dest, kBufTagXy, MPI_COMM_WORLD, &send_request_z);
        int receiving_task_id_y = proc_pos.by_y - 1;
        if (receiving_task_id_y < 0) receiving_task_id_y = g.by_y - 1;
        dest = proc_pos.by_z + receiving_task_id_y * g.by_z + proc_pos.by_x * g.by_y * g.by_z;
        MPI_Isend(buf_xz, sizes.by_x * sizes.by_z, MPI_DOUBLE, dest, kBufTagXz, MPI_COMM_WORLD, &send_request_y);
        int receiving_task_id_x = proc_pos.by_x - 1;
        if (receiving_task_id_x < 0) receiving_task_id_x = g.by_x - 1;
        dest = proc_pos.by_z + proc_pos.by_y * g.by_z + receiving_task_id_x * g.by_y * g.by_z;
        MPI_Isend(buf_yz, sizes.by_y * sizes.by_z, MPI_DOUBLE, dest, kBufTagYz, MPI_COMM_WORLD, &send_request_x);

        // Recv buffers
        MPI_Status request_z, request_y, request_x;
        int src;
        int sender_task_id_z = proc_pos.by_z + 1;
        if (sender_task_id_z >= g.by_z) sender_task_id_z = 0;
        src = sender_task_id_z + proc_pos.by_y * g.by_z + proc_pos.by_x * g.by_y * g.by_z;
        MPI_Recv(xy_boarder_u, sizes.by_x * sizes.by_y, MPI_DOUBLE, src, kBufTagXy, MPI_COMM_WORLD, &request_z);
        int sender_task_id_y = proc_pos.by_y + 1;
        if (sender_task_id_y >= g.by_y) sender_task_id_y = 0;
        src = proc_pos.by_z + sender_task_id_y * g.by_z + proc_pos.by_x * g.by_y * g.by_z;
        MPI_Recv(xz_boarder_u, sizes.by_x * sizes.by_z, MPI_DOUBLE, src, kBufTagXz, MPI_COMM_WORLD, &request_y);
        int sender_task_id_x = proc_pos.by_x + 1;
        if (sender_task_id_x >= g.by_x) sender_task_id_x = 0;
        src = proc_pos.by_z + proc_pos.by_y * g.by_z + sender_task_id_x * g.by_y * g.by_z;
        MPI_Recv(yz_boarder_u, sizes.by_z * sizes.by_y, MPI_DOUBLE, src, kBufTagYz, MPI_COMM_WORLD, &request_x);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    void SendRecvUpperBuffers(double*** grid) {
        // Copy buffers
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < sizes.by_x; ++i) {
            for (int j = 0; j < sizes.by_y; ++j) {
                buf_xy[i * sizes.by_y + j] = grid[i][j][sizes.by_z - 1];
            }
        }
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < sizes.by_x; ++i) {
            for (int j = 0; j < sizes.by_z; ++j) {
                buf_xz[i* sizes.by_z + j] = grid[i][sizes.by_y - 1][j];
            }
        }
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < sizes.by_y; ++i) {
            for (int j = 0; j < sizes.by_z; ++j) {
                buf_yz[i * sizes.by_z + j] = grid[sizes.by_x - 1][i][j];
            }
        }

        // ISend buffers
        MPI_Request send_request_z, send_request_y, send_request_x;
        int dest;
        int receiving_task_id_z = proc_pos.by_z + 1;
        if (receiving_task_id_z >= g.by_z) receiving_task_id_z = 0;
        dest = receiving_task_id_z + proc_pos.by_y * g.by_z + proc_pos.by_x * g.by_y * g.by_z;
        MPI_Isend(buf_xy, sizes.by_x * sizes.by_y, MPI_DOUBLE, dest, kBufTagXy, MPI_COMM_WORLD, &send_request_z);
        int receiving_task_id_y = proc_pos.by_y + 1;
        if (receiving_task_id_y >= g.by_y) receiving_task_id_y = 0;
        dest = proc_pos.by_z + receiving_task_id_y * g.by_z + proc_pos.by_x * g.by_y * g.by_z;
        MPI_Isend(buf_xz, sizes.by_x * sizes.by_z, MPI_DOUBLE, dest, kBufTagXz, MPI_COMM_WORLD, &send_request_y);
        int receiving_task_id_x = proc_pos.by_x + 1;
        if (receiving_task_id_x >= g.by_x) receiving_task_id_x = 0;
        dest = proc_pos.by_z + proc_pos.by_y * g.by_z + receiving_task_id_x * g.by_y * g.by_z;
        MPI_Isend(buf_yz, sizes.by_y * sizes.by_z, MPI_DOUBLE, dest, kBufTagYz, MPI_COMM_WORLD, &send_request_x);

        // Recv buffers
        MPI_Status request_z, request_y, request_x;
        int src;
        int sender_task_id_z = proc_pos.by_z - 1;
        if (sender_task_id_z < 0) sender_task_id_z = g.by_z - 1;
        src = sender_task_id_z + proc_pos.by_y * g.by_z + proc_pos.by_x * g.by_y * g.by_z;
        MPI_Recv(xy_boarder_d, sizes.by_x * sizes.by_y, MPI_DOUBLE, src, kBufTagXy, MPI_COMM_WORLD, &request_z);
        int sender_task_id_y = proc_pos.by_y - 1;
        if (sender_task_id_y < 0) sender_task_id_y = g.by_y - 1;
        src = proc_pos.by_z + sender_task_id_y * g.by_z + proc_pos.by_x * g.by_y * g.by_z;
        MPI_Recv(xz_boarder_d, sizes.by_x * sizes.by_z, MPI_DOUBLE, src, kBufTagXz, MPI_COMM_WORLD, &request_y);
        int sender_task_id_x = proc_pos.by_x - 1;
        if (sender_task_id_x < 0) sender_task_id_x = g.by_x - 1;
        src = proc_pos.by_z + proc_pos.by_y * g.by_z + sender_task_id_x * g.by_y * g.by_z;
        MPI_Recv(yz_boarder_d, sizes.by_z * sizes.by_y, MPI_DOUBLE, src, kBufTagYz, MPI_COMM_WORLD, &request_x);
        MPI_Barrier(MPI_COMM_WORLD);
    }

    void SendBuffers(bool is_grid_n0) {
        double ***grid = grid_n1;
        if (is_grid_n0) {
            grid = grid_n0;
        }
        SendRecvLowerBuffers(grid);
        SendRecvUpperBuffers(grid);
    }

    int SubstitutePeriodicIndex(int index) const {
        if (index < 0) return N - 1;
        if (index >= N) return 0; 
        return index;
    }

    double GetLaplasOperator(double*** grid, int i, int j, int k) const {
        double x_component = -2 * grid[i][j][k];
        if (i > 0) {
            x_component += grid[i - 1][j][k];
        }
        if (i < sizes.by_x - 1) {
            x_component += grid[i + 1][j][k];
        }
        if (i == 0 && proc_pos.by_x > 0) {
            x_component += yz_boarder_d[j * sizes.by_z + k];
        }
        if (i == sizes.by_x - 1 && proc_pos.by_x < g.by_x - 1) {
            x_component += yz_boarder_u[j * sizes.by_z + k];
        }
        x_component /= (h.x * h.x); 

        double y_component = -2 * grid[i][j][k];
        if (j > 0) {
            y_component += grid[i][j - 1][k];
        }
        if (j < sizes.by_y - 1) {
            y_component += grid[i][j + 1][k];
        }
        if (j == 0) {
            y_component += xz_boarder_d[i * sizes.by_z + k];
        }
        if (j == sizes.by_y - 1) {
            y_component += xz_boarder_u[i * sizes.by_z + k];
        }
        y_component /= (h.y * h.y);

        double z_component = -2 * grid[i][j][k];
        if (k > 0) {
            z_component += grid[i][j][k - 1];
        }
        if (k < sizes.by_z - 1) {
            z_component += grid[i][j][k + 1];
        }
        if (k == 0) {
            z_component += xy_boarder_d[i * sizes.by_y + j];
        }
        if (j == sizes.by_y - 1) {
            z_component += xy_boarder_u[i * sizes.by_y + j];
        }
        z_component /= (h.z * h.z);
        return x_component + y_component + z_component;
    }

    void FillZerothLayer(const Uanalytical& u) {
        int x_offset = sizes.by_x * proc_pos.by_x;
        int y_offset = sizes.by_y * proc_pos.by_y;
        int z_offset = sizes.by_z * proc_pos.by_z;
        #pragma omp parallel for collapse(3)
        for (int i = x_offset; i < x_offset + sizes.by_x; ++i) {
            for (int j = y_offset; j < y_offset + sizes.by_y; ++j) {
                for (int k = z_offset; k < z_offset + sizes.by_z; ++k) {
                    grid_n0[i - x_offset][j - y_offset][k - z_offset] = u.Value(h.x * i, h.y * j, h.z * k, 0);
                }
            }
        }
    }

    void FillFirstLayer(const double a2, const Uanalytical& u) {
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < sizes.by_x; ++i) {
            for (int j = 0; j < sizes.by_y; ++j) {
                for (int k = 0; k < sizes.by_z; ++k) {
                    double delta_uh = GetLaplasOperator(grid_n0, i, j, k);
                    grid_n1[i][j][k] = grid_n0[i][j][k] + 0.5 * a2 * tau * tau * delta_uh;
                }
            }
        }
    }

    void FillNextLayer(const Uanalytical& u, const double a2) {
        #pragma omp parallel for collapse(3)
        for (int i = 0; i < sizes.by_x; ++i) {
            for (int j = 0; j < sizes.by_y; ++j) {
                for (int k = 0; k < sizes.by_z; ++k) {
                    double delta_uh = GetLaplasOperator(grid_n1, i, j, k);
                    grid_n2[i][j][k] = a2 * tau * tau * delta_uh - grid_n0[i][j][k] + 2 * grid_n1[i][j][k];
                }
            }
        }
    }

    double CountMaxError(const Uanalytical& u, int step) const {
        int x_offset = sizes.by_x * proc_pos.by_x;
        int y_offset = sizes.by_y * proc_pos.by_y;
        int z_offset = sizes.by_z * proc_pos.by_z;
        double max_error = 0;
        #pragma omp parallel for collapse(3)
        for (int i = x_offset; i < x_offset + sizes.by_x; ++i) {
            for (int j = y_offset; j < y_offset + sizes.by_y; ++j) {
                for (int k = z_offset; k < z_offset + sizes.by_z; ++k) {
                    double true_val = u.Value(i * h.x, j * h.y, k * h.z, tau * step);
                    double error = std::fabs(true_val - grid_n2[i - x_offset][j - y_offset][k - z_offset]);
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
        double*** tmp = grid_n0;
        grid_n0 = grid_n1;
        grid_n1 = grid_n2;
        grid_n2 = tmp;
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

    void PrintMetaElement() const {
        int i = 6;
        int j = 27;
        int k = 35;
        std::cout << "Element pos: (" << i + proc_pos.by_x * sizes.by_x 
                  << ", " << j + proc_pos.by_y * sizes.by_y << ", " 
                  << k + proc_pos.by_z * sizes.by_z << "). "
                  << "Value: " << grid_n2[i][j][k] << std::endl;
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

void DumpMetaInfo(const std::string& path_prefix, const double error, double memory_time, double diffur_time, int Nthreads) {
    std::ofstream out(path_prefix + "_meta.json");
    out << std::setprecision(15);
    out << "{\n";
    out << "    \"errors\": [";
    out << error;
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

    int task_id, ntasks;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    // Setup number of threads
    if (!task_id) {
        std::cout << "Number of threads: " << Nthreads << std::endl;
        std::cout << "Number of processes: " << ntasks << std::endl;
    }
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
    Uapproximated solution(N, K, L, T, ntasks, task_id);
    memory_time += GetTimeDiff(t_memory);
    // Count first two steps of calculation
    t_diffur = std::chrono::high_resolution_clock::now();
    solution.FillZerothLayer(u);
    solution.SendBuffers(true);
    solution.FillFirstLayer(a2, u);
    diffur_time += GetTimeDiff(t_diffur);
    // Count the rest steps
    for (int t = 2; t < K; ++t) {
        t_diffur = std::chrono::high_resolution_clock::now();
        solution.SendBuffers(true);
        solution.FillNextLayer(u, a2);
        MPI_Barrier(MPI_COMM_WORLD);
        diffur_time += GetTimeDiff(t_diffur);
        double max_error = solution.CountMaxError(u, t);
        errors.push_back(max_error);
        t_memory = std::chrono::high_resolution_clock::now();
        solution.MoveGridLayer();
        memory_time += GetTimeDiff(t_memory);
    }
    double max_global_error = errors[errors.size() - 1];
    MPI_Allreduce(&errors[errors.size() - 1], &max_global_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    if (!task_id) {
        std::cout << "Worst case: " << max_global_error << std::endl;
        DumpMetaInfo(kDumpPrefixFolder, max_global_error, memory_time, diffur_time, omp_get_num_threads());
    }
    MPI_Finalize();
    return 0;
}