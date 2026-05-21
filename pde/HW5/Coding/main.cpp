// 对流：二维数组 u[i][j]，首行 n m；波动：u_tt - u_xx = f，首行 nt nx。

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "PDEwave.h"

// ---------- 对流问题 f（u_t + u_x = f）----------
static double f_advection(double x, double t) {
    return 2.0 * t * std::sin(x) + t * t * std::cos(x);
}

// ---------- 波动问题 f(x,t) = (t^2 - x^2) sin(x t)，0<x,t<1 ----------
static double f_wave(double x, double t) {
    return (t * t - x * x) * std::sin(x * t);
}

static void write_grid(std::ostream& os, int n, int m,
                       const std::vector<std::vector<double>>& u) {
    os << n << ' ' << m << '\n';
    os << std::scientific << std::setprecision(16);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            if (j) {
                os << ' ';
            }
            os << u[i][j];
        }
        os << '\n';
    }
}

static void apply_scheme(int scheme, std::vector<std::vector<double>>& u, int n, int m,
                         double h, double g, double a) {
    switch (scheme) {
        case 1:
            pdewave::solve_ut_plus_aux(u, n, m, h, g, a, f_advection, 0.0, 0.0);
            break;
        case 2:
            pdewave::solve_ut_plus_aux_lax_wendroff(u, n, m, h, g, a, f_advection, 0.0, 0.0);
            break;
        case 3:
            pdewave::solve_ut_plus_aux_backward_euler(u, n, m, h, g, a, f_advection, 0.0, 0.0);
            break;
        default:
            std::cerr << "scheme must be 1 (upwind), 2 (Lax–Wendroff), or 3 (backward Euler)\n";
            std::exit(1);
    }
}

static void print_usage(const char* prog) {
    std::cerr << "Usage:\n"
              << "  Advection: " << prog << " n m [scheme outfile]\n"
              << "    scheme: 1=upwind, 2=Lax–Wendroff, 3=backward Euler; outfile => no stdout dump\n"
              << "    Explicit (1,2) need CFL<=1: with dt=5/n,dx=10/m,a=1 ⇒ n>=ceil(m/2).\n"
              << "    Or no args: read n m from stdin; scheme=1; write u_data.txt + stdout\n"
              << "  Wave:      " << prog << " wave nt nx [method outfile]\n"
              << "    u_tt-u_xx=(t^2-x^2)sin(xt), 0<t,x<1; u(x,0)=0, u_t(x,0)=x;\n"
              << "    u(0,t)=0, u(1,t)=sin(t). Grid: x_j=j/(nx-1), t_i=i/(nt-1).\n"
              << "    method: 1=explicit central (leapfrog), 2=theta-implicit (theta=0.5),\n"
              << "            3=compact spatial + leapfrog (PDEwave.h wave1d::*).\n"
              << "    outfile default wave_u.txt; if only one extra arg and no digit, treated as outfile.\n";
}

static int run_advection(int argc, char** argv) {
    int n = 0;
    int m = 0;
    int scheme = 1;
    std::string out_name = "u_data.txt";
    bool write_stdout = true;

    if (argc >= 5) {
        n = std::atoi(argv[1]);
        m = std::atoi(argv[2]);
        scheme = std::atoi(argv[3]);
        out_name = argv[4];
        write_stdout = false;
        if (n <= 0 || m <= 0) {
            print_usage(argv[0]);
            return 1;
        }
    } else if (argc == 1) {
        if (!(std::cin >> n >> m) || n <= 0 || m <= 0) {
            std::cerr << "Expected two positive integers: n m\n";
            return 1;
        }
    } else {
        print_usage(argv[0]);
        return 1;
    }

    const double g = 10.0 / static_cast<double>(m);
    const double h = 5.0 / static_cast<double>(n);
    constexpr double a = 1.0;

    const double cfl = pdewave::cfl_number(a, h, g);
    pdewave::warn_if_unstable_cfl(cfl, std::cerr, n, m);

    std::vector<std::vector<double>> u(n, std::vector<double>(m));
    apply_scheme(scheme, u, n, m, h, g, a);

    if (write_stdout) {
        write_grid(std::cout, n, m, u);
    }

    std::ofstream fout(out_name);
    if (!fout) {
        std::cerr << "Failed to open file for writing: " << out_name << '\n';
        return 1;
    }
    write_grid(fout, n, m, u);
    std::cerr << "Written: " << out_name << " (n=" << n << ", m=" << m << ", scheme=" << scheme
              << ")\n";
    return 0;
}

static int run_wave(int argc, char** argv) {
    if (argc < 4) {
        print_usage(argv[0]);
        return 1;
    }
    const int nt = std::atoi(argv[2]);
    const int nx = std::atoi(argv[3]);
    int method = 1;
    std::string out_name = "wave_u.txt";

    if (argc >= 6) {
        method = std::atoi(argv[4]);
        out_name = argv[5];
    } else if (argc == 5) {
        const char* a4 = argv[4];
        if (a4[0] >= '1' && a4[0] <= '3' && a4[1] == '\0') {
            method = std::atoi(a4);
        } else {
            out_name = a4;
        }
    }

    if (nt < 2 || nx < 2) {
        std::cerr << "wave: need nt>=2, nx>=2\n";
        return 1;
    }
    if (method < 1 || method > 3) {
        std::cerr << "wave: method must be 1 (central), 2 (theta=0.5), or 3 (compact)\n";
        return 1;
    }

    const double dx = 1.0 / static_cast<double>(nx - 1);
    const double dt = 1.0 / static_cast<double>(nt - 1);

    if (method == 1 || method == 3) {
        pdewave::wave1d::warn_if_unstable_wave_cfl(pdewave::wave1d::wave_cfl(dt, dx));
    }

    std::function<double(double)> bc_left = [](double /*t*/) { return 0.0; };
    std::function<double(double)> bc_right = [](double t) { return std::sin(t); };
    std::function<double(double)> v0 = [](double x) { return x; };

    std::vector<std::vector<double>> u(nt, std::vector<double>(nx, 0.0));

    for (int j = 0; j < nx; ++j) {
        u[0][j] = 0.0;
    }
    u[0][0] = bc_left(0.0);
    u[0][nx - 1] = bc_right(0.0);

    pdewave::wave1d::init_second_time_level(u, nx, dt, dx, f_wave, v0, bc_left, bc_right);

    for (int i = 0; i < nt; ++i) {
        const double t = static_cast<double>(i) * dt;
        u[i][0] = bc_left(t);
        u[i][nx - 1] = bc_right(t);
    }

    constexpr double k_theta = 0.5;
    switch (method) {
        case 1:
            pdewave::wave1d::march_central_explicit(u, nt, nx, dt, dx, f_wave, bc_left, bc_right);
            break;
        case 2:
            pdewave::wave1d::march_theta_implicit(u, nt, nx, dt, dx, k_theta, f_wave, bc_left,
                                                  bc_right);
            break;
        case 3:
            pdewave::wave1d::march_compact_leapfrog(u, nt, nx, dt, dx, f_wave, bc_left, bc_right);
            break;
        default:
            break;
    }

    std::ofstream fout(out_name);
    if (!fout) {
        std::cerr << "Failed to open file for writing: " << out_name << '\n';
        return 1;
    }
    write_grid(fout, nt, nx, u);
    std::cerr << "Written: " << out_name << " (wave method=" << method << ", nt=" << nt
              << ", nx=" << nx << ", dt=" << dt << ", dx=" << dx << ")\n";
    return 0;
}

int main(int argc, char** argv) {
    if (argc >= 2 && std::string(argv[1]) == "wave") {
        return run_wave(argc, argv);
    }
    return run_advection(argc, argv);
}
