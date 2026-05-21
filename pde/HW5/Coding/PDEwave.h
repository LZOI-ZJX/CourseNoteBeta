#pragma once

#include <iostream>
#include <vector>

namespace pdewave {

// ---------------------------------------------------------------------------
// 公共：u_t + a u_x = rhs(x, t)，a > 0；网格 t_i = i*dt, x_j = j*dx，
// u[i][j] ≈ u(t_i, x_j)，i = 0..n-1，j = 0..m-1。
// ---------------------------------------------------------------------------

inline double cfl_number(double a, double dt, double dx) { return a * dt / dx; }

// 对 u_t + a u_x = … 的显式迎风 / LW：稳定性需 CFL = a*dt/dx <= 1。
// 若 main 中取 dt=5/n、dx=10/m（t∈[0,5], x∈[0,10]），则 CFL = a*m/(2n)；a=1 时需 n >= ceil(m/2)。
inline void warn_if_unstable_cfl(double cfl, std::ostream& err = std::cerr, int n = 0, int m = 0) {
    if (cfl > 1.0 + 1e-12) {
        err << "Warning: CFL = a*dt/dx = " << cfl
            << " > 1; explicit upwind / Lax-Wendroff are unstable; solution blows up (not a plot bug).\n";
        if (n > 0 && m > 0) {
            const int n_min = (m + 1) / 2;
            err << "  For dt=5/n, dx=10/m, a=1: need n >= " << n_min
                << " (with m=" << m << "), or reduce m, or use backward Euler (scheme 3).\n";
        }
    }
}

// ---------------------------------------------------------------------------
// 显式迎风格式（一阶）
// (u[i,j]-u[i-1,j])/dt + a*(u[i-1,j]-u[i-1,j-1])/dx = rhs(x_j, t_i)
// ---------------------------------------------------------------------------
template <typename Rhs>
void solve_ut_plus_aux(std::vector<std::vector<double>>& u, int n, int m, double dt,
                       double dx, double a, const Rhs& rhs, double u_ic = 0.0,
                       double u_bc = 0.0) {
    const double lambda = cfl_number(a, dt, dx);
    for (int j = 0; j < m; ++j) {
        u[0][j] = u_ic;
    }
    for (int i = 0; i < n; ++i) {
        u[i][0] = u_bc;
    }
    for (int i = 1; i < n; ++i) {
        for (int j = 1; j < m; ++j) {
            const double x = static_cast<double>(j) * dx;
            const double t = static_cast<double>(i) * dt;
            u[i][j] = u[i - 1][j] - lambda * (u[i - 1][j] - u[i - 1][j - 1]) + dt * rhs(x, t);
        }
    }
}

// ---------------------------------------------------------------------------
// Lax–Wendroff（显式、二阶）：对 u_t + a u_x = rhs，
// u^{n+1}_j = u^n_j - (ν/2)(u^n_{j+1}-u^n_{j-1}) + (ν²/2)(u^n_{j+1}-2u^n_j+u^n_{j-1})
//           + dt * rhs(x_j, t_{n+1})，ν = a*dt/dx。
// 内点 j = 1..m-2；j = m-1 处用一阶迎风；m < 3 时整行用迎风。
// ---------------------------------------------------------------------------
template <typename Rhs>
void solve_ut_plus_aux_lax_wendroff(std::vector<std::vector<double>>& u, int n, int m,
                                    double dt, double dx, double a, const Rhs& rhs,
                                    double u_ic = 0.0, double u_bc = 0.0) {
    const double nu = cfl_number(a, dt, dx);
    for (int j = 0; j < m; ++j) {
        u[0][j] = u_ic;
    }
    for (int i = 0; i < n; ++i) {
        u[i][0] = u_bc;
    }
    for (int i = 1; i < n; ++i) {
        const double t = static_cast<double>(i) * dt;
        if (m < 3) {
            for (int j = 1; j < m; ++j) {
                const double x = static_cast<double>(j) * dx;
                u[i][j] = u[i - 1][j] - nu * (u[i - 1][j] - u[i - 1][j - 1]) + dt * rhs(x, t);
            }
            continue;
        }
        for (int j = 1; j < m - 1; ++j) {
            const double x = static_cast<double>(j) * dx;
            const double um = u[i - 1][j - 1];
            const double u0 = u[i - 1][j];
            const double up = u[i - 1][j + 1];
            u[i][j] = u0 - 0.5 * nu * (up - um) + 0.5 * nu * nu * (up - 2.0 * u0 + um) + dt * rhs(x, t);
        }
        {
            const int j = m - 1;
            const double x = static_cast<double>(j) * dx;
            u[i][j] = u[i - 1][j] - nu * (u[i - 1][j] - u[i - 1][j - 1]) + dt * rhs(x, t);
        }
    }
}

// ---------------------------------------------------------------------------
// Backward Euler（隐式）+ 空间迎风：在每一时间层对 x 前代求解
// (u[i,j]-u[i-1,j])/dt + a*(u[i,j]-u[i,j-1])/dx = rhs(x_j, t_i)
// => (1+λ) u[i,j] - λ u[i,j-1] = u[i-1,j] + dt*rhs, λ = a*dt/dx。
// 隐式迎风通常对纯对流无条件稳定（仍受离散误差限制）。
// ---------------------------------------------------------------------------
template <typename Rhs>
void solve_ut_plus_aux_backward_euler(std::vector<std::vector<double>>& u, int n, int m,
                                      double dt, double dx, double a, const Rhs& rhs,
                                      double u_ic = 0.0, double u_bc = 0.0) {
    const double lambda = cfl_number(a, dt, dx);
    const double denom = 1.0 + lambda;
    for (int j = 0; j < m; ++j) {
        u[0][j] = u_ic;
    }
    for (int i = 0; i < n; ++i) {
        u[i][0] = u_bc;
    }
    for (int i = 1; i < n; ++i) {
        for (int j = 1; j < m; ++j) {
            const double x = static_cast<double>(j) * dx;
            const double t = static_cast<double>(i) * dt;
            const double b = u[i - 1][j] + dt * rhs(x, t);
            u[i][j] = (b + lambda * u[i][j - 1]) / denom;
        }
    }
}

// ===========================================================================
// 一维波动方程 u_tt - u_xx = F(x,t)（可推广系数，此处单位波速 c=1）
// 存储与对流部分一致：u[i][j] ≈ u(t_i, x_j)，t_i = i*dt，x_j = j*dx。
// 与手写下标 u_j^n 对应关系：u_j^n ≡ u[n][j]。式中 g≡dt，h≡dx。
// ===========================================================================
namespace wave1d {

inline double wave_cfl(double dt, double dx) { return dt / dx; }

inline void warn_if_unstable_wave_cfl(double cfl, std::ostream& err = std::cerr) {
    if (cfl > 1.0 + 1e-12) {
        err << "Warning: wave CFL = dt/dx = " << cfl
            << " > 1; explicit central / compact leapfrog may be unstable.\n";
    }
}

namespace detail {

inline void thomas_solve_inplace(const std::vector<double>& a_sub,
                                 const std::vector<double>& b_main,
                                 const std::vector<double>& c_sup,
                                 const std::vector<double>& d_rhs,
                                 std::vector<double>& x_out) {
    const int n = static_cast<int>(d_rhs.size());
    std::vector<double> b_mod = b_main;
    std::vector<double> d_mod = d_rhs;
    for (int i = 1; i < n; ++i) {
        const double w = a_sub[i] / b_mod[i - 1];
        b_mod[i] -= w * c_sup[i - 1];
        d_mod[i] -= w * d_mod[i - 1];
    }
    x_out.assign(n, 0.0);
    x_out[n - 1] = d_mod[n - 1] / b_mod[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x_out[i] = (d_mod[i] - c_sup[i] * x_out[i + 1]) / b_mod[i];
    }
}

}  // namespace detail

// ---------------------------------------------------------------------------
// 中心差分（显式 leapfrog）：单步推进 n-1,n → n+1，写 u[n+1][*] 内点。
// (u^{n+1}_j - 2u^n_j + u^{n-1}_j)/dt^2 - (u^n_{j+1}-2u^n_j+u^n_{j-1})/dx^2 = F(j*dx, n*dt)
// ---------------------------------------------------------------------------
template <typename F>
void step_central_explicit(const std::vector<std::vector<double>>& u, int n, int nx, double dt,
                           double dx, const F& force, std::vector<double>& u_np1_row) {
    const double r2 = (dt / dx) * (dt / dx);
    u_np1_row.resize(nx);
    for (int j = 0; j < nx; ++j) {
        u_np1_row[j] = u[n][j];
    }
    for (int j = 1; j < nx - 1; ++j) {
        const double x = static_cast<double>(j) * dx;
        const double t = static_cast<double>(n) * dt;
        u_np1_row[j] = 2.0 * u[n][j] - u[n - 1][j] +
                       r2 * (u[n][j + 1] - 2.0 * u[n][j] + u[n][j - 1]) + dt * dt * force(x, t);
    }
}

// 完整时间迭代：需已给定 u[0][*]、u[1][*]；每步后由 bc_* 写 u[n+1][0]、u[n+1][nx-1]。
template <typename F, typename Bc>
void march_central_explicit(std::vector<std::vector<double>>& u, int nt, int nx, double dt,
                            double dx, const F& force, const Bc& bc_left, const Bc& bc_right) {
    std::vector<double> buf(nx);
    for (int n = 1; n < nt - 1; ++n) {
        step_central_explicit(u, n, nx, dt, dx, force, buf);
        const double tp1 = static_cast<double>(n + 1) * dt;
        u[n + 1][0] = bc_left(tp1);
        u[n + 1][nx - 1] = bc_right(tp1);
        for (int j = 1; j < nx - 1; ++j) {
            u[n + 1][j] = buf[j];
        }
    }
}

// 由 u(0)、初速 v0(x) 生成第二层 u[1]（二阶）：u^1 = u^0 + dt*v0 + (dt^2/2)(u_xx^0 + F^0)。
template <typename F, typename V0, typename Bc>
void init_second_time_level(std::vector<std::vector<double>>& u, int nx, double dt, double dx,
                            const F& force, const V0& v0, const Bc& bc_left, const Bc& bc_right) {
    for (int j = 1; j < nx - 1; ++j) {
        const double x = static_cast<double>(j) * dx;
        const double dxx_u0 = u[0][j + 1] - 2.0 * u[0][j] + u[0][j - 1];
        u[1][j] = u[0][j] + dt * v0(x) + 0.5 * dt * dt * (dxx_u0 / (dx * dx) + force(x, 0.0));
    }
    u[1][0] = bc_left(dt);
    u[1][nx - 1] = bc_right(dt);
}

// ---------------------------------------------------------------------------
// 加权隐式（θ）：δ_t^2 u_j^n - [θ δ_x^2 u^{n+1} + (1-2θ) δ_x^2 u^n + θ δ_x^2 u^{n-1}] = F_j^n
// 每层对内点 j=1..nx-2 解三对角。bc 在时刻 t 的左右 Dirichlet 值。
// ---------------------------------------------------------------------------
template <typename F, typename Bc>
void step_theta_implicit(std::vector<std::vector<double>>& u, int n, int nx, double dt, double dx,
                         double theta, const F& force, const Bc& bc_left, const Bc& bc_right) {
    const double r2 = (dt / dx) * (dt / dx);
    const double lam = theta * r2;
    const int ni = nx - 2;
    if (ni <= 0) {
        return;
    }
    const double t_np1 = static_cast<double>(n + 1) * dt;
    const double t_n = static_cast<double>(n) * dt;
    const double u_left = bc_left(t_np1);
    const double u_right = bc_right(t_np1);

    std::vector<double> a_sub(ni, -lam);
    a_sub[0] = 0.0;
    std::vector<double> b_main(ni, 1.0 + 2.0 * lam);
    std::vector<double> c_sup(ni, -lam);
    c_sup[ni - 1] = 0.0;
    std::vector<double> d_rhs(ni);

    for (int k = 0; k < ni; ++k) {
        const int j = k + 1;
        const double x = static_cast<double>(j) * dx;
        const double dxx_n = u[n][j + 1] - 2.0 * u[n][j] + u[n][j - 1];
        const double dxx_nm1 = u[n - 1][j + 1] - 2.0 * u[n - 1][j] + u[n - 1][j - 1];
        d_rhs[k] = 2.0 * u[n][j] - u[n - 1][j] + (1.0 - 2.0 * theta) * r2 * dxx_n +
                   theta * r2 * dxx_nm1 + dt * dt * force(x, t_n);
    }
    d_rhs[0] += lam * u_left;
    d_rhs[ni - 1] += lam * u_right;

    std::vector<double> x_sol;
    detail::thomas_solve_inplace(a_sub, b_main, c_sup, d_rhs, x_sol);

    u[n + 1][0] = u_left;
    u[n + 1][nx - 1] = u_right;
    for (int k = 0; k < ni; ++k) {
        u[n + 1][k + 1] = x_sol[k];
    }
}

template <typename F, typename Bc>
void march_theta_implicit(std::vector<std::vector<double>>& u, int nt, int nx, double dt, double dx,
                          double theta, const F& force, const Bc& bc_left, const Bc& bc_right) {
    for (int n = 1; n < nt - 1; ++n) {
        step_theta_implicit(u, n, nx, dt, dx, theta, force, bc_left, bc_right);
    }
}

// ---------------------------------------------------------------------------
// 紧致格式（空间）：Padé 型 δ_x^2
// (1/12) u''_{j-1} + (5/6) u''_j + (1/12) u''_{j+1} = (u_{j+1}-2u_j+u_{j-1})/dx^2
// 再显式 leapfrog：u^{n+1}_j = 2u^n_j - u^{n-1}_j + dt^2 (u_xx)_j^n + dt^2 F^n_j
// ---------------------------------------------------------------------------
inline void solve_compact_laplacian_row(const std::vector<double>& u_row, int nx, double dx,
                                        std::vector<double>& uxx_out) {
    uxx_out.assign(nx, 0.0);
    if (nx < 3) {
        return;
    }
    uxx_out[0] = (u_row[2] - 2.0 * u_row[1] + u_row[0]) / (dx * dx);
    uxx_out[nx - 1] =
        (u_row[nx - 1] - 2.0 * u_row[nx - 2] + u_row[nx - 3]) / (dx * dx);

    const int ni = nx - 2;
    std::vector<double> a_sub(ni, 1.0 / 12.0);
    a_sub[0] = 0.0;
    std::vector<double> b_main(ni, 5.0 / 6.0);
    std::vector<double> c_sup(ni, 1.0 / 12.0);
    c_sup[ni - 1] = 0.0;
    std::vector<double> d_rhs(ni);
    for (int k = 0; k < ni; ++k) {
        const int j = k + 1;
        d_rhs[k] = (u_row[j + 1] - 2.0 * u_row[j] + u_row[j - 1]) / (dx * dx);
    }
    d_rhs[0] -= (1.0 / 12.0) * uxx_out[0];
    d_rhs[ni - 1] -= (1.0 / 12.0) * uxx_out[nx - 1];

    std::vector<double> x_mid;
    detail::thomas_solve_inplace(a_sub, b_main, c_sup, d_rhs, x_mid);
    for (int k = 0; k < ni; ++k) {
        uxx_out[k + 1] = x_mid[k];
    }
}

template <typename F>
void step_compact_leapfrog(const std::vector<std::vector<double>>& u, int n, int nx, double dt,
                           double dx, const F& force, std::vector<double>& u_np1_row,
                           std::vector<double>& work_uxx) {
    u_np1_row.resize(nx);
    solve_compact_laplacian_row(u[n], nx, dx, work_uxx);
    for (int j = 0; j < nx; ++j) {
        u_np1_row[j] = u[n][j];
    }
    for (int j = 1; j < nx - 1; ++j) {
        const double x = static_cast<double>(j) * dx;
        const double t = static_cast<double>(n) * dt;
        u_np1_row[j] = 2.0 * u[n][j] - u[n - 1][j] + dt * dt * work_uxx[j] + dt * dt * force(x, t);
    }
}

template <typename F, typename Bc>
void march_compact_leapfrog(std::vector<std::vector<double>>& u, int nt, int nx, double dt,
                            double dx, const F& force, const Bc& bc_left, const Bc& bc_right) {
    std::vector<double> buf(nx);
    std::vector<double> uxx(nx);
    for (int n = 1; n < nt - 1; ++n) {
        step_compact_leapfrog(u, n, nx, dt, dx, force, buf, uxx);
        const double tp1 = static_cast<double>(n + 1) * dt;
        u[n + 1][0] = bc_left(tp1);
        u[n + 1][nx - 1] = bc_right(tp1);
        for (int j = 1; j < nx - 1; ++j) {
            u[n + 1][j] = buf[j];
        }
    }
}

}  // namespace wave1d

}  // namespace pdewave
