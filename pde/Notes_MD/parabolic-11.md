# 第五节 多维问题

## 二维问题

初边值问题：

$$
\begin{cases}
\displaystyle \frac{\partial u}{\partial t} = a \left(\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}\right), & 0 < x, y < 1 \\[0.6em]
u(x, y, 0) = f(x, y), & 0 < x, y < 1 \\[0.6em]
u(0, y, t) = u(1, y, t) = u(x, 0, t) = u(x, 1, t) = 0, & t \ge 0
\end{cases}
$$

取 $\Delta t = \tau$，$\Delta x = \Delta y = h$，形成网格：

$$
D_h = \left\{ (x_j, y_l, t_n) \,\middle|\, \begin{aligned} x_j &= jh,\quad j = 0, 1, \dots, J,\quad Jh = 1 \\ y_l &= lh,\quad l = 0, 1, \dots, J,\quad Jh = 1 \\ t_n &= n\tau,\quad n \ge 0 \end{aligned} \right\}
$$

并记：$\delta_x^2 u_{jl}^n = u_{j+1, l}^n - 2u_{jl}^n + u_{j-1, l}^n$，

$$
\delta_y^2 u_{jl}^n = u_{j, l+1}^n - 2u_{jl}^n + u_{j, l-1}^n
$$

## 5.1 一维格式的直接推广

一维扩散方程古典显式格式推广：

$$
\frac{u_{jl}^{n+1} - u_{jl}^n}{\tau} = \frac{a}{h^2}\left(\delta_x^2 u_{jl}^n + \delta_y^2 u_{jl}^n\right) ， \tag{5.1}
$$

利用 Taylor 展开式分析截断误差为 $O(\tau + h^2)$。

计算格式：

$$
u_{jl}^{n+1} = u_{jl}^n + a\frac{\tau}{h^2}\left(\delta_x^2 u_{jl}^n + \delta_y^2 u_{jl}^n\right) ， \tag{5.2}
$$

可用 Fourier 法分析稳定性：

设 $\lambda = \tau / h^2$，$u_{jl}^n = v^n e^{\mathrm{i}k_1 jh} e^{\mathrm{i}k_2 lh}$，代入 (5.2)，有

$$
v^{n+1} = \bigl(1 + 2a\lambda(\cos k_1 h - 1) + 2a\lambda(\cos k_2 h - 1)\bigr) v^n
$$

传播因子：$G(\lambda, k_1, k_2) = 1 - 4a\lambda\left(\sin^2\frac{k_1 h}{2} + \sin^2\frac{k_2 h}{2}\right)$。

解得当 $a\lambda \le \frac{1}{4}$ 时，$|G| \le 1$，格式稳定。

一般 $p$ 维显式格式稳定当且仅当 $a\lambda \le \frac{1}{2p}$。

一维古典隐式格式的推广：

$$
\frac{u_{jl}^{n+1} - u_{jl}^n}{\tau} = \frac{a}{h^2}\left(\delta_x^2 u_{jl}^{n+1} + \delta_y^2 u_{jl}^{n+1}\right) ， \tag{5.3}
$$

格式的截断误差仍为 $O(\tau + h^2)$，但增长因子：

$$
G(\lambda, k_1, k_2) = \frac{1}{1 + 4a\lambda\left(\sin^2\frac{k_1 h}{2} + \sin^2\frac{k_2 h}{2}\right)}
$$

$|G(\lambda, k_1, k_2)| \le 1$，格式无条件稳定。

Crank--Nicolson 格式：

$$
\frac{u_{jl}^{n+1} - u_{jl}^n}{\tau} = \frac{a}{2h^2}\left[\delta_x^2\left(u_{jl}^{n+1} + u_{jl}^n\right) + \delta_y^2\left(u_{jl}^{n+1} + u_{jl}^n\right)\right] ， \tag{5.4}
$$

格式具有二阶精度，$E = O(\tau^2 + h^2)$。

$$
G(\lambda, k_1, k_2) = \frac{1 - 2a\lambda\sin^2\frac{k_1 h}{2} - 2a\lambda\sin^2\frac{k_2 h}{2}}{1 + 2a\lambda\sin^2\frac{k_1 h}{2} + 2a\lambda\sin^2\frac{k_2 h}{2}}
$$

对任意 $\lambda, k_{1,2}$，有 $|G(\lambda, k_1, k_2)| \le 1$，格式绝对稳定。

## 5.2 交替方向隐式格式

### PR（Peaceman-Rachford）格式

$$
\begin{cases}
\displaystyle \frac{u_{jl}^{n+\frac{1}{2}} - u_{jl}^n}{\tau/2} = \frac{a}{h^2}\left(\delta_x^2 u_{jl}^{n+\frac{1}{2}} + \delta_y^2 u_{jl}^n\right) \\[0.8em]
\displaystyle \frac{u_{jl}^{n+1} - u_{jl}^{n+\frac{1}{2}}}{\tau/2} = \frac{a}{h^2}\left(\delta_x^2 u_{jl}^{n+\frac{1}{2}} + \delta_y^2 u_{jl}^{n+1}\right) ，
\end{cases} \tag{5.5}
$$

注：求解时每层需用两次追赶法。

### PR 格式的特点

#### 截断误差

消去 $u_{jl}^{n+\frac{1}{2}}$、$\delta_x^2 u_{jl}^{n+\frac{1}{2}}$，可得

$$
\begin{aligned}
&\left(1 + \frac{1}{4}\frac{\tau^2 a^2}{h^4}\delta_x^2\delta_y^2\right)\frac{u_{jl}^{n+1} - u_{jl}^n}{\tau} \\
&= \frac{a}{h^2}\left(\delta_x^2 + \delta_y^2\right)\frac{u_{jl}^{n+1} + u_{jl}^n}{2} ，
\end{aligned} \tag{5.6}
$$

Taylor 展开分析截断误差：$O(\tau^2 + h^2)$。

所以 PR 格式是二阶精度的格式。

#### 稳定性

改写 (5.6) 为

$$
\begin{aligned}
&\left(1 - \frac{a\lambda}{2}\delta_x^2\right)\left(1 - \frac{a\lambda}{2}\delta_y^2\right) u_{jl}^{n+1} \\
&= \left(1 + \frac{a\lambda}{2}\delta_x^2\right)\left(1 + \frac{a\lambda}{2}\delta_y^2\right) u_{jl}^n ，
\end{aligned} \tag{5.7}
$$

Fourier 方法分析稳定性，得增长因子

$$
G(\lambda, k_1, k_2) = \frac{\bigl(1 - 2a\lambda\sin^2\frac{k_1 h}{2}\bigr)\bigl(1 - 2a\lambda\sin^2\frac{k_2 h}{2}\bigr)}{\bigl(1 + 2a\lambda\sin^2\frac{k_1 h}{2}\bigr)\bigl(1 + 2a\lambda\sin^2\frac{k_2 h}{2}\bigr)}
$$

所以 PR 格式是无条件稳定的。

### 附：PR 格式的另一种推导方法

也可以通过 Crank--Nicolson 格式变形的途径得到 Peaceman-Rachford 格式。首先将 Crank--Nicolson 格式变形为

$$
u_{jl}^{n+1} - \frac{a\lambda}{2}\left[\delta_x^2 u_{jl}^{n+1} + \delta_y^2 u_{jl}^{n+1}\right] = u_{jl}^n + \frac{a\lambda}{2}\left[\delta_x^2 u_{jl}^n + \delta_y^2 u_{jl}^n\right]
$$

$$
\begin{aligned}
&u_{jl}^{n+1} - \frac{a\lambda}{2}\left[\delta_x^2 u_{jl}^{n+1} + \delta_y^2 u_{jl}^{n+1}\right] + \frac{a^2\lambda^2}{4}\delta_x^2 \delta_y^2 u_{jl}^{n+1} \\
&= u_{jl}^n + \frac{a\lambda}{2}\left[\delta_x^2 u_{jl}^n + \delta_y^2 u_{jl}^n\right] + \frac{a^2\lambda^2}{4}\delta_x^2 \delta_y^2\left[u_{jl}^n + O(\tau)\right]
\end{aligned}
$$

（注：因 $\lambda=\tau/h^2$，$\frac{a^2\lambda^2}{4}\delta_x^2\delta_y^2$ 作用于 $O(\tau)$ 的增量上整体为 $O(\tau^3)$；与式子两边同除以 $\tau$ 的尺度相配时，对截断误差贡献 $O(\tau^2)$，从而与 PR 格式的二阶截断相容。）

忽略高阶小项，并因式分解得到

$$
\left(1 - \frac{a\lambda}{2}\delta_x^2\right)\left(1 - \frac{a\lambda}{2}\delta_y^2\right) u_{jl}^{n+1} = \left(1 + \frac{a\lambda}{2}\delta_x^2\right)\left(1 + \frac{a\lambda}{2}\delta_y^2\right) u_{jl}^n 。
$$

令 $\left(1 - \frac{a\lambda}{2}\delta_x^2\right)u_{jl}^{n+\frac{1}{2}} = \left(1 + \frac{a\lambda}{2}\delta_y^2\right)u_{jl}^n$，

则 $\left(1 - \frac{a\lambda}{2}\delta_y^2\right)u_{jl}^{n+1} = \left(1 + \frac{a\lambda}{2}\delta_x^2\right)u_{jl}^{n+\frac{1}{2}}$。

所以有
$$
\begin{cases}
\left(1 - \frac{a\lambda}{2}\delta_x^2\right) u_{jl}^{n+\frac{1}{2}} = \left(1 + \frac{a\lambda}{2}\delta_y^2\right) u_{jl}^n \\
\left(1 - \frac{a\lambda}{2}\delta_y^2\right) u_{jl}^{n+1} = \left(1 + \frac{a\lambda}{2}\delta_x^2\right) u_{jl}^{n+\frac{1}{2}}
\end{cases}
$$

注：PR 格式与 CN 格式有相同的精度，但计算量小很多。

### 其他形式的 ADI

$$
\begin{cases}
\left(1 - \frac{a\lambda}{2}\delta_x^2\right) u_{jl}^{n+\frac{1}{2}} = \left(1 + \frac{a\lambda}{2}\delta_x^2\right)\left(1 + \frac{a\lambda}{2}\delta_y^2\right) u_{jl}^n \\[0.6em]
\left(1 - \frac{a\lambda}{2}\delta_y^2\right) u_{jl}^{n+1} = u_{jl}^{n+\frac{1}{2}}
\end{cases}
$$

格式特点：1．二阶精度；2．无条件稳定。

### Douglas 格式

Peaceman-Rachford 格式的缺点是无法推广到三维：它不能如二维格式那样保持对称形式的传播因子，无条件稳定性不再成立。下面的格式可以推广到三维。

改写 (5.6) 为

$$
\begin{aligned}
&\left(1 - \frac{a\lambda}{2}\delta_x^2\right)\left(1 - \frac{a\lambda}{2}\delta_y^2\right)\frac{u_{jl}^{n+1} - u_{jl}^n}{\tau} \\
&= \frac{a}{h^2}\left(\delta_x^2 + \delta_y^2\right) u_{jl}^n ，
\end{aligned} \tag{5.8}
$$

把 (5.8) 分裂为

$$
\begin{cases}
\left(1 - \frac{a\lambda}{2}\delta_x^2\right)\dfrac{u_{jl}^{n+\frac{1}{2}} - u_{jl}^n}{\tau} = \dfrac{a}{h^2}\left(\delta_x^2 + \delta_y^2\right) u_{jl}^n \\[0.9em]
\left(1 - \frac{a\lambda}{2}\delta_y^2\right)\dfrac{u_{jl}^{n+1} - u_{jl}^n}{\tau} = \dfrac{u_{jl}^{n+\frac{1}{2}} - u_{jl}^n}{\tau}
\end{cases}
$$

还可以表示为

$$
\begin{cases}
\left(1 - \frac{a\lambda}{2}\delta_x^2\right)\left(u_{jl}^{n+\frac{1}{2}} - u_{jl}^n\right) = a\lambda\left(\delta_x^2 + \delta_y^2\right) u_{jl}^n \\[0.6em]
u_{jl}^{n+1} - u_{jl}^{n+\frac{1}{2}} = \frac{a\lambda}{2}\delta_y^2\left(u_{jl}^{n+1} - u_{jl}^n\right)
\end{cases}
$$

此格式称为 Douglas 格式；其稳定性、局部截断误差与 PR 格式相同，且易于向高维推广。

### 三维 Douglas 格式

（注：在三维格式中，记空间网格节点为 $(x_i, y_j, z_l)$，对应解记为 $u_{ijl}^n$。）

$$
\begin{cases}
\left(1 - \frac{a\lambda}{2}\delta_x^2\right)\left(u_{ijl}^{n+\frac{1}{3}} - u_{ijl}^n\right) = a\lambda\left(\delta_x^2 + \delta_y^2 + \delta_z^2\right) u_{ijl}^n \\[0.7em]
u_{ijl}^{n+\frac{2}{3}} - u_{ijl}^{n+\frac{1}{3}} = \dfrac{a\lambda}{2}\delta_y^2\left(u_{ijl}^{n+\frac{2}{3}} - u_{ijl}^n\right) \\[0.7em]
u_{ijl}^{n+1} - u_{ijl}^{n+\frac{2}{3}} = \dfrac{a\lambda}{2}\delta_z^2\left(u_{ijl}^{n+1} - u_{ijl}^n\right)
\end{cases}
$$

同书 (5.35)。

## 5.3 局部一维格式

对方程建立如下的局部一维格式：

$$
\begin{cases}
\displaystyle \frac{u_{jl}^{n+\frac{1}{2}} - u_{jl}^n}{\tau} = \frac{a}{h^2}\delta_x^2\left(\frac{u_{jl}^{n+\frac{1}{2}} + u_{jl}^n}{2}\right) \\[0.9em]
\displaystyle \frac{u_{jl}^{n+1} - u_{jl}^{n+\frac{1}{2}}}{\tau} = \frac{a}{h^2}\delta_y^2\left(\frac{u_{jl}^{n+1} + u_{jl}^{n+\frac{1}{2}}}{2}\right)
\end{cases}
$$

计算时，可化为 (5.9)：

$$
\begin{cases}
\left(1 - \frac{a\lambda}{2}\delta_x^2\right) u_{jl}^{n+\frac{1}{2}} = \left(1 + \frac{a\lambda}{2}\delta_x^2\right) u_{jl}^n \\[0.7em]
\left(1 - \frac{a\lambda}{2}\delta_y^2\right) u_{jl}^{n+1} = \left(1 + \frac{a\lambda}{2}\delta_y^2\right) u_{jl}^{n+\frac{1}{2}} ，
\end{cases} \tag{5.9}
$$

用两次追赶法即可求得 $\{u_{jl}^{n+1}\}$。

此格式与 PR 格式等价，是二阶绝对稳定格式；且可推广至三维。

## 5.4 预测--校正格式

预测--校正联立格式 (5.10) 为

$$
\begin{cases}
\displaystyle \frac{u_{jl}^{n+\frac{1}{4}} - u_{jl}^n}{\tau/2} = \frac{a}{h^2}\delta_x^2 u_{jl}^{n+\frac{1}{4}} \\[0.8em]
\displaystyle \frac{u_{jl}^{n+\frac{1}{2}} - u_{jl}^{n+\frac{1}{4}}}{\tau/2} = \frac{a}{h^2}\delta_y^2 u_{jl}^{n+\frac{1}{2}} \\[0.8em]
\displaystyle \frac{u_{jl}^{n+1} - u_{jl}^n}{\tau} = \frac{a}{h^2}\left(\delta_x^2 + \delta_y^2\right) u_{jl}^{n+\frac{1}{2}} ，
\end{cases} \tag{5.10}
$$

消去 $u_{jl}^{n+\frac{1}{4}}$ 与 $u_{jl}^{n+\frac{1}{2}}$，可得 (5.8)。**所以此格式与 Douglas 格式等价**，是二阶绝对稳定格式。
