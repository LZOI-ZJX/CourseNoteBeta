# 第五章 椭圆型方程的差分方法

## 主要内容

- Poisson 方程
- 差分方程的性质
- 边界条件的处理
- 变系数方程
- 双调和方程
- 特征值问题

## 第一节 Poisson 方程

### 控制方程与离散区域

考虑方程

$$
-\Delta u = -\left(\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}\right) = f(x,y), \quad (x,y) \in D \tag{1}
$$

其中 $D$ 是平面内有界区域，边界用 $\partial D$ 表示。

不妨设 $D$ 为矩形区域，$D = \{(x,y) \mid 0 < x < 1,\ 0 < y < 1\}$，方程满足边值条件（Dirichlet 条件）

$$
u(x,y) = a(x,y), \quad (x,y) \in \partial D.
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/c763622a-db61-41ee-8d7f-a09a4a0bdbf5/12f560103c8120f483193eb256d0bcac24b45e0ec3e9c138a1cdff86f4bb5c4d.jpg)

#### 区域 $D$ 的内点

$$
D_h = \{(x_i, y_j) \mid x_i = i h_1,\ 1 \le i \le I;\ y_j = j h_2,\ 1 \le j \le J\}.
$$

#### 区域 $D$ 的边界点

$$
\partial D_h = \{(x_i, y_j) \mid x_i = i h_1,\ y_j = j h_2,\ i = 0,\ I+1,\ 0 \le j \le J+1;\ j = 0,\ J+1,\ 0 \le i \le I+1\}.
$$

### 1.1 五点差分格式

现在假定 $(x_i, y_j)$ 为内点。记**离散 Laplace 算子**（不含负号）

$$
\Delta_h u_{ij} = \frac{u_{i+1,j} - 2u_{ij} + u_{i-1,j}}{h_1^2} + \frac{u_{i,j+1} - 2u_{ij} + u_{i,j-1}}{h_2^2},
$$

则用二阶中心差商逼近 $-\Delta u$ 后，差分方程为

$$
-\Delta_h u_{ij} = f_{ij} \tag{1.1}
$$

差分格式为：

$$
\begin{cases}
-\Delta_h u_{ij} = f_{ij}, & (x_i, y_j) \in D_h \\
u_{ij} = a_{ij}, & (x_i, y_j) \in \partial D_h
\end{cases}
$$

利用 Taylor 展式

$$
\begin{aligned}
\frac{u(x_{i+1}, y_j) - 2u(x_i, y_j) + u(x_{i-1}, y_j)}{h_1^2} &= \frac{\partial^2 u(x_i, y_j)}{\partial x^2} + \frac{h_1^2}{12}\frac{\partial^4 u(x_i, y_j)}{\partial x^4} + \frac{h_1^4}{360}\frac{\partial^6 u(x_i, y_j)}{\partial x^6} + O(h_1^6), \\
\frac{u(x_i, y_{j+1}) - 2u(x_i, y_j) + u(x_i, y_{j-1})}{h_2^2} &= \frac{\partial^2 u(x_i, y_j)}{\partial y^2} + \frac{h_2^2}{12}\frac{\partial^4 u(x_i, y_j)}{\partial y^4} + \frac{h_2^4}{360}\frac{\partial^6 u(x_i, y_j)}{\partial y^6} + O(h_2^6),
\end{aligned}
$$

可得逼近 $-\Delta u$ 的截断误差（与 (1.1) 中左端一致）

$$
\begin{aligned}
R_{ij}(u) &= -\Delta_h u(x_i, y_j) - (-\Delta u(x_i, y_j)) \\
&= -\frac{1}{12}\left[ h_1^2 \frac{\partial^4 u(x_i, y_j)}{\partial x^4} + h_2^2 \frac{\partial^4 u(x_i, y_j)}{\partial y^4} \right] + O(h^4) \\
&= O(h_1^2 + h_2^2)\quad\text{（若记 } h = \max(h_1, h_2)\text{，则亦为 } O(h^2)\text{）}.
\end{aligned}
$$

由于差分方程 (1.1) 中只出现 $u$ 在 $(i,j)$ 及其四个邻点上的值，故称为五点差分格式。

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/c763622a-db61-41ee-8d7f-a09a4a0bdbf5/d4ccc17fc4461b6bc596059c0e18176e8a028f3c020394927ed3c1ca347f77cf.jpg)

特别取正方形网格：$h_1 = h_2 = h$，则差分方程 (1.1) 简化为

$$
4u_{ij} - (u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1}) = h^2 f_{ij}.
$$

若 $f \equiv 0$（Laplace 方程），则有

$$
4u_{ij} = u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1}.
$$

单位正方形内部结点上的 $J^2$ 个线性方程写成矩阵形式为

$$
AU = K,
$$

其中 $A$ 是 $J^2$ 阶方阵：

$$
A = \left[\begin{array}{ccccc}
B & -I & & & \\
-I & B & -I & & \\
 & \ddots & \ddots & \ddots & \\
 & & -I & B & -I \\
 & & & -I & B
\end{array}\right].
$$

$I$ 是 $J$ 阶单位方阵；$B$ 是 $J$ 阶方阵：

$$
B = \left[\begin{array}{ccccc}
4 & -1 & & & \\
-1 & 4 & -1 & & \\
 & \ddots & \ddots & \ddots & \\
 & & -1 & 4 & -1 \\
 & & & -1 & 4
\end{array}\right].
$$

### 1.2 九点差分格式

**注** 若将 1.1 节中 Taylor 展式关于 $x$、$y$ 的两式相加，则得

$$
\begin{aligned}
\Delta_h u(x_i, y_j) &= \frac{u(x_{i+1}, y_j) - 2u(x_i, y_j) + u(x_{i-1}, y_j)}{h_1^2} + \frac{u(x_i, y_{j+1}) - 2u(x_i, y_j) + u(x_i, y_{j-1})}{h_2^2} \\
&= \Delta u(x_i, y_j) + \frac{1}{12}\left[ h_1^2 \frac{\partial^4 u(x_i, y_j)}{\partial x^4} + h_2^2 \frac{\partial^4 u(x_i, y_j)}{\partial y^4} \right] + O(h^4) \\
&= \Delta u(x_i, y_j) + \frac{1}{12}\left(h_1^2 \frac{\partial^2}{\partial x^2} + h_2^2 \frac{\partial^2}{\partial y^2}\right)\left(\frac{\partial^2 u(x_i, y_j)}{\partial x^2} + \frac{\partial^2 u(x_i, y_j)}{\partial y^2}\right) \\
&\quad - \frac{h_1^2 + h_2^2}{12}\frac{\partial^4 u(x_i, y_j)}{\partial x^2 \partial y^2} + O(h^4)
\end{aligned}
$$

注意到 $-\Delta u(x_i, y_j) = f(x_i, y_j)$，即 $\Delta u(x_i, y_j) = -f(x_i, y_j)$，代入上式得：

$$
\Delta_h u(x_i, y_j) = -f(x_i, y_j) - \frac{1}{12}\left(h_1^2 \frac{\partial^2 f(x_i, y_j)}{\partial x^2} + h_2^2 \frac{\partial^2 f(x_i, y_j)}{\partial y^2}\right) - \frac{h_1^2 + h_2^2}{12}\frac{\partial^4 u(x_i, y_j)}{\partial x^2 \partial y^2} + O(h^4).
$$

又

$$
\begin{aligned}
\frac{\partial^4 u(x_i, y_j)}{\partial x^2 \partial y^2} &= \frac{u_{xx}''(x_i, y_{j+1}) - 2u_{xx}''(x_i, y_j) + u_{xx}''(x_i, y_{j-1})}{h_2^2} + O(h_2^2) \\
&= \frac{1}{h_1^2 h_2^2} \Bigl[ u(x_{i+1}, y_{j+1}) - 2u(x_i, y_{j+1}) + u(x_{i-1}, y_{j+1}) \\
&\quad - 2\bigl(u(x_{i+1}, y_j) - 2u(x_i, y_j) + u(x_{i-1}, y_j)\bigr) \\
&\quad + u(x_{i+1}, y_{j-1}) - 2u(x_i, y_{j-1}) + u(x_{i-1}, y_{j-1}) \Bigr] + O(h^2).
\end{aligned}
$$

因此

$$
\begin{aligned}
&\Delta_h u(x_i, y_j) + \frac{1}{12}\frac{h_1^2 + h_2^2}{h_1^2 h_2^2} \Bigl[ 4u(x_i, y_j) - 2\bigl(u(x_{i-1}, y_j) + u(x_i, y_{j-1}) + u(x_{i+1}, y_j) + u(x_i, y_{j+1})\bigr) \\
&\quad + u(x_{i-1}, y_{j-1}) + u(x_{i+1}, y_{j-1}) + u(x_{i+1}, y_{j+1}) + u(x_{i-1}, y_{j+1}) \Bigr] \\
&= -f(x_i, y_j) - \frac{1}{12}\left(h_1^2 \frac{\partial^2 f(x_i, y_j)}{\partial x^2} + h_2^2 \frac{\partial^2 f(x_i, y_j)}{\partial y^2}\right) + O(h^4).
\end{aligned}
$$

舍去截断误差项，便得到逼近方程的九点差分格式：

$$
\begin{aligned}
&-\Delta_h u_{ij} - \frac{1}{12}\frac{h_1^2 + h_2^2}{h_1^2 h_2^2} \Bigl[ 4u_{ij} - 2(u_{i-1,j} + u_{i,j-1} + u_{i+1,j} + u_{i,j+1}) \\
&\quad + u_{i-1,j-1} + u_{i+1,j-1} + u_{i+1,j+1} + u_{i-1,j+1} \Bigr] \\
&= f_{ij} + \frac{1}{12}\left(h_1^2 \frac{\partial^2 f(x_i, y_j)}{\partial x^2} + h_2^2 \frac{\partial^2 f(x_i, y_j)}{\partial y^2}\right).
\end{aligned}
$$

其截断误差的阶为 $O(h^4)$。

#### 附：九点差分格式推导（方法二）

取 $\xi = h\dfrac{\partial}{\partial x}$，$\eta = h\dfrac{\partial}{\partial y}$，$D = \dfrac{\partial^2}{\partial x \partial y}$，

$$
\xi^2 + \eta^2 = h^2 \Delta, \quad \xi\eta = h^2 D, \quad \xi^4 + \eta^4 = h^4(\Delta^2 - 2D^2).
$$

由 Taylor 公式有：

$$
u(x+h) = \left(1 + h\frac{d}{dx} + \frac{h^2}{2}\frac{d^2}{dx^2} + \cdots + \frac{h^n}{n!}\frac{d^n}{dx^n}\right) u(x) = e^{h\frac{d}{dx}} u(x).
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/c763622a-db61-41ee-8d7f-a09a4a0bdbf5/42a2892245a6651eaf7bbb133dfe3eb492d76c62e5fd70ae1dbf26dca3c0e15b.jpg)

令 $S_1 = u_1 + u_2 + u_3 + u_4$，$S_2 = u_5 + u_6 + u_7 + u_8$，

$$
\begin{aligned}
S_1 &= u_1 + u_2 + u_3 + u_4 \\
&= \left(e^\xi + e^\eta + e^{-\xi} + e^{-\eta}\right) u_0 \\
&= 4u_0 + h^2 \Delta u_0 + \frac{1}{12} h^4(\Delta^2 - 2D^2) u_0 + O(h^6).
\end{aligned}
$$

$$
S_2 = 4u_0 + 2h^2 \Delta u_0 + \frac{1}{6} h^4(\Delta^2 + 4D^2) u_0 + O(h^6).
$$

从上面两式消去 $D^2 u_0$，并忽略高阶项，得

$$
4S_1 + S_2 - 20u_0 = 6h^2 \Delta u_0 + \frac{1}{2} h^4 \Delta^2 u_0.
$$

注意到 $-\Delta u_0 = f_{ij}$，且 $\Delta^2 u_0 = -\Delta f_{ij}$，代入上式得：

$$
20u_0 - 4S_1 - S_2 = 6h^2 f_{ij} + \frac{1}{2} h^4 \Delta f_{ij}.
$$

#### 积分法

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/c763622a-db61-41ee-8d7f-a09a4a0bdbf5/ca1d0bb6344df58c5c56cb30472c75b62945f8af5ebf29b84dd922863b80c450.jpg)

$$
-\int_D \Delta u \,\mathrm{d}x\,\mathrm{d}y = \int_D f \,\mathrm{d}x\,\mathrm{d}y.
$$

由 Gauss 公式，$\displaystyle \int_D \Delta u \,\mathrm{d}x\,\mathrm{d}y = \int_{\partial D} \frac{\partial u}{\partial n} \,\mathrm{d}s$，故

$$
-\int_{\partial D} \frac{\partial u}{\partial n} \,\mathrm{d}s = \int_D f \,\mathrm{d}x\,\mathrm{d}y,
$$

亦即

$$
-\left(\int_{\ell_1} + \int_{\ell_2} + \int_{\ell_3} + \int_{\ell_4}\right) \frac{\partial u}{\partial n} \,\mathrm{d}s = \int_D f \,\mathrm{d}x\,\mathrm{d}y.
$$

将上式左端边界积分在围绕内点的矩形控制体上逐边近似，于是积分分成四段进行，转化为定积分。

对于边 $\ell_1$，因为在矩形边 $\ell_1$ 上的外法向就是 $y$ 的负方向，弧长的微分 $\mathrm{d}s = \mathrm{d}x$，于是对此定积分用矩形公式近似计算，并且用差商代替微商，得到

$$
\begin{aligned}
\int_{\ell_1} \frac{\partial u}{\partial n} \,\mathrm{d}s &= \int_{\ell_1} \left(-\frac{\partial u}{\partial y}\right) \mathrm{d}x = -\frac{\partial u}{\partial y}\left(x_i, y_{j-\frac{1}{2}}\right) h + O(h^3) \\
&= -\frac{u_{ij} - u_{i,j-1}}{h}\,h + O(h^3).
\end{aligned}
$$

对于边 $\ell_2$，因为在矩形边 $\ell_2$ 上的外法向就是 $x$ 的正方向，弧长的微分 $\mathrm{d}s = \mathrm{d}y$，于是对此定积分用矩形公式近似计算，并且用差商代替微商，得到

$$
\begin{aligned}
\int_{\ell_2} \frac{\partial u}{\partial n} \,\mathrm{d}s &= \int_{\ell_2} \frac{\partial u}{\partial x} \,\mathrm{d}y = \frac{\partial u}{\partial x}\left(x_{i+\frac{1}{2}}, y_j\right) h + O(h^3) \\
&= \frac{u_{i+1,j} - u_{ij}}{h}\,h + O(h^3).
\end{aligned}
$$

类似地算出另外两条边上的积分如下：

$$
\begin{aligned}
\int_{\ell_3} \frac{\partial u}{\partial n} \,\mathrm{d}s &= \int_{\ell_3} \frac{\partial u}{\partial y} \,\mathrm{d}x = \frac{\partial u}{\partial y}\left(x_i, y_{j+\frac{1}{2}}\right) h + O(h^3) \\
&= \frac{u_{i,j+1} - u_{ij}}{h}\,h + O(h^3),
\end{aligned}
$$

$$
\begin{aligned}
\int_{\ell_4} \frac{\partial u}{\partial n} \,\mathrm{d}s &= \int_{\ell_4} \left(-\frac{\partial u}{\partial x}\right) \mathrm{d}y = -\frac{\partial u}{\partial x}\left(x_{i-\frac{1}{2}}, y_j\right) h + O(h^3) \\
&= -\frac{u_{ij} - u_{i-1,j}}{h}\,h + O(h^3).
\end{aligned}
$$

#### 积分法得到的五点差分格式

将四边近似代入 $\displaystyle -\int_{\partial D} \frac{\partial u}{\partial n} \,\mathrm{d}s \approx h^2 f_{ij}$（控制体上 $\displaystyle \int f\,\mathrm d x\,\mathrm d y$ 的近似），得

$$
4u_{ij} - (u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1}) = h^2 f_{ij},
$$

与 (1.1) 在正方形网格下的形式一致。

**注：** (1) 关于积分区域的选取：方法不同，得到的格式也有所不同；在此方法中积分区域是围绕节点的矩形。也可以直接用差分的网格作为积分区域，或者在非矩形网格上（例如三角形网格）分割积分区域。如下图：

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/c763622a-db61-41ee-8d7f-a09a4a0bdbf5/8321ebab10b339435ceb0117b8c9078d1f8381e7bc3074d266f03af4de6f3313.jpg)

(2) 采用不同的数值积分（例如梯形公式等）也可以得到不同的差分格式。

### 1.3 极坐标下的差分格式

如果求解域是圆环、环形域或扇形域，则采用极坐标是方便的，此时 Poisson 方程形如

$$
-\Delta_{r,\theta} u = -\left[ \frac{1}{r}\frac{\partial}{\partial r}\left(r\frac{\partial u}{\partial r}\right) + \frac{1}{r^2}\frac{\partial^2 u}{\partial \theta^2} \right] = f(r,\theta). \tag{1.2}
$$

方程 (1.2) 的系数于 $r = 0$ 奇异，因此只当 $r > 0$ 时有意义。为了定出我们感兴趣的解，需补充 $u$ 于 $r = 0$ 有界的条件，从而 $u$ 满足

$$
\lim_{r \to 0^+} r\frac{\partial u}{\partial r} = 0. \tag{1.3}
$$

关于变量 $r$、$\theta$ 分别取等步长 $h_r$ 和 $h_\theta$。令

$$
\begin{aligned}
r_i &= (i + 0.5) h_r, \quad i = 0, 1, 2, \ldots, \\
\theta_j &= (j + 1) h_\theta, \quad j = 0, 1, \ldots, J-1, \quad h_\theta = 2\pi/J.
\end{aligned}
$$

对任一点 $(r_i, \theta_j)$（$i > 0$），用中心差商公式

$$
\begin{aligned}
\left[ \frac{1}{r}\frac{\partial}{\partial r}\left(r\frac{\partial u}{\partial r}\right) \right]_{(r_i, \theta_j)} &\approx \frac{1}{r_i}\frac{r_{i+\frac{1}{2}} u_{i+1,j} - (r_{i+\frac{1}{2}} + r_{i-\frac{1}{2}}) u_{ij} + r_{i-\frac{1}{2}} u_{i-1,j}}{h_r^2}, \\
\left[ \frac{1}{r^2}\frac{\partial^2 u}{\partial \theta^2} \right]_{(r_i, \theta_j)} &\approx \frac{1}{r_i^2}\frac{u_{i,j+1} - 2u_{ij} + u_{i,j-1}}{h_\theta^2}, \\
(u_{i0} = u_{iJ}, &\quad u_{i,-1} = u_{i,J-1}).
\end{aligned}
$$

代到 (1.2)，则就 $i \ge 1$ 得到逼近它的差分方程

$$
\begin{aligned}
&-\left[ \frac{1}{r_i}\frac{r_{i+\frac{1}{2}} u_{i+1,j} - (r_{i+\frac{1}{2}} + r_{i-\frac{1}{2}}) u_{ij} + r_{i-\frac{1}{2}} u_{i-1,j}}{h_r^2} \right. \\
&\quad \left. + \frac{1}{r_i^2}\frac{u_{i,j+1} - 2u_{ij} + u_{i,j-1}}{h_\theta^2} \right] = f(r_i, \theta_j).
\end{aligned}
$$

今用积分插值法导出点 $(r_0, \theta_j)$ 的差分方程：以 $r$ 乘方程 (1.2)，并对 $r$ 由 $\varepsilon$ 到 $h_r$ 积分，对 $\theta$ 由 $\theta_{j-\frac{1}{2}}$ 到 $\theta_{j+\frac{1}{2}}$ 积分，然后令 $\varepsilon \to 0$，且注意条件 (1.3)，则得

$$
\begin{aligned}
&-\left[ h_r \int_{\theta_{j-\frac{1}{2}}}^{\theta_{j+\frac{1}{2}}} \frac{\partial u}{\partial r}(h_r, \theta) \,\mathrm{d}\theta + \int_0^{h_r} \frac{1}{r}\left( \frac{\partial u}{\partial \theta}(r, \theta_{j+\frac{1}{2}}) - \frac{\partial u}{\partial \theta}(r, \theta_{j-\frac{1}{2}}) \right) \mathrm{d}r \right] \\
&= \int_0^{h_r} r \,\mathrm{d}r \int_{\theta_{j-\frac{1}{2}}}^{\theta_{j+\frac{1}{2}}} f(r,\theta) \,\mathrm{d}\theta .
\end{aligned}
$$

用中矩形公式代替上述积分，则

$$
\begin{aligned}
&-\left[ h_r h_\theta \frac{\partial u}{\partial r}(h_r, \theta_j) + 2\left(\frac{\partial u}{\partial \theta}\left(\frac{h_r}{2}, \theta_{j+\frac{1}{2}}\right) - \frac{\partial u}{\partial \theta}\left(\frac{h_r}{2}, \theta_{j-\frac{1}{2}}\right) \right) \right] \\
&\approx \frac{1}{2} h_r^2 h_\theta f\left(\frac{h_r}{2}, \theta_j\right).
\end{aligned}
$$

再用中心差商代替微商，就得出点 $(r_0, \theta_j)$ 的差分方程：

$$
-\left[ h_r h_\theta \frac{u_{1j} - u_{0j}}{h_r} + 2\frac{u_{0,j+1} - 2u_{0j} + u_{0,j-1}}{h_\theta} \right] = \frac{1}{2} h_r^2 h_\theta f_{0j}.
$$

或（两端除以 $\frac{1}{2} h_r^2 h_\theta$）

$$
-\left[ \frac{2}{h_r}\frac{u_{1j} - u_{0j}}{h_r} + \frac{4}{h_r^2}\frac{u_{0,j+1} - 2u_{0j} + u_{0,j-1}}{h_\theta^2} \right] = f_{0j}.
$$
