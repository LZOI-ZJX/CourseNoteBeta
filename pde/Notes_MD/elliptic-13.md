# 第五章 椭圆型方程的差分方法（续）

## 第二节 差分格式的性质

### 2.1 矩形区域上的 Poisson 方程的五点差分逼近的敛速分析

考虑泊松方程的第一边值问题

$$
\left\{ \begin{array}{l l} \frac {\partial^ {2} u}{\partial x ^ {2}} + \frac {\partial^ {2} u}{\partial y ^ {2}} = - f (x, y) & (x, y) \in \Omega \\ u (x, y) = g (x, y) & (x, y) \in \partial \Omega \end{array} \right. \tag {2.1}
$$

其中 $\Omega = \{(x, y) \mid 0 < x < 1, 0 < y < 1\}, \Omega$ 中的网格点集合为 $\Omega_{h}$ , 而

$$
\Omega_ {h} = \left\{\left(x _ {l}, y _ {m}\right) \mid x _ {l} = l h, y _ {m} = m h; l, m = 1, 2, \dots , M - 1 \right\}
$$

边界上网格点集合为 $\partial \Omega_{h}$ , 而

$$
\overline {{\Omega}} = \Omega \cup \partial \Omega , \quad \overline {{\Omega}} _ {h} = \Omega_ {h} \cup \partial \Omega_ {h}
$$

网格点 $\left(x_{l}, y_{m}\right) \quad (l, m) \in \Omega_{h}$ 上的五点差分逼近是

$$
\frac {1}{h ^ {2}} \left(U _ {l + 1, m} + U _ {l - 1, m} + U _ {l, m + 1} + U _ {l, m - 1} - 4 U _ {l, m}\right) = - f _ {l, m}
$$

写成 $L U_{l,m} = -f_{l,m}\quad f_{l,m} = f\big(x_l,y_m\big);\quad (l,m)\in \Omega_h$ (2.2)

在边界网格点 $(l, m)$ 上满足边值条件

$$
U _ {l, m} = g _ {l, m} \quad g _ {l, m} = g (x _ {l}, y _ {m}); (l, m) \in \partial \Omega_ {h} (2. 3)
$$

为了讨论差分方程解 $U_{l, m}$ 与微分方程解 $u\left(x_{l}, y_{m}\right)$ 之间的逼近程度，令

$$
\boldsymbol {e} _ {l, m} = \boldsymbol {u} _ {l, m} - \boldsymbol {U} _ {l, m}
$$

于是

$$
L e _ {l, m} = L u _ {l, m} - L U _ {l, m} = L u _ {l, m} + f _ {l, m} = T _ {l, m}
$$

由Taylor展开，设 $u \in C^{4,4}(\overline{\Omega})$ ，于是

$$
T _ {l, m} = \frac {1}{1 2} h ^ {2} \left\{\left(\frac {\partial^ {4} u}{\partial x ^ {4}}\right) _ {\xi , y _ {m}} + \left(\frac {\partial^ {4} u}{\partial y ^ {4}}\right) _ {x _ {l}, \eta} \right\} \tag {2.4}
$$

其中 $x_{l} - h < \xi < x_{l} + h,\quad y_{m} - h < \eta < y_{m} + h$ ，并假定 $|\partial_x^4 u|$ 与 $|\partial_y^4 u|$ 在 $\overline\Omega$ 上有共同上界 $M_{4}$ ，则由 (2.4) 与 $M_4$ 的定义得

$$
\left| L e _ {l, m} \right| = \left| T _ {l, m} \right| \leq \frac {1}{12} h^ {2} \left(M _ {4} + M _ {4}\right) = \frac {1}{6} h^ {2} M _ {4}
$$

于是 $\max\limits_{\Omega_h}\left|Le_{l,m}\right| = \max\limits_{\Omega_h}\left|T_{l,m}\right| \leq \dfrac{1}{6} h^2 M_4$

为了估计 $\max_{\Omega_h} \left| e_{l,m} \right|$ , 建立它和 $\max_{\Omega_h} \left| L e_{l,m} \right|$ 之间的关系。

#### 定理 2.1（极值原理）

假设

(1) $V_{l, m}$ 为定义在网格点区域 $\overline{\Omega}_{h} = \Omega_{h} \cup \partial \Omega_{h}$ 上的一组值，

(2) $V_{l, m}$ 不恒为常数, $(l, m) \in \overline{\Omega}_{h}$ ,

(3) $L V_{l, m} \geq 0, (l, m) \in \Omega_{h}$ . 

则 $V_{l, m}$ 不可能在内部结点即 $\Omega_ {h}$ 上达到最大值。同样，如果 $L V_{l, m} \leq 0$ ，则 $V_{l, m}$ 在 $\partial \Omega_ {h}$ 上达到最小值。

#### 估计式：$\max_{\Omega_h}|V|$ 与 $|LV|$ 的关系

若 $V$ 是定义在矩形区域 $0 \leq x \leq 1, 0 \leq y \leq 1$ 上的网格点集 $\overline{\Omega}_{h} = \Omega_{h} \cup \partial \Omega_{h}$ 上的函数，则有

$$
\max _ {\Omega_ {h}} | V | \leq \max _ {\partial \Omega_ {h}} | V | + \frac {1}{2} \max _ {\Omega_ {h}} | L V |
$$

其中 $L V_{l,m} = \frac{1}{h^2}\left(V_{l+1,m} + V_{l-1,m} + V_{l,m+1} + V_{l,m-1} - 4V_{l,m}\right)$ 

应用定理 2.1 的结论，则差分方程解与微分方程解之差 $e_{l, m}$ 满足

$$
\max _ {\Omega_ {h}} \left| e _ {l, m} \right| \leq \max _ {\partial \Omega_ {h}} \left| e _ {l, m} \right| + \frac {1}{2} \max _ {\Omega_ {h}} \left| L e _ {l, m} \right|
$$

在边界网格点 $\partial \Omega_{h}$ 上 $U_{l, m} = u_{l, m} = g_{l, m}$ ，于是由 $e_{l,m}|_{\partial\Omega_h}=0$ 与上式得到

$$
\max  _ {\Omega_ {h}} \left| e _ {l, m} \right| \leq \frac {1}{1 2} h ^ {2} M _ {4} \tag {2.5}
$$

它证明了当网距 $h \rightarrow 0$ 差分方程解收敛到微分方程解, 式(2.5)是敛速估计, 我们有定理：

**定理 2.3** 若 Poisson 方程的第一边值问题的解 $u(x, y) \in C^{4,4}(\overline{\Omega})$ ，则五点差分公式的解一致收敛到 $u(x, y)$ ，且有敛速估计

$$
\max _ {\Omega_ {h}} \left| u \left(x _ {l}, y _ {m}\right) - U _ {l, m} \right| \leq \frac {1}{1 2} M _ {4} h ^ {2}
$$

其中 $M_{4}$ 为与 $h$ 无关的常数。

## 第三节 边界条件的处理

### 3.1 正方形区域中的 Laplace 方程 Dirichlet 边值问题的差分模拟

考虑Laplace方程

$$
\frac {\partial^ {2} u}{\partial x ^ {2}} + \frac {\partial^ {2} u}{\partial y ^ {2}} = 0 (x, y) \in \Omega \tag {3.1}
$$

设 $\Omega$ 为正方形区域, $0 < x < 1, 0 < y < 1$ , 求方程 (3.1)满足边值条件

$$
u (x, y) = f (x, y) \quad (x, y) \in \partial \Omega \tag {3.2}
$$

的解。

#### Laplace 方程的五点差分格式

$$
\frac {1}{h ^ {2}} \left(U _ {l + 1, m} + U _ {l - 1, m} + U _ {l, m + 1} + U _ {l, m - 1} - 4 U _ {l, m}\right) = 0 \tag {3.3}
$$

它具有截断误差: $\frac{1}{12} h^{2}\left(\frac{\partial^{4} u}{\partial x^{4}} + \frac{\partial^{4} u}{\partial y^{4}}\right)_{l, m} + \cdots$ 

我们引进记号◇，有

$$
\diamond U _ {l, m} = \frac {1}{h ^ {2}} \left(U _ {l + 1, m} + U _ {l - 1, m} + U _ {l, m + 1} + U _ {l, m - 1} - 4 U _ {l, m}\right) \tag {3.4}
$$

因此差分方程 (3.3) 即 $\diamond U_{l, m} = 0$ 。

#### 内部结点上的差分方程（图 3.1）

$$
\left(l = 1, \dots , M - 1; m = 1, \dots , M - 1\right)
$$

由此在区域 $\Omega$ 内部 $(M - 1)^{2}$ 个点上建立 $(M - 1)^{2}$ 个方程。

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/5deb179f-0b3a-4964-99cc-26921ab26250/ef9f7e293b153970f7e2be1bf90471a604c5c75a71e137b5621a1d0ef251724a.jpg)



图3.1


$$
\diamond U _ {l, m} = 0 \quad l, m = 1, \dots , M - 1
$$

(3.8) 

#### 向量记号与矩阵方程 $AU = K$

$$
U = \left[ U _ {1, 1}, U _ {2, 1}, \dots , U _ {M - 1, 1}; U _ {1, 2}, U _ {2, 2}, \dots , U _ {M - 1, 2}; U _ {1, M - 1}, U _ {2, M - 1}, \dots , U _ {M - 1, M - 1} \right] ^ {T}
$$

单位正方形中的内部结点上的 $(M - 1)^{2}$ 个线性方程 (3.8)写成矩阵形式为

$$
A U = K \tag {3.9}
$$

其中， $A$ 是 $(M - 1)^2$ 阶方阵

$$
A = \left[ \begin{array}{c c c c c} B & - I & & & \\ - I & B & - I & & \\ & \ddots & \ddots & \ddots & \\ & & - I & B & - I \\ & & & - I & B \end{array} \right]
$$

$I$ 是 $(M - 1)$ 阶单位方阵; $B$ 是 $(M - 1)$ 阶方阵。

$$
B = \left[ \begin{array}{c c c c c} 4 & - 1 & & & \\ - 1 & 4 & - 1 & & \\ & \ddots & \ddots & \ddots & \\ & & - 1 & 4 & - 1 \\ & & & - 1 & 4 \end{array} \right]
$$

其中 $K$ 是由边界条件 $f(x,y)$ 决定的 $(M-1)^2$ 维已知列向量。它的非零元素主要分布在靠近边界 $\partial \Omega$ 的节点所对应的位置上。

### 3.2 Neumann 边值问题的差分模拟

现在我们考虑 Laplace 方程 Neumann 边值问题，即

$$
\left\{ \begin{array}{l} \frac {\partial^ {2} u}{\partial x ^ {2}} + \frac {\partial^ {2} u}{\partial y ^ {2}} = 0 \quad (x, y) \in \Omega ; \Omega = \{(x, y) | 0 <   x <   1, 0 <   y <   1 \} \\ \frac {\partial u}{\partial n} | _ {\partial \Omega} = g (x, y) \end{array} \right. \tag {3.10}
$$

$\frac{\partial u}{\partial n}\bigg|_{\partial \Omega}$ 表示函数 $u$ 沿着边界的外法线方向导数。

#### Neumann 问题 (3.10) 的可解性条件

$$
\int_ {\partial \Omega} g (x, y) d l = 0
$$

且除了一个任意常数外，解唯一。因为容易看到，如果 $u(x, y)$ 是式(3.10)的解，于是， $u(x, y)+C(C$ 是一个任意常数)也是其解。为了唯一性，需要规定 $u(x, y)$ 在区域中某一点上的值。

#### Neumann 边值问题的差分模拟

先在区域 $\Omega$ 中给定一个正方形网格区域，步长为 $h, M h = 1$ , 于是必须确定解的结点为 $(M + 1)^{2}$ 个, 网格点上的差分方程的解为

$$
U _ {l, m} \left(0 \leq l, m \leq M\right) \quad (3. 1 2)
$$

在内点上，Laplace 方程由差分方程 (3.3) 代替：

$$
\frac {1}{h ^ {2}} \left(\delta_ {x} ^ {2} + \delta_ {y} ^ {2}\right) U _ {l, m} = 0
$$

$$
l, m = 1, 2, \dots , M - 1
$$

#### $x = 0$ 上的导数边值条件的差分模拟

$$
\frac {1}{2 h} \left(U _ {- 1, m} - U _ {1, m}\right) = g _ {0, m} \quad m = 1, 2, \dots M - 1 \tag {3.13}
$$

这里 $g_{0,m} = g(0,mh)$。 

#### 在五点差分格式 (3.12) 中令 $l = 0$

$$
\left(\delta_ {x} ^ {2} + \delta_ {y} ^ {2}\right) U _ {0, m} = 0
$$

即 $U_{-1,m} = 2U_{0,m} - U_{1,m} - U_{0,m-1} + 2U_{0,m} - U_{0,m+1}$ 

代入式 (3.13)，则

$$
4 U _ {0, m} - 2 U _ {1, m} - U _ {0, m - 1} - U _ {0, m + 1} = 2 h g _ {0, m} \quad m = 1, \dots , M - 1
$$

( 3.14 ) 

#### $x = 1,\ y = 0,\ y = 1$ 上的边界方程

$$
4 U _ {M, m} - 2 U _ {M - 1, m} - U _ {M, m - 1} - U _ {M, m + 1} = 2 h g _ {M, m} \quad m = 1, \dots , M - 1 \tag {3.15}
$$

$$
4 U _ {l, 0} - 2 U _ {l, 1} - U _ {l - 1, 0} - U _ {l + 1, 0} = 2 h g _ {l, 0} \quad l = 1, \dots , M - 1 \tag {3.16}
$$

$$
4 U _ {l, M} - 2 U _ {l, M - 1} - U _ {l - 1, M} - U _ {l + 1, M} = 2 h g _ {l, M} \quad l = 1, \dots , M - 1 \tag {3.17}
$$

#### 四个角点上的方程

$$
4 U _ {0, 0} - 2 U _ {1, 0} - 2 U _ {0, 1} = 4 h g _ {0, 0}
$$

$$
4 U _ {0, M} - 2 U _ {1, M} - 2 U _ {0, M - 1} = 4 h g _ {0, M}
$$

$$
4 U _ {M, 0} - 2 U _ {M, 1} - 2 U _ {M - 1, 0} = 4 h g _ {M, 0}
$$

$$
4 U _ {M, M} - 2 U _ {M - 1, M} - 2 U _ {M, M - 1} = 4 h g _ {M, M}
$$

由此, 正方形 $0 \leq x, y \leq 1$ 区域的 $(M + 1)^{2}$ 个结点上差分方程解 $U_{lm}$ 满足线性方程组

$$
A U = 2 h g \tag {3.18}
$$

这里 $A$ 是 $(M + 1)^{2}$ 阶方阵

$$
A = \left[ \begin{array}{c c c c c} B & - 2 I & & & \\ - I & B & - I & & \\ & \ddots & \ddots & \ddots & \\ & & - I & B & - I \\ & & & - 2 I & B \end{array} \right]
$$

$I$ 是 $(M + 1)$ 阶单位方阵； $B$ 是如下 $(M + 1)$ 阶方阵：

$$
B = \left[ \begin{array}{c c c c c} 4 & - 2 & & & \\ - 1 & 4 & - 1 & & \\ & \ddots & \ddots & \ddots & \\ & & - 1 & 4 & - 1 \\ & & & - 2 & 4 \end{array} \right]
$$

#### 方程组 (3.18) 中的向量 $U$ 和 $g$

$$
U = \left[ U _ {0, 0}, U _ {1, 0}, \dots , U _ {M, 0}; U _ {0, 1}, \dots , U _ {M, 1}; \dots ; U _ {0, M}, \dots , U _ {M, M} \right] ^ {T}
$$

$$
\begin{array}{l} g = \left[ 2 g _ {0, 0}, g _ {1, 0}, \dots , g _ {M - 1, 0}, 2 g _ {M, 0}; g _ {0, 1}, 0, \dots , 0, g _ {M, 1}; \right. \\ \dots ; g _ {0, M - 1}, 0, \dots , 0, g _ {M, M - 1}; 2 g _ {0, M}, g _ {1, M}, \dots , g _ {M - 1, M}, 2 g _ {M, M} \bigg ] ^ {T} \\ \end{array}
$$

$$
g _ {l, m} = g \big (l h, m h \big)
$$

例 3.1 在单位正方形区域 $\Omega$ 上解Laplace方程Neumann问题

$$
\left\{ \begin{array}{l} \frac {\partial^ {2} u}{\partial x ^ {2}} + \frac {\partial^ {2} u}{\partial y ^ {2}} = 0 \quad 0 \leq x, y \leq 1 \\ \left. \frac {\partial u}{\partial n} \right| _ {\partial \Omega} = g (x, y) \end{array} \right.
$$

解 令 $h = 1 / 2$ , 应用图3.2中结点次序, 则方程(3.18)为

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/5deb179f-0b3a-4964-99cc-26921ab26250/494f846fdfbf1955e694a2ad8677d241fcd4e551cf78aba4860fc9c292f15791.jpg)



图3.2


$$
\begin{bmatrix}
4 & -2 & 0 & -2 & 0 & 0 & 0 & 0 & 0 \\
-1 & 4 & -1 & 0 & -2 & 0 & 0 & 0 & 0 \\
0 & -2 & 4 & 0 & 0 & -2 & 0 & 0 & 0 \\
-1 & 0 & 0 & 4 & -2 & 0 & -1 & 0 & 0 \\
0 & -1 & 0 & -1 & 4 & -1 & 0 & -1 & 0 \\
0 & 0 & -1 & 0 & -2 & 4 & 0 & 0 & -1 \\
0 & 0 & 0 & -2 & 0 & 0 & 4 & -2 & 0 \\
0 & 0 & 0 & 0 & -2 & 0 & -1 & 4 & -1 \\
0 & 0 & 0 & 0 & 0 & -2 & 0 & -2 & 4
\end{bmatrix}
\left[ \begin{array}{l} U _ {1} \\ U _ {2} \\ U _ {3} \\ U _ {4} \\ U _ {5} \\ U _ {6} \\ U _ {7} \\ U _ {8} \\ U _ {9} \end{array} \right] = 2 h \left[ \begin{array}{l} 2 g _ {1} \\ g _ {2} \\ 2 g _ {3} \\ g _ {4} \\ 0 \\ g _ {6} \\ 2 g _ {7} \\ g _ {8} \\ 2 g _ {9} \end{array} \right]
$$

(3.19) 

或简写成 $A U = 2 h g$ 。显然 $A$ 是一奇异矩阵。

### 3.3 混合边值问题

在 $xy$ 平面的某区域 $\Omega$ 中, 未知函数 $u$ 满足Laplace方程, 将边界 $\partial \Omega$ 分成若干弧段, 要求 $u$ 在每一弧段上满足不同类型的边界条件。讨论此类定解问题的差分模拟。

例如，求解如下定解问题：

$$
\left\{ \begin{array}{l} \frac {\partial^ {2} u}{\partial x ^ {2}} + \frac {\partial^ {2} u}{\partial y ^ {2}} = 0 \quad (x, y) \in \Omega ; \Omega = \{(x, y) | 0 <   x, y <   1 \} \\ \frac {\partial u}{\partial x} - p (y) u = f _ {0} (y) \quad x = 0; 0 \leq y \leq 1 \\ \frac {\partial u}{\partial y} - q (x) u = f _ {1} (x) \quad y = 0; 0 \leq x \leq 1 \\ u = g (x, y) \quad \left\{ \begin{array}{l} x = 1; 0 \leq y \leq 1 \\ y = 1; 0 \leq x \leq 1 \end{array} \right. \end{array} \right.
$$

(3.21) 

$$
p, q, f _ {0}, f _ {1}, g \text{ 是给定的函数。}
$$

#### 内部结点：五点差分公式

在求解区域 $\Omega$ 内由逼近 Laplace 方程的五点差分公式

$$
\frac {1}{h ^ {2}} \left(\delta_ {x} ^ {2} + \delta_ {y} ^ {2}\right) U _ {l, m} = 0, l, m = 1, 2, \dots , M - 1
$$

给出函数 $u$ 在结点 $(lh, mh)$ 的近似值 $U_{lm}$ 所满足的差分方程。

#### $x = 0$ 边界：边值条件的离散

对于在 $x = 0$ 上的结点 $(0, m h)$ ，应用边值条件

$$
\frac {\partial u}{\partial x} - p (y) u = f _ {0} (y)
$$

联立差分模拟与五点差分公式，即

$$
U _ {- 1, m} + U _ {1, m} + U _ {0, m + 1} + U _ {0, m - 1} - 4 U _ {0, m} = 0
$$

$$
\left(U _ {1, m} - U _ {- 1, m}\right) / (2 h) - p _ {m} U _ {0, m} = f _ {0, m} \quad p _ {m} = p (m h); f _ {0, m} = f _ {0} (m h) \tag {3.22}
$$

消去 $U_{-1, m}$ , 得

$$
2 U _ {1, m} + U _ {0, m + 1} + U _ {0, m - 1} - \left(4 + 2 h p _ {m}\right) U _ {0, m} = 2 h f _ {0, m} \quad m = 1, 2, \dots , M - 1
$$

相似地对于 $y = 0$ 上的结点 $(l h, 0)$ , 我们有

$$
2 U _ {l, 1} + U _ {l + 1, 0} + U _ {l - 1, 0} - \left(4 + 2 h q _ {l}\right) U _ {l, 0} = 2 h f _ {1, l} \quad l = 1, 2, \dots , M - 1 \tag {3.23}
$$

其中 $q_{l} = q(lh), f_{1,l} = f_{1}(lh)$ 。

#### 原点 $(0,0)$ 处两边界相遇

$$
\begin{array}{l} U _ {- 1, 0} + U _ {1, 0} + U _ {0, - 1} + U _ {0, 1} - 4 U _ {0, 0} = 0 \\ U _ {1, 0} - U _ {- 1, 0} - 2 h p _ {0} U _ {0, 0} = 2 h f _ {0, 0} \\ U _ {0, 1} - U _ {0, - 1} - 2 h q _ {0} U _ {0, 0} = 2 h f _ {1, 0} \\ \end{array}
$$

消去 $U_{-1,0}$ 和 $U_{0,-1}$ , 则

$$
2 U _ {1, 0} + 2 U _ {0, 1} - \left(4 + 2 h p _ {0} + 2 h q _ {0}\right) U _ {0, 0} = 2 h \left(f _ {0, 0} + f _ {1, 0}\right) \tag {3.24}
$$

且对 $l = 0$ 和 $m = 0$ 上成立的方程(3.22)，(3.23)用 $1 / 2$ 乘之，对 $l = m = 0$ 上的方程(3.24)用 $1 / 4$ 乘之。

#### 整体差分方程组

$$
\begin{array}{l} U _ {l - 1, m} + U _ {l + 1, m} + U _ {l, m + 1} + U _ {l, m - 1} - 4 U _ {l, m} = 0 \quad l, m = 1, 2, \dots , M - 1 \\ U _ {l, 1} + \left. U _ {l + 1, 0} \right/ _ {2} + \left. U _ {l - 1, 0} \right/ _ {2} - (2 + h q _ {l}) U _ {l, 0} = h f _ {1, l} \quad l = 1, 2, \dots , M - 1 \\ U _ {1, m} + \left. U _ {0, m + 1} \right/ _ {2} + \left. U _ {0, m - 1} \right/ _ {2} - (2 + h p _ {m}) U _ {0, m} = h f _ {0, m} \quad m = 1, 2, \dots , M - 1 \\ U _ {1, 0} / 2 + U _ {0, 1} / 2 - \left(2 + h p _ {0} + h q _ {0}\right) / 2 U _ {0, 0} = h \left(f _ {0, 0} + f _ {1, 0}\right) / 2 \\ \end{array}
$$

令

$$
U = \left[ U _ {0, 0}, U _ {1, 0}, \dots , U _ {M - 1, 0}; U _ {0, 1}, \dots , U _ {M - 1, 1}; \dots , U _ {0, M - 1}, \dots , U _ {M - 1, M - 1} \right] ^ {T}
$$

下面把所有差分方程写成矩阵形式，于是 $U$ 满足方程

$$
A U = h g \tag {3.25}
$$

#### 矩阵 $A$（$M^{2}$ 阶对称方阵）

$$
A = \left[ \begin{array}{l l l l l l} E _ {0} & K & & & & \\ K & E _ {1} & K & & & \\ & K & E _ {2} & K & & \\ & & \ddots & \ddots & \ddots & \\ & & & K & E _ {M - 2} & K \\ & & & & K & E _ {M - 1} \end{array} \right]
$$

$$
E _ {m} = \left[ \begin{array}{c c c c c} - (2 + h p _ {m}) & 1 & & & \\ 1 & - 4 & 1 & & \\ & \ddots & \ddots & \ddots & \\ & & 1 & - 4 & 1 \\ & & & 1 & - 4 \end{array} \right] m = 1, 2, \dots M - 1
$$

$E_{0}, E_{1}, \dots, E_{M-1}, K$ 是 $M$ 阶方阵。其中

$$
E _ {0} = \left[ \begin{array}{c c c c c} - \left[ 1 + \frac {1}{2} h \left(p _ {0} + q _ {0}\right) \right] & \frac {1}{2} & & & \\ \frac {1}{2} & - (2 + h q _ {1}) & \frac {1}{2} & & \\ & \ddots & \ddots & \ddots & \\ & & \frac {1}{2} & - (2 + h q _ {M - 2}) & \frac {1}{2} \\ & & & \frac {1}{2} & - (2 + h q _ {M - 1}) \end{array} \right]
$$

$$
K = \left[ \begin{array}{l l l l} \frac {1}{2} & & & \\ & 1 & & \\ & & \ddots & \\ & & & 1 \end{array} \right]
$$

而 $g$ 依赖于函数 $f_{0}, f_{1}$ 和 $g$ 在边界结点上的值。

### 3.4 一般区域边界条件的处理

#### 第一类边界条件

方法一：直接转移，找最近的节点代替边界节点。

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/5deb179f-0b3a-4964-99cc-26921ab26250/25a1ec81afa073ed84b6f5f35ecd67082b5a704c575f4a472c0fdd1fdf2d0280.jpg)


$$
u (S) = f (S)
$$

$$
u (P) \approx f (T)
$$

方法二：线性插值法。

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/5deb179f-0b3a-4964-99cc-26921ab26250/2fe50439c523ce8b3d4d6276f80401a793c115577ba2de660cc84cb976050c05.jpg)


在 $P$ 点可以用 $T, Q$ 做两点线性插值, 设 $T P = \delta$ 

则 $u_{P} = \frac{h}{h + \delta} u_{T} + \frac{\delta}{h + \delta} u_{Q}$ 

#### 第三（Robin）边界条件

**情形一：** 边界节点在求解区域的边界上。

见图 5.7、5.8。

如果外法线与坐标轴平行, 即

$$
\frac {\partial u}{\partial n} = \pm \frac {\partial u}{\partial x}, \quad \text{或} \quad \frac {\partial u}{\partial n} = \pm \frac {\partial u}{\partial y}
$$

取 $\frac{u(P) - u(Q)}{h} + \gamma u(P) = \beta(P)$ 

**外法线与坐标轴不平行（见图 5.9）。**

$$
\begin{array}{l} \frac {\partial u (P)}{\partial n} = \frac {\partial u}{\partial x} \cos (n, x) + \frac {\partial u}{\partial y} \cos (n, y) \\ = - \frac {\partial u}{\partial x} \cos \theta_ {1} - \frac {\partial u}{\partial y} \cos \theta_ {2} \\ \end{array}
$$

$$
\begin{array}{l} \frac {u (P) - u (Q)}{h} \cos (\theta_ {1}) + \frac {u (P) - u (R)}{h} \cos (\theta_ {2}) \\ + \gamma (P) u (P) = \beta (P) \\ \end{array}
$$

**情形二：** 边界节点不在求解区域的边界上（见图 5.10）。

$$
S T = \delta_ {3} h, T R = \delta_ {2} h, P T = \delta_ {1} h
$$

用边界结点处的方向导数近似外法向导数：

$$
\frac {u (P) - u _ {T}}{\delta_ {1} h} + \gamma (P) u (P) = \beta (P)
$$

其中

$$
u _ {T} = \delta_ {3} u _ {R} + \delta_ {2} u _ {S}
$$