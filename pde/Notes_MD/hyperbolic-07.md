## 第五节 双曲型方程及方程组的初边值问题

### 5.1 二阶双曲型方程的边界处理

#### 定解问题与网格剖分

波动方程：

$$
\frac {\partial^ {2} \boldsymbol {u}}{\partial t ^ {2}} = \boldsymbol {a} ^ {2} \frac {\partial^ {2} \boldsymbol {u}}{\partial x ^ {2}} \qquad \quad \boldsymbol {0} <   \boldsymbol {x} <   l, t > \boldsymbol {0}
$$

边界条件:

$$
\begin{array}{l} \boldsymbol {\mathrm{I}}: \boldsymbol {u} (\mathbf {0}, t) = \boldsymbol {\mu} _ {0} (t), \quad \boldsymbol {u} (l, t) = \boldsymbol {\mu} _ {1} (t) \\ \boldsymbol {\mathrm{II}}: \left[ \frac {\partial u}{\partial x} + \alpha_ {0} (t) u \right] | _ {x = 0} = \eta_ {0} (t), \left[ \frac {\partial u}{\partial x} + \alpha_ {1} (t) u \right] | _ {x = l} = \eta_ {1} (t) \\ \end{array}
$$

若取网格:

$$
x _ {j} = j h, j = 0, 1, \dots , J, t _ {n} = n \tau , h = l / J
$$

#### 第一类边界条件与零阶离散

$$
I: \boldsymbol {u} _ {0} ^ {n} = \mu_ {0} (\boldsymbol {n} \tau) = \mu_ {0} ^ {n}, \quad \boldsymbol {u} _ {J} ^ {n} = \mu_ {1} (\boldsymbol {n} \tau) = \mu_ {1} ^ {n}
$$

#### 第二类边界条件（一阶差商）

若用向前和向后差商近似 $\left[\frac{\partial u}{\partial x}\right]_{x=0}$ 和 $\left[\frac{\partial u}{\partial x}\right]_{x=l}$ 

$$
I I: \left\{ \begin{array}{l} \frac {u _ {1} ^ {n} - u _ {0} ^ {n}}{h} + \alpha_ {0} ^ {n} u _ {0} ^ {n} = \eta_ {0} ^ {n} \\ \frac {u _ {J} ^ {n} - u _ {J - 1} ^ {n}}{h} + \alpha_ {1} ^ {n} u _ {J} ^ {n} = \eta_ {1} ^ {n} \end{array} \right.
$$

#### 第二类边界条件（二阶精度）：左端点 $x=0$

在边界用三点组合逼近 $u_x(0)$：

$$
\frac {\partial u}{\partial x} \bigg | _ {x = 0} = \frac {a u _ {0} + b u _ {1} + c u _ {2}}{h} + o (h ^ {2})
$$

确定 $a, b, c$

$$
u (h) = u (0) + h \frac {\partial u}{\partial x} \big | _ {x = 0} + \frac {h ^ {2}}{2} \frac {\partial^ {2} u}{\partial x ^ {2}} \big | _ {x = 0} + O \left(h ^ {3}\right)
$$

$$
u (2 h) = u (0) + 2 h \frac {\partial u}{\partial x} \big | _ {x = 0} + \frac {(2 h) ^ {2}}{2} \frac {\partial^ {2} u}{\partial x ^ {2}} \big | _ {x = 0} + O \left(h ^ {3}\right)
$$

得方程组

$$
\left\{ \begin{array}{l} a + b + c = 0 \\ b + 2 c = 1; \\ b / 2 + 2 c = 0 \end{array} \right.
$$

$$
a = - \frac {3}{2}, b = 2, c = - \frac {1}{2}
$$

#### 第二类边界条件（二阶精度）：右端点 $x=l$

在边界用向内三点组合逼近 $u_x(l)$（$u_J=u(l)$，$u_{J-1}=u(l-h)$，$u_{J-2}=u(l-2h)$）：

$$
\frac {\partial u}{\partial x} \big | _ {x = l} = \frac {a u _ {J} + b u _ {J - 1} + c u _ {J - 2}}{h} + o \left(h ^ {2}\right)
$$

将 $u(l-h)$、$u(l-2h)$ 在 $x=l$ 处 Taylor 展开并代入 $\frac{1}{h}\bigl(a u_J + b u_{J-1} + c u_{J-2}\bigr)$，要求 $u(l)$、$u_{xx}(l)$ 的系数为 $0$、$u_x(l)$ 的系数为 $1$，即得下列方程组。

$$
u (l - h) = u (l) - h \frac {\partial u}{\partial x} \big | _ {x = l} + \frac {h ^ {2}}{2} \frac {\partial^ {2} u}{\partial x ^ {2}} \big | _ {x = l} + O (h ^ {3})
$$

$$
u (l - 2 h) = u (l) - 2 h \frac {\partial u}{\partial x} \big | _ {x = l} + \frac {(2 h) ^ {2}}{2} \frac {\partial^ {2} u}{\partial x ^ {2}} \big | _ {x = l} + O (h ^ {3})
$$

得方程组

$$
\left\{ \begin{array}{l} a + b + c = 0 \\ - b - 2 c = 1; \\ b / 2 + 2 c = 0 \end{array} \right.
$$

$$
a = \frac {3}{2}, b = - 2, c = \frac {1}{2}
$$

#### 第二类边界条件（二阶精度）：离散式汇总

右边界：$u_x|_{x=l} = \dfrac{3 u_J - 4 u_{J-1} + u_{J-2}}{2h} + O(h^2)$，与 $\dfrac{1}{h}\bigl(\frac{3}{2} u_J - 2 u_{J-1} + \frac{1}{2} u_{J-2}\bigr)$ 等价。分子须与左端对称：边界上 $u_J$（或 $u_0$）系数为 $+3$，向内相邻点为 $-4$，再向内为 $+1$。**不可**写成 $3 u_{J-2} - 4 u_{J-1} + u_J$（否则把 $u_J$ 与 $u_{J-2}$ 的系数弄反）。

$$
I I: \left\{ \begin{array}{l} \frac {- 3 u _ {0} ^ {n} + 4 u _ {1} ^ {n} - u _ {2} ^ {n}}{2 h} + \alpha_ {0} ^ {n} u _ {0} ^ {n} = \eta_ {0} ^ {n} \\ \frac {3 u _ {J} ^ {n} - 4 u _ {J - 1} ^ {n} + u _ {J - 2} ^ {n}}{2 h} + \alpha_ {1} ^ {n} u _ {J} ^ {n} = \eta_ {1} ^ {n} \end{array} \right.
$$

### 5.2 一阶线性双曲型方程（组）的边界条件

#### 标量一阶双曲型方程：特征方向与传入边界

设有限区域内的对流方程：

$$
\frac {\partial \boldsymbol {u}}{\partial t} + \boldsymbol {a} \frac {\partial \boldsymbol {u}}{\partial x} = \boldsymbol {0}
$$

$$
\boldsymbol {x} \in (\mathbf {0}, l), t > 0
$$

初值：

$$
\boldsymbol {u} (\boldsymbol {x}, \boldsymbol {0}) = \boldsymbol {f} (\boldsymbol {x})
$$

$$
\boldsymbol {x} \in (0, l)
$$

如何正确给出边界条件？

若 $a > 0$ , 区域 $\{(x, t) | 0 < x < l, t \geq 0\}$ 

的右边界不能给出边界条件,

左边需给出边界条件：

$$
\boldsymbol {u} (\boldsymbol {0}, t) = \boldsymbol {g} (t), \quad t > 0 \text {且} \boldsymbol {g} (\boldsymbol {0}) = \boldsymbol {f} (\boldsymbol {0}) 。
$$

若 $a < 0$ , 完全相反,

需在右侧给出边界条件：

$$
\boldsymbol {u} (\boldsymbol {l}, t) = \boldsymbol {g} (\boldsymbol {t}), \quad \boldsymbol {t} > \mathbf {0} \text {且} \boldsymbol {g} (\mathbf {0}) = \boldsymbol {f} (\boldsymbol {l})
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/1cc501f3-2e11-4cac-83db-e365e24c70f1/b59eb22c0c662422c82c61f137f79c5284b40a08b4470185f9ef178fbfff0248.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/1cc501f3-2e11-4cac-83db-e365e24c70f1/62a20c0bb34c72800a784e712878886bde43646ec14e4804aea3ed3158524fe2.jpg)


#### 算例：内点格式所需的数值边界

例：对流方程 $\left\{ \begin{array}{l} \frac{\partial u}{\partial x} + a \frac{\partial u}{\partial t} = 0, (a < 0) 0 < x < 1, t > 0 \\ u(x,0) = g(x), 0 < x < 1 \\ u(1,t) = \varphi(x), t > 0 \end{array} \right.$ ，

考虑跳蛙格式: $u_{j}^{n+1} = u_{j}^{n-1} - a r\left(u_{j+1}^{n} - u_{j-1}^{n}\right)$ 

由于 $a < 0$ , 给的边界条件是 $x = 1$ 的边界, 然而在计算 $u_{1}^{2}$ 的时候需要用到 $u_{0}^{1}$ , 在计算 $u_{0}^{3}$ 的时候需要 $u_{0}^{2}$ 的值.

可是方程本身并没有给出在 $x = 0$ 的值. 需要我们增补这些值. 比如采用如下可行的方法:

##### 人工数值边界条件（格式）

(1) $u_{0}^{n + 1} = u_{0}^{n} + a\lambda \left(u_{1}^{n} - u_{0}^{n}\right)$ 

(2) $u_{0}^{n + 1} = u_{0}^{n - 1} + a\lambda \left[u_{1}^{n} - \frac{1}{2}\left(u_{0}^{n + 1} + u_{0}^{n - 1}\right)\right]$ 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/1cc501f3-2e11-4cac-83db-e365e24c70f1/1e2fa80348ba3ddae9b9b8921ec02dce5652539aa73aafc244fad3da3eb6d2ef.jpg)


解释 (1) : 根据特征线的性质有 $u(P) = u(Q)$ ,对于 $u(Q)$ 可以采用插值给出, 比如由 $x_{0}, x_{1}$ 进行线性插值得到, 即

$$
\begin{array}{l} u _ {0} ^ {n + 1} = u (Q) = u _ {1} ^ {n} + (1 - a r) u _ {0} ^ {n} \\ = u _ {0} ^ {n} + a r \left(u _ {1} ^ {n} - u _ {0} ^ {n}\right) \\ \end{array}
$$

实际上是迎风格式

注：采用插值法构造边界条件要用内插公式，使用外推方法往往是不行。即要用稳定的格式构造边界条件。

例如：下面的两个不可用的边界条件

用 $u_{1}, u_{2}$ 两点的值作线性插值, 外推得 $u_{0}$ 的值

$$
u _ {0} ^ {n} = 2 u _ {1} ^ {n} - u _ {2} ^ {n}
$$

再如（对边值不稳定）注、改进即为上页（2）

$$
u _ {0} ^ {n + 1} = u _ {0} ^ {n - 1} + 2 \lambda \left(u _ {1} ^ {n} - u _ {0} ^ {n}\right)
$$

#### 一阶双曲型方程组：特征块分解与边界定解

双曲型方程组：

$$
\frac {\partial \bar {\boldsymbol {u}}}{\partial t} + A \frac {\partial \bar {\boldsymbol {u}}}{\partial \boldsymbol {x}} = \boldsymbol {0}, \quad \boldsymbol {x} \in (\boldsymbol {0}, l), t > 0
$$

其中: $\bar{\pmb{u}} = \left(\pmb{u}_{1},\dots,\pmb{u}_{p}\right)^{T}, \quad A = \left(\pmb{a}_{ij}\right)_{p \times p} = \text {const}$ 

初始条件: $\overrightarrow{u}(x,0) = \overrightarrow{f}(x) = (f_1, \ldots, f_p)^T$ 

不失一般性，可设

$$
\boldsymbol {A} = \boldsymbol {\Lambda} = \left( \begin{array}{c c c} \boldsymbol {\Lambda} _ {1} & & \\ & \mathbf {0} & \\ & & \boldsymbol {\Lambda} _ {3} \end{array} \right) \quad \boldsymbol {\Lambda} _ {1} = \left( \begin{array}{c c c} \lambda_ {1} & & \\ & \lambda_ {2} & \\ & & \ddots \\ & & & \lambda_ {r} \end{array} \right) <   \mathbf {0}
$$

$$
\boldsymbol {\Lambda} _ {3} = \left( \begin{array}{c c c c} \lambda_ {m} & & & \\ & \lambda_ {m + 1} & & \\ & & \ddots & \\ & & & \lambda_ {p} \end{array} \right) > \mathbf {0}
$$

（各 $\lambda_j>0$；部分教材把同一事实写成「$- \Lambda_3$ 负定」，勿与特征值的符号混淆。）

0表示 $(m - r + 1) \times (m - r + 1)$ 的零矩阵。

这样 $\bar{\pmb{u}} = \left(\bar{\pmb{u}}_{I}, \bar{\pmb{u}}_{II}, \bar{\pmb{u}}_{III}\right)^{T}$ , 方程可分裂为:

$$
\frac {\partial \vec {u} _ {I}}{\partial t} + \Lambda_ {1} \frac {\partial \vec {u} _ {I}}{\partial x} = 0 \quad \frac {\partial \vec {u} _ {I I}}{\partial t} = 0
$$

$$
\frac {\partial \vec {u} _ {_ {I I I}}}{\partial t} + \varLambda_ {3} \frac {\partial \vec {u} _ {_ {I I I}}}{\partial x} = 0
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/1cc501f3-2e11-4cac-83db-e365e24c70f1/2212b2d7976e371358a6b3aab34c63f6f426ae1747573eaa97a8340ea2f2b8ff.jpg)

**法则（与特征值符号统一表述）**：对流方程 $\frac{\partial u}{\partial t}+\lambda\frac{\partial u}{\partial x}=0$ 的特征线满足 $\frac{dx}{dt}=\lambda$。**特征值为正**的分量，信息沿 $x$ 增大方向自**左边界**流入区域，须在 $x=0$ 处给定边界条件；**特征值为负**的分量，信息沿 $x$ 减小方向自**右边界**流入区域，须在 $x=l$ 处给定边界条件。下设分块 $\bar{u}_I$ 满足 $\partial_t \bar{u}_I + \Lambda_1 \partial_x \bar{u}_I = 0$（$\Lambda_1$ 的特征值均 $<0$），$\bar{u}_{III}$ 满足 $\partial_t \bar{u}_{III} + \Lambda_3 \partial_x \bar{u}_{III} = 0$（$\Lambda_3$ 的特征值均 $>0$）。阅读示意图时以「信息从哪一条边界进入计算区域」为准，不单用「向左/向右」以免与坐标轴正向混淆。

对上述分块：$\bar{\pmb{u}}_{I}$ 对应 $\Lambda_1<0$，信息从 $x=l$ 一侧进入区域，故在 $\pmb{x} = \pmb{l}$ 处给出边界条件

$$
\bar {u} _ {I} (l, t) = \bar {g} _ {I} (t), \quad \text {且}: \bar {g} _ {I} (\mathbf {0}) = \bar {u} _ {I} (l, \mathbf {0}) = \bar {f} _ {I} (l)
$$

$\bar{u}_{II}$ 满足的方程为 $\partial_t \bar{u}_{II}=0$，不涉及 $x$ 方向特征传播，不需要 $x=0,l$ 上的边界条件。

$\overline{\boldsymbol{u}}_{III}$ 对应 $\Lambda_3>0$，信息从 $x=0$ 一侧进入区域，故在 $\boldsymbol{x} = \mathbf{0}$ 处给出边界条件：

$$
\bar {\boldsymbol {u}} _ {I I I} (\mathbf {0}, t) = \bar {\boldsymbol {g}} _ {I I I} (t), \text {且}: \bar {\boldsymbol {g}} _ {I I I} (\mathbf {0}) = \bar {\boldsymbol {u}} _ {I I I} (\mathbf {0}, \mathbf {0}) = \bar {\boldsymbol {f}} _ {I I I} (\mathbf {0})
$$

否则错误地给出边界条件，问题会不适定。

例：考虑微分方程组（半无界问题）

$$
\frac {\partial \vec {u}}{\partial t} + A \frac {\partial \vec {u}}{\partial x} = 0, \quad x > 0, t > 0
$$

其中： $\vec{u} = (\mathbf{u}, \mathbf{v})$ ， $A = \begin{pmatrix} 0 & -1 \\ -1 & 0 \end{pmatrix}$ 

定解条件为:

$$
\left\{ \begin{array}{l} u (0, t) = 0 \\ u (x, 0) = \varphi (x) \\ v (x, 0) = g (x) \end{array} \right.
$$

思考: 不能给出 $v(x, t)$ 在 $x = 0$ 处的边界条件, 否则问题不适定。

##### 采用 Lax–Wendroff 格式

$$
\left\{ \begin{array}{l} u _ {j} ^ {n + 1} = u _ {j} ^ {n} + \frac {1}{2} \lambda (v _ {j + 1} ^ {n} - v _ {j - 1} ^ {n}) + \frac {1}{2} \lambda^ {2} (u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}) \\ v _ {j} ^ {n + 1} = v _ {j} ^ {n} + \frac {1}{2} \lambda (u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}) + \frac {1}{2} \lambda^ {2} (v _ {j + 1} ^ {n} - 2 v _ {j} ^ {n} + v _ {j - 1} ^ {n}) \end{array} \right.
$$

需要附加边界条件

$$
v _ {0} ^ {n + 1}
$$

##### 方法一、从特征形式出发

特征型: $\frac{\partial \vec{v}}{\partial t} + \Lambda \frac{\partial \vec{v}}{\partial x} = 0, x > 0, t > 0$ 

$$
\begin{array}{l} \vec {v} = S ^ {- 1} u = \left(\xi , \eta\right) ^ {T}, \quad \Lambda = S ^ {- 1} A S = \left[ \begin{array}{l l} 1 & 0 \\ 0 & - 1 \end{array} \right] \\ S = \frac {1}{2} \left[ \begin{array}{l l} 1 & 1 \\ - 1 & 1 \end{array} \right], \quad S ^ {- 1} = \frac {1}{2} \left[ \begin{array}{l l} 1 & - 1 \\ 1 & 1 \end{array} \right] \\ \end{array}
$$

计算 $v_{0}^{n + 1}$ , 就转化为计算 $\eta_{0}^{n + 1}$ 

采用迎风格式 $\eta_0^{n + 1} = \eta_0^n +\lambda (\eta_1^n -\eta_0^n)$ 

##### 方法二、从方程本身出发

已知边界条件 $u(0, t) = 0$ 

$\frac {\partial u (0 , t)}{\partial t} = 0, t > 0$ 

利用第一个方程: $\frac{\partial u(0, t)}{\partial t} = \frac{\partial v(0, t)}{\partial x} = 0, \quad t > 0$ 

(1)一阶差商 $\nu_{0}^{n + 1} = \nu_{1}^{n + 1}$ 

(2)二阶差商 $3 v_{0}^{n+1} - 4 v_{1}^{n+1} + v_{2}^{n+1} = 0$ 

利用第二个方程: $\frac{\partial v(0, t)}{\partial t} = \frac{\partial u(0, t)}{\partial x},\quad t > 0$ 

$$
v _ {0} ^ {n + 1} = v _ {0} ^ {n} + \lambda \left(u _ {1} ^ {n} - u _ {0} ^ {n}\right)
$$

## 第六节 二维问题

### 6.1 一阶双曲型方程

初值问题： $\left\{ \begin{array}{l} \frac{\partial \boldsymbol{u}}{\partial t} + \boldsymbol{a} \frac{\partial \boldsymbol{u}}{\partial x} + \boldsymbol{b} \frac{\partial \boldsymbol{u}}{\partial y} = \boldsymbol{0} \\ \boldsymbol{u}(x, y, 0) = \boldsymbol{u}_0(x, y), -\infty < x, y < +\infty \end{array} \right.$ 

其解为: $u(x, y, t) = u_{0}(x - at, y - bt)$ 

一般设 $\Delta x = \Delta y = h$ , 有:

(1) Lax-Friedrichs格式:

$$
\begin{array}{l} \frac {\boldsymbol {u} _ {j, m} ^ {n + 1} - \frac {1}{4} \left(\boldsymbol {u} _ {j, m + 1} ^ {n} + \boldsymbol {u} _ {j, m - 1} ^ {n} + \boldsymbol {u} _ {j + 1, m} ^ {n} + \boldsymbol {u} _ {j - 1, m} ^ {n}\right)}{\tau} \\ + a \frac {\boldsymbol {u} _ {j + 1 , m} ^ {n} - \boldsymbol {u} _ {j - 1 , m} ^ {n}}{2 h} + b \frac {\boldsymbol {u} _ {j , m + 1} ^ {n} - \boldsymbol {u} _ {j , m - 1} ^ {n}}{2 h} = \mathbf {0} \\ \end{array}
$$

此格式是一阶精度的。网格比 $\lambda = \tau / h$ 

下面讨论稳定性:

设 $\boldsymbol{u}_{j, m}^{n} = v^{n} e^{i\left(k_{1} j h + k_{2} m h\right)}$ 代入差分格式，有：

$$
v ^ {n + 1} = \left[ \frac {1}{2} \left(\cos k _ {1} h + \cos k _ {2} h\right) - i \lambda \left(a \sin k _ {1} h + b \sin k _ {2} h\right) \right] v ^ {n}
$$

$$
G \left(\lambda , k _ {1}, k _ {2}\right) = \frac {1}{2} \left(\cos k _ {1} h + \cos k _ {2} h\right) - i \lambda \left(a \sin k _ {1} h + b \sin k _ {2} h\right)
$$

$$
\begin{array}{l} \left| G \right| ^ {2} = \frac {1}{4} \left(\cos k _ {1} h + \cos k _ {2} h\right) ^ {2} + \lambda^ {2} \left(a \sin k _ {1} h + b \sin k _ {2} h\right) ^ {2} \\ = 1 - \left(\sin^ {2} k _ {1} h + \sin^ {2} k _ {2} h\right) \left[ \frac {1}{2} - \lambda^ {2} \left(a ^ {2} + b ^ {2}\right) \right] - \frac {1}{4} \left(\cos k _ {1} h - \cos k _ {2} h\right) ^ {2} \\ - \lambda^ {2} (a \sin k _ {1} h - b \sin k _ {2} h) ^ {2} \\ \end{array}
$$

$$
\leq 1 - \left(\sin^ {2} k _ {1} h + \sin^ {2} k _ {2} h\right) \left[ \frac {1}{2} - \lambda^ {2} \left(a ^ {2} + b ^ {2}\right) \right]
$$

如果 $\lambda^2 (a^2 + b^2) \leq \frac{1}{2}$ 

即 $\sqrt{\left(a^{2}+b^{2}\right)} \lambda \leq \frac{1}{\sqrt{2}}$ 

那么Von Neumann条件条件满足，格式稳定。

若 $\mathbf{b} = \mathbf{a}$ , 则条件为 $|\mathbf{a}| \lambda \leq \frac{1}{2}$ 

#### (2) Lax-Wendroff格式

设 $u(x, y, t)$ 是方程 $\frac{\partial u}{\partial t} + a \frac{\partial u}{\partial x} + b \frac{\partial u}{\partial y} = 0$ 的解，那么：

$$
\begin{array}{l} \frac {\partial^ {2} \boldsymbol {u}}{\partial t ^ {2}} = \frac {\partial}{\partial t} \left(- \boldsymbol {a} \frac {\partial \boldsymbol {u}}{\partial x} - \boldsymbol {b} \frac {\partial \boldsymbol {u}}{\partial y}\right) \\ = - \boldsymbol {a} \frac {\partial}{\partial \boldsymbol {x}} \left(\frac {\partial \boldsymbol {u}}{\partial t}\right) - \boldsymbol {b} \frac {\partial}{\partial \boldsymbol {y}} \left(\frac {\partial \boldsymbol {u}}{\partial t}\right) \\ = a ^ {2} \frac {\partial^ {2} u}{\partial x ^ {2}} + 2 a b \frac {\partial^ {2} u}{\partial x \partial y} + b ^ {2} \frac {\partial^ {2} u}{\partial y ^ {2}} \\ \end{array}
$$

又由Taylor展开，有

$$
\begin{array}{l} u \left(x _ {j}, y _ {m}, t + \tau\right) = \left[ I - \tau \left(a \frac {\partial}{\partial x} + b \frac {\partial}{\partial y}\right) + \frac {\tau^ {2}}{2} \left(a ^ {2} \frac {\partial^ {2}}{\partial x ^ {2}} + 2 a b \frac {\partial^ {2}}{\partial x \partial y} + b ^ {2} \frac {\partial^ {2}}{\partial y ^ {2}}\right) \right] \\ u (x _ {j}, y _ {m}, t _ {n}) + O (\tau^ {3}) \\ \end{array}
$$

用二阶中心差商代替偏导数，有Lax-Wendroff格式：

$$
\begin{array}{l} \boldsymbol {u} _ {j, m} ^ {n + 1} = \boldsymbol {u} _ {j, m} ^ {n} - \frac {\lambda}{2} \left(\boldsymbol {a} \Delta_ {0} ^ {x} \boldsymbol {u} _ {j, m} ^ {n} + \boldsymbol {b} \Delta_ {0} ^ {y} \boldsymbol {u} _ {j, m} ^ {n}\right) \\ + \frac {1}{2} \lambda^ {2} \left(a ^ {2} \Delta_ {+} ^ {x} \Delta_ {-} ^ {x} u _ {j, m} ^ {n} + \frac {1}{2} a b \Delta_ {0} ^ {x} \Delta_ {0} ^ {y} u _ {j, m} ^ {n} + b ^ {2} \Delta_ {+} ^ {y} \Delta_ {-} ^ {y} u _ {j, m} ^ {n}\right) \\ \end{array}
$$

其中: $\Delta_{0}^{x} u_{j, m}^{n} = u_{j+1, m}^{n} - u_{j-1, m}^{n}, \quad \Delta_{+}^{x} u_{j, m}^{n} = u_{j+1, m}^{n} - u_{j, m}^{n}$ 

$$
\Delta_ {-} ^ {x} \boldsymbol {u} _ {j, m} ^ {n} = \boldsymbol {u} _ {j, m} ^ {n} - \boldsymbol {u} _ {j - 1, m} ^ {n}, \quad \Delta_ {0} ^ {y}, \Delta_ {+} ^ {y}, \Delta_ {-} ^ {y} \text {类 似 。}
$$

这是二阶格式, 稳定性条件 $\sqrt{\left(a^{2} + b^{2}\right)} \lambda \leq \frac{1}{\sqrt{2}}$ 

#### (3) 分数步长法（局部一维格式）

为放宽稳定性条件而引入的技巧。

方法是：

第一步由 $\mathbf{x}$ 方向的差分把 $t_{n}$ 推进到 $t_{n} + \frac{\tau}{2}$ ;

第二步由 $\mathbf{y}$ 方向的差分把 $t_{n} + \frac{\tau}{2}$ 推进到 $t_{n+1}$ 。一般形式：

$$
\left\{ \begin{array}{l} u _ {j, m} ^ {n + \frac {1}{2}} = u _ {j, m} ^ {n} + \tau D _ {1} u _ {j, m} ^ {n} \\ u _ {j, m} ^ {n + 1} = u _ {j, m} ^ {n + \frac {1}{2}} + \tau D _ {2} u _ {j, m} ^ {n + \frac {1}{2}} \end{array} \right.
$$

例如:

以Lax-Wendroff格式来完成二步法

$$
\begin{array}{l} D _ {1} = - \boldsymbol {a} \cdot \frac {1}{2 h} \Delta_ {0} ^ {x} + \frac {\tau}{2} \boldsymbol {a} ^ {2} \cdot \frac {1}{h ^ {2}} \Delta_ {+} ^ {x} \Delta_ {-} ^ {x} \\ D _ {2} = - \boldsymbol {b} \cdot \frac {1}{2 h} \Delta_ {0} ^ {y} + \frac {\tau}{2} \boldsymbol {b} ^ {2} \cdot \frac {1}{h ^ {2}} \Delta_ {+} ^ {y} \Delta_ {-} ^ {y} \\ \end{array}
$$

这样计算的稳定性条件是：

$$
| a | \lambda <   1, | b | \lambda <   1
$$

#### (4) 隐式格式

全隐格式 精度: $O(\tau + h)$ ; 无条件稳定

$$
\frac {u _ {j , m} ^ {n + 1} - u _ {j , m} ^ {n}}{\tau} + a \frac {u _ {j + 1 , m} ^ {n + 1} - u _ {j - 1 , m} ^ {n + 1}}{2 h} + b \frac {u _ {j , m + 1} ^ {n + 1} - u _ {j , m - 1} ^ {n + 1}}{2 h}
$$

Crank-Nicolson格式 精度: $O\left(\tau^{2} + h^{2}\right)$ ; 无条件稳定

$$
\begin{array}{l} \frac {u _ {j , m} ^ {n + 1} - u _ {j , m} ^ {n}}{\tau} + \frac {1}{2} \{a \frac {u _ {j + 1 , m} ^ {n + 1} - u _ {j - 1 , m} ^ {n + 1}}{2 h} + b \frac {u _ {j , m + 1} ^ {n + 1} - u _ {j , m - 1} ^ {n + 1}}{2 h} \\ + a \frac {u _ {j + 1 , m} ^ {n} - u _ {j - 1 , m} ^ {n}}{2 h} + b \frac {u _ {j , m + 1} ^ {n} - u _ {j , m - 1} ^ {n}}{2 h} \} = 0 \\ \end{array}
$$

#### (5) ADI（Alternate Direction Implicit）交替方向隐式格式

##### ADI-1：局部一维格式

$$
\left\{ \begin{array}{l} {\left(1 + \frac {a}{2} \lambda \Delta_ {0} ^ {x}\right) u _ {j m} ^ {n + \frac {1}{2}} = u _ {j m} ^ {n}} \\ {\left(1 + \frac {b}{2} \lambda \Delta_ {0} ^ {y}\right) u _ {j m} ^ {n + 1} = u _ {j m} ^ {n + \frac {1}{2}}} \end{array} \right.
$$

##### 二维 Beam–Warming 格式

$$
\left(1 + \frac {a \lambda}{4} \Delta_ {0} ^ {x}\right) \left(1 + \frac {b \lambda}{4} \Delta_ {0} ^ {y}\right) u _ {j m} ^ {n + 1} = \left(1 - \frac {a \lambda}{4} \Delta_ {0} ^ {x}\right) \left(1 - \frac {b \lambda}{4} \Delta_ {0} ^ {y}\right) u _ {j m} ^ {n}
$$

##### ADI-2

$$
\left\{ \begin{array}{l} \left(1 + \frac {a \lambda}{4} \Delta_ {0} ^ {x}\right) u _ {j m} ^ {*} = \left(1 - \frac {a \lambda}{4} \Delta_ {0} ^ {x}\right) \left(1 - \frac {b \lambda}{4} \Delta_ {0} ^ {y}\right) u _ {j m} ^ {n} \\ \left(1 + \frac {b \lambda}{4} \Delta_ {0} ^ {y}\right) u _ {j m} ^ {n + 1} = u _ {j m} ^ {*} \end{array} \right.
$$

### 6.2 一阶双曲型方程组

设方程: $\frac{\partial \bar{u}}{\partial t} + A \frac{\partial \bar{u}}{\partial x} + B \frac{\partial \bar{u}}{\partial y} = 0$ 

其中 $\bar{\pmb{u}} = \left(\pmb{u}_{1},\dots ,\pmb{u}_{p}\right)^{T}, A, B$ 为实的 $p^{*} p$ 矩阵

如果对所有的 $\alpha, \beta, \alpha + \beta = 1$ , 有非奇异的矩阵 $S$ 使 $\Lambda = S\left(\alpha A + \beta B\right) S^{-1}$ 为实对角矩阵。

则称方程组是双曲型方程组。

如果A,B为实对称阵, 则方程组是双曲型方程组,也称为对称双曲型方程组。

下面以Lax-Wendroff格式为例，讨论差分方程：利用多元Taylor展开，有

$$
\begin{array}{l} \vec {u} (x, y, t + \tau) \\ = \vec {u} (x, y, t) + \tau \frac {\partial \vec {u}}{\partial t} (x, y, t) \\ + \frac {\tau^ {2}}{2} \frac {\partial^ {2} \vec {u}}{\partial t ^ {2}} (x, y, t) + O (\tau^ {3}) \\ = \vec {u} (x, y, t) - \tau \left(A \frac {\partial \vec {u}}{\partial x} + B \frac {\partial \vec {u}}{\partial y}\right) \\ + \frac {\tau^ {2}}{2} \left[ A ^ {2} \frac {\partial^ {2} \vec {u}}{\partial x ^ {2}} + (A B + B A) \frac {\partial^ {2} \vec {u}}{\partial x \partial y} + B ^ {2} \frac {\partial^ {2} \vec {u}}{\partial y ^ {2}} \right] + O \left(\tau^ {3}\right) \\ \end{array}
$$

故有:

$$
\begin{array}{l} \vec {\mathcal {U}} _ {j, m} ^ {n + 1} = L _ {h} \vec {\mathcal {U}} _ {j, m} ^ {n} \\ = [ I - \frac {1}{2} \lambda (A \Delta_ {0} ^ {x} + B \Delta_ {0} ^ {y}) + \frac {1}{2} \lambda^ {2} (A ^ {2} \Delta_ {+} ^ {x} \Delta_ {-} ^ {x} + B ^ {2} \Delta_ {+} ^ {y} \Delta_ {-} ^ {y}) \\ + \frac {1}{8} \lambda^ {2} \left(A B + B A\right) \Delta_ {0} ^ {x} \Delta_ {0} ^ {y} ] \vec {u} _ {j, m} ^ {n} \\ \end{array}
$$

利用Fourier方法可讨论上式的稳定性:

$$
\text {令} \bar {u} _ {j, m} ^ {n} = \bar {\nu} ^ {n} e ^ {i \left(k _ {1} j h + k _ {2} m h\right)}
$$

可得增长矩阵:

$$
\begin{array}{l} G \big (\lambda , k _ {1}, k _ {2} \big) \\ = I + i \lambda \left(A \sin k _ {1} h + B \sin k _ {2} h\right) + \lambda^ {2} \left[ A ^ {2} \left(\cos k _ {1} h - 1\right) \right. \\ + B ^ {2} \left(\cos k _ {2} h - 1\right) ] - \frac {1}{2} \lambda^ {2} \left(A B + B A\right) \sin k _ {1} h \sin k _ {2} h \\ \end{array}
$$

如果A,B是对称阵,

可以证明Lax-Wendroff格式的稳定性条件是:

$$
\lambda \rho (A) \leq \frac {1}{2 \sqrt {2}}, \lambda \rho (B) \leq \frac {1}{2 \sqrt {2}}
$$

为放宽稳定性条件，有许多改进技术。

如: Strang方法

### 6.3 高阶双曲型方程

考虑波动方程

$$
\frac {\partial^ {2} u}{\partial t ^ {2}} - a ^ {2} \left(\frac {\partial^ {2} u}{\partial x ^ {2}} + \frac {\partial^ {2} u}{\partial y ^ {2}}\right) = 0, a > 0
$$

由于显隐格式均不适用于多个空间变量的情况,下面构造交替方向隐格式, 这种格式有良好的稳定性以及得到的差分格式方程组的系数矩阵是三对角的, 用追赶法求解。

#### 6.3.1 ADI 格式（1）

$$
\begin{array}{l} \left[ \frac {u _ {j} ^ {n + 1 / 2} - 2 u _ {j} ^ {n} + u _ {j} ^ {n - 1}}{\tau^ {2}} = \frac {a ^ {2}}{h ^ {2}} \left[ \delta_ {x} ^ {2} \left(\theta u _ {j} ^ {n + 1 / 2} + (1 - 2 \theta) u _ {j} ^ {n} + \theta u _ {j} ^ {n - 1}\right) \right. \right. \\ + \delta_ {y} ^ {2} \left(2 \theta u _ {j} ^ {n - 1} + (1 - 2 \theta) u _ {j} ^ {n}\right) ] \\ \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n + 1 / 2}}{\tau / 2} = \frac {a ^ {2}}{h ^ {2}} \delta_ {y} ^ {2} (\theta u _ {j} ^ {n + 1} - \theta u _ {j} ^ {n - 1}) \\ \end{array}
$$

$\theta \geq \frac{1}{4}$ 是无条件稳定的。

#### 6.3.2 ADI 格式（2）

$$
\begin{array}{l} \left[ \frac {u _ {j} ^ {n + 1 / 2} - 2 u _ {j} ^ {n} + u _ {j} ^ {n - 1}}{\tau^ {2}} \right. \\ = \frac {a ^ {2}}{h ^ {2}} \left[ \delta_ {x} ^ {2} \left(\theta u _ {j} ^ {n + 1 / 2} + (1 - 2 \theta) u _ {j} ^ {n} + \theta u _ {j} ^ {n - 1}\right) + \delta_ {y} ^ {2} u _ {j} ^ {n} \right] \\ \frac {\delta_ {t} ^ {2} u _ {j} ^ {n}}{\tau^ {2}} = \frac {a ^ {2}}{h ^ {2}} \left[ \delta_ {x} ^ {2} \left(\theta \left(u _ {j} ^ {n + 1 / 2} + u _ {j} ^ {n - 1}\right) + (1 - 2 \theta) u _ {j} ^ {n}\right) \right. \\ + \delta_ {y} ^ {2} \left(\theta \left(u _ {j} ^ {n - 1} + u _ {j} ^ {n + 1}\right) + \left(1 - 2 \theta\right) \left. u _ {j} ^ {n}\right) \right] \\ \end{array}
$$

$\theta \geq \frac{1}{4}$ 是无条件稳定的。

#### 6.3.3 高维波动方程的等价一阶双曲方程组

引入变量化为2维一阶方程组

引入: $v = \frac{\partial u}{\partial t}, w = a \frac{\partial u}{\partial x}, z = a \frac{\partial u}{\partial y}$ 

令 $U = (v, w, z)$ ，则有

$$
\frac {\partial U}{\partial t} + A \frac {\partial U}{\partial x} + B \frac {\partial U}{\partial y} = 0
$$

$$
A = \left( \begin{array}{c c c} 0 & - a & 0 \\ - a & 0 & 0 \\ 0 & 0 & 0 \end{array} \right)
$$

$$
B = \left( \begin{array}{c c c} 0 & 0 & - a \\ 0 & 0 & 0 \\ - a & 0 & 0 \end{array} \right)
$$

一维的双曲方程（组）的所有差分格式都可以推广到高维方程（组）。下面以LW为例。

$$
\frac {\partial U}{\partial t} = - A \frac {\partial U}{\partial x} - B \frac {\partial U}{\partial y}
$$

$$
\begin{array}{l} \frac {\partial^ {2} U}{\partial t ^ {2}} = - A \frac {\partial}{\partial t} \frac {\partial U}{\partial x} - B \frac {\partial}{\partial t} \frac {\partial U}{\partial y} \\ = A ^ {2} \frac {\partial^ {2} U}{\partial x ^ {2}} + (A B + B A) \frac {\partial^ {2} U}{\partial x \partial y} + B ^ {2} \frac {\partial^ {2} U}{\partial y ^ {2}} \\ \end{array}
$$

#### Lax–Wendroff 格式

$$
\begin{array}{l} L _ {h} = I - \frac {\lambda}{2} [ A (T _ {x} - T _ {x} ^ {- 1}) + B (T _ {y} - T _ {y} ^ {- 1}) ] \\ + \frac {1}{2} \lambda^ {2} [ A ^ {2} (T _ {x} - 2 I + T _ {x} ^ {- 1}) + B ^ {2} (T _ {y} - 2 I + T _ {y} ^ {- 1}) ] \\ + \frac {1}{8} \lambda^ {2} [ (A B + B A) (T _ {x} - T _ {x} ^ {- 1}) (T _ {y} - T _ {y} ^ {- 1}) ] \\ \end{array}
$$

其中

$$
\begin{array}{l} T _ {x} U _ {i j} = U _ {i + 1, j}, T _ {x} ^ {- 1} U _ {i j} = U _ {i - 1, j} \\ T _ {y} U _ {i j} = U _ {i, j + 1}, T _ {y} ^ {- 1} U _ {i j} = U _ {i, j - 1} \\ \end{array}
$$

#### 6.3.4 Strang 算子分裂格式

$$
U _ {i j} ^ {k + 1} = L _ {h} (\tau) U _ {i j} ^ {k}
$$

$$
L _ {h} (\tau) = L _ {h} ^ {x} (\frac {\tau}{2}) L _ {h} ^ {y} (\tau) L _ {h} ^ {x} (\frac {\tau}{2})
$$

$$
L _ {h} ^ {x} (\tau) = I - \frac {\lambda}{2} A (T _ {x} - T _ {x} ^ {- 1}) + \frac {1}{2} \lambda^ {2} A ^ {2} (T _ {x} - 2 I + T _ {x} ^ {- 1})
$$

$$
L _ {h} ^ {y} (\tau) = I - \frac {\lambda}{2} B (T _ {y} - T _ {y} ^ {- 1}) + \frac {1}{2} \lambda^ {2} B ^ {2} (T _ {y} - 2 I + T _ {y} ^ {- 1})
$$

此格式是等价方程组的差分格式且是二阶精度，因为

$$
L ^ {x} (\tau) = I - \frac {\tau}{2} A \frac {\partial}{\partial x} + \frac {1}{2} \tau^ {2} A ^ {2} \frac {\partial^ {2}}{\partial x ^ {2}}
$$

$$
L ^ {y} (\tau) = I - \frac {\tau}{2} B \frac {\partial}{\partial y} + \frac {1}{2} \tau^ {2} B ^ {2} \frac {\partial^ {2}}{\partial y ^ {2}}
$$

稳定性条件：各一步 Lax–Wendroff 型一维算子仍须满足一维 CFL 条件；Strang 分裂中 $x$ 方向虽分两个半步 $\tau/2$，**并不能**把 $x$ 方向的允许 Courant 数放大到 2。标准取法为

$$
\lambda \rho (A) \leq 1, \quad \lambda \rho (B) \leq 1
$$
