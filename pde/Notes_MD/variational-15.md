# 第七章 数学物理方程的变分原理

## 变分问题介绍

## 1.1 古典变分问题

**例 1.1（最速降线问题）** 这是 Bernoulli 于 1696 年提出的问题。如图 7.1，设点 $A(0,0)$ 与 $B(x_1,y_1)$（取 $y$ 轴垂直向下）不在同一条与 $y$ 轴平行的直线上。一质点受重力从 $A$ 沿光滑曲线路径自由下滑至 $B$，求下降最快的路径；不计摩擦。

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/42856561caafd102ef68a18de98143b951ae3dbadf55d5ad8dac9d1bc6944906.jpg)

“下降最快”即所用时间最短。设从 $A$ 到 $B$ 的任一光滑曲线为

$$
l:\ y=y(x),\quad 0\le x\le x_1.
$$

从 $A$ 到曲线上任一点 $P(x,y)$，弧长元为 $\mathrm{d}s$。速率 $v=\dfrac{\mathrm{d}s}{\mathrm{d}t}$，动能为 $\tfrac12 mv^2$。质点相对 $A$ 下降高度 $y$ 时，势能减少 $mgy$。由机械能守恒

$$
\frac{1}{2}m\left(\frac{\mathrm{d}s}{\mathrm{d}t}\right)^2 = mgy.
$$

故 $v=\dfrac{\mathrm{d}s}{\mathrm{d}t}=\sqrt{2gy}$。又 $\mathrm{d}s=\sqrt{1+\left(y'\right)^{2}}\,\mathrm{d}x$，因此

$$
\sqrt{1+\left(y'\right)^{2}}\,\frac{\mathrm{d}x}{\mathrm{d}t}=\sqrt{2gy},
\qquad
\mathrm{d}t=\sqrt{\frac{1+\left(y'\right)^{2}}{2gy}}\,\mathrm{d}x.
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/8b811ef1c7428d8500c86951dd3f80b7889b72ce065fb4268a5050770cc612e2.jpg)

从 $A$ 沿 $l$ 下滑至 $B$ 的总时间为

$$
T=\int_0^{x_1}\sqrt{\frac{1+\left(y'\right)^{2}}{2gy}}\,\mathrm{d}x. \tag{1.1}
$$

$T$ 依赖于函数 $y=y(x)$ 的整体选取。

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/207d87f85b913fc374175b73d60aa0b41c56ba6ee750d50d6eefa4edfec9404f.jpg)

当 $y$ 在某一函数集合 $K$ 中取定后，（1.1）给出一个确定的实数 $T$。称 $T$ 为 $K$ 上的**泛函**，记

$$
T(y)=\int_0^{x_1}\sqrt{\frac{1+\left(y'\right)^{2}}{2gy}}\,\mathrm{d}x.
$$

即 $T=T(y)$ 定义在 $K$ 上。例如

$$
K=\{y\mid y\in C^1[0,x_1],\ y(0)=0,\ y(x_1)=y_1\}.
$$

最速降线问题：在 $K$ 上求 $y_0$ 使 $T(y_0)=\min\limits_{y\in K}T(y)$。

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/a576be9eefdf8d6181a7c470cfde47945b3aef996415a09349f6d0ede7e22c48.jpg)

例1.2（最小曲面问题）如图7.2：在 $xy$ 平面上给定有界开区域 $\Omega$ 与其边界 $\partial\Omega$，在 $\partial\Omega$ 上给定 $u\big|_{\partial\Omega}=\varphi(x,y)$ ，其中 $\varphi$ 为已知函数。

这样就给出了三维空间中的一条封闭曲线C，最小曲面问题就是求张紧在曲线C的曲面中，其面积最小的曲面。

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/a5f9adb1c76b7816d8f67b1239f919369d32331917a13158aedc3b1d5fea3490.jpg)

记曲面 $u=u(x,y)$，$(x,y)\in\overline\Omega$。

对应 $u$ 的曲面面积为

$$
S (u) = \iint_ {\Omega} \sqrt {1 + \left(\frac {\partial u}{\partial x}\right) ^ {2} + \left(\frac {\partial u}{\partial y}\right) ^ {2}} d x d y \tag{1.2}
$$

例1.3（等周问题）在长度为一定的所有平面光滑封闭曲线中求所围面积为最大的曲线（图7.3）。

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/8d2d6e2a6556db3ca2369dceccb893542bba3ed99b2c6034e5a4c00326c0e251.jpg)

设曲线以弧长参数表示为 $x=x(s),\ y=y(s)$，$s_1\le s\le s_2$。**闭合条件**为

$$
x(s_1)=x(s_2),\quad y(s_1)=y(s_2).
$$

**周长约束**（给定弧长 $l$）为

$$
\int_{s_1}^{s_2}\sqrt{\left(\frac{\mathrm{d}x}{\mathrm{d}s}\right)^2+\left(\frac{\mathrm{d}y}{\mathrm{d}s}\right)^2}\,\mathrm{d}s=l.
$$

所围面积（Green 公式）为

$$
A(x,y)=\frac12\int_{s_1}^{s_2}\left(x\frac{\mathrm{d}y}{\mathrm{d}s}-y\frac{\mathrm{d}x}{\mathrm{d}s}\right)\mathrm{d}s.
$$

等周问题：在满足上述约束的 $(x(s),y(s))$ 集合上求 $A(x,y)$ 的极大。

---

**变分问题**：在给定的函数（或曲线）集合 $K$ 上，求某泛函的极小或极大。

- **例 1.1**：$K$ 为一元函数集合；最小化 $T(y)$。
- **例 1.2**：$K$ 为边界上取定 $\varphi$ 的多元函数集合；最小化（1.2）。
- **例 1.3**：$K$ 为满足周长与闭合一族参数曲线；最大化 $A$。

$$
\left\{ \begin{array}{l} {\text{求 } u _ {0} \in K, \text{使得}} \\ {J \left(u _ {0}\right) \leq J \left(u\right),\ \forall u \in K.} \end{array} \right. \tag{1.3}
$$

## 1.2 变分问题解的必要条件

为了推导变分问题解的必要条件，我们首先定义端点取值为零的连续可微函数空间：

$$
C _ {0} ^ {1} [ a, b ] = \left\{v \mid v \in C ^ {1} [ a, b ], v (a) = 0, v (b) = 0 \right\}.
$$

**引理 1.1（变分法基本引理）** 设 $u\in C[a,b]$，且

$$
\int_a^b u(x)v(x)\,\mathrm{d}x=0,\quad \forall v\in C_0^1[a,b],
$$

则 $u(x)\equiv 0$（$a\le x\le b$）。

**证** 若存在 $\xi\in(a,b)$ 使 $u(\xi)\neq 0$，不妨设 $u(\xi)>0$。由 $u$ 的连续性，存在 $(\xi_1,\xi_2)\subset(a,b)$，使在 $[\xi_1,\xi_2]$ 上 $u(x)>0$。取

$$
v(x)=\begin{cases}
0, & a\le x\le \xi_1,\\
(x-\xi_1)^2(x-\xi_2)^2, & \xi_1<x<\xi_2,\\
0, & \xi_2\le x\le b
\end{cases} \tag{1.4}
$$

则 $v\in C_0^1[a,b]$，且

$$
\int_a^b u(x)v(x)\,\mathrm{d}x=\int_{\xi_1}^{\xi_2} u(x)(x-\xi_1)^2(x-\xi_2)^2\,\mathrm{d}x>0,
$$

与假设矛盾。故 $u\equiv 0$。

若将条件改为 $\forall v\in C_0^n[a,b]$，结论仍成立，其中

$$
C _ {0} ^ {n} [ a, b ] = \left\{v \mid v \in C ^ {n} [ a, b ], v^ {(k)} (a) = v^ {(k)} (b) = 0, k = 0, 1, \dots , n - 1 \right\}.
$$

证明时只需将（1.4）中多项式的幂次提高至 $n+1$（即取 $v(x) = (x-\xi_1)^{n+1}(x-\xi_2)^{n+1}$），即可保证函数在边界处直到 $n$ 阶导数皆为零，从而满足 $v \in C_0^n[a,b]$。

讨论最简单的变分问题，推导其解满足的必要条件。函数集合

$$
K = \left\{y \mid y \in C ^ {1} [ a, b ], y (a) = y _ {a}, y (b) = y _ {b} \right\} \tag{1.5}
$$

考虑与 $y$ 及 $y'$ 有关的泛函

$$
J (y) = \int_ {a} ^ {b} F \left(x, y, y ^ {\prime}\right) d x \tag{1.6}
$$

其中 $F$ 对各自变量的偏导数均连续。

**变分问题（Ⅰ）** 求 $y\in K$，使得 $J(y)\le J(w)$，$\forall w\in K$。

设 $y$ 是变分问题（Ⅰ）的解，则 $y\in K$。对一切 $\eta\in C_0^1[a,b]$，因 $\eta(a)=\eta(b)=0$，有 $y+\alpha\eta\in K$，其中 $\alpha\in\mathbb{R}$ 任意。

因 $y$ 使 $J(y)$ 达到极小，故

$$
J (y) \leq J (y + \alpha \eta), \quad \forall \alpha \in \mathbb{R}. \tag{1.7}
$$

把 $J(y+\alpha\eta)$ 看成 $\alpha$ 的一元函数，记

$$
\varphi (\alpha) = J (y + \alpha \eta) = \int_ {a} ^ {b} F (x, y (x) + \alpha \eta (x), y ^ {\prime} (x) + \alpha \eta^ {\prime} (x)) d x.
$$

$\varphi(\alpha)$ 在 $\alpha=0$ 附近的导数给出泛函的一阶变分；当 $y\in C^1[a,b]$ 时，

$$
\left. \frac{d}{d\alpha} J(y+\alpha\eta) \right|_{\alpha=0}
= \frac {d \varphi}{d \alpha} \bigg | _ {\alpha = 0} = \int_ {a} ^ {b} \left(\frac {\partial F (x , y , y ^ {\prime})}{\partial y} \eta + \frac {\partial F (x , y , y ^ {\prime})}{\partial y ^ {\prime}} \eta^ {\prime}\right) d x.
$$

由（1.7），$\alpha=0$ 是 $\varphi(\alpha)$ 的极小点，故 $\varphi'(0)=0$，从而

$$
\int_ {a} ^ {b} \left(\frac {\partial F}{\partial y} \eta + \frac {\partial F}{\partial y ^ {\prime}} \eta^ {\prime}\right) d x = 0, \quad \forall \eta \in C _ {0} ^ {1} [ a, b ].
$$

第二项分部积分，得

$$
\int_ {a} ^ {b} \left[ \frac {\partial F}{\partial y} - \frac {d}{d x} \left(\frac {\partial F}{\partial y ^ {\prime}}\right) \right] \eta\, d x = 0, \quad \forall \eta \in C _ {0} ^ {1} [ a, b ]. \tag{1.8}
$$

根据引理 1.1，得到 **Euler 方程**

$$
\frac {\partial F}{\partial y} - \frac {d}{d x} \left(\frac {\partial F}{\partial y ^ {\prime}}\right) = 0. \tag{1.9}
$$

**例 1.4** 设 $K$ 仍如（1.5）式所示，且

$$
J (y) = \frac {1}{2} \int_ {a} ^ {b} \left[ p (x) \left(\frac {d y}{d x}\right) ^ {2} + q (x) y ^ {2} - 2 f (x) y \right] d x,
$$

求 $J(y)$ 在 $K$ 上的极值。若极值在 $y\in C^2[a,b]$ 处取得，则 $y$ 满足

$$
-\frac{d}{dx}\left(p(x)\frac{dy}{dx}\right)+q(x)y=f(x),
$$

连同 $y(a)=y_a,\ y(b)=y_b$ 构成定解问题。

在例 1.4 中，$F$ 关于 $y$ 与 $y'$ 的头两项为平方项，$J(y)$ 为**二次泛函**。

$$
K _ {1} = \left\{y \mid y \in C ^ {1} [ a, b ], \quad y (a) = y _ {a} \right\} \tag{1.10}
$$

同样求泛函（1.6）的极小值。

**变分问题（Ⅱ）** 求 $y\in K_1$，使得 $J(y)\le J(w)$，$\forall w\in K_1$。

在必要条件推导中，只要求

$$
\eta \in \left\{\nu \mid \nu \in C ^ {1} [ a, b ], \quad \nu (a) = 0 \right\}. \tag{1.11}
$$

此时

$$
\delta J = \left. \frac {d \varphi}{d \alpha} \right| _ {\alpha = 0} = \int_ {a} ^ {b} \left[ \frac {\partial F}{\partial y} - \frac {d}{d x} \left(\frac {\partial F}{\partial y ^ {\prime}}\right) \right] \eta\, d x + \left[ \eta \frac {\partial F}{\partial y ^ {\prime}} \right] _ {x = b} = 0 \tag{1.12}
$$

对一切满足（1.11）的 $\eta$ 成立。同上仍得（1.8）与（1.9）；代回（1.12）得

$$
\left[ \eta \frac {\partial F}{\partial y ^ {\prime}} \right] _ {x = b} = 0, \quad \forall \eta \in \left\{v \mid v \in C ^ {1} [ a, b ], \quad v (a) = 0 \right\}.
$$

由 $\eta$ 的任意性：

$$
\left. \frac {\partial F}{\partial y ^ {\prime}} \right| _ {x = b} = 0. \tag{1.13}
$$

（1.13）是在右端点 $x=b$ 上的**自然边界条件**。

**变分问题（Ⅱ）** 的解 $y$ 满足

$$
\left\{ \begin{array}{l} \dfrac {\partial F}{\partial y} - \dfrac {d}{d x} \left(\dfrac {\partial F}{\partial y ^ {\prime}}\right) = 0, \\ y (a) = y _ {a}, \quad \left. \dfrac{\partial F}{\partial y ^ {\prime}} \right| _ {x = b} = 0. \end{array} \right.
$$

将变分问题（Ⅱ）用于例 1.4，（1.13）化为 $\bigl[p(x)\frac{dy}{dx}\bigr]_{x=b}=0$。若 $p(x)\ge p_0>0$（$a\le x\le b$），在 $x=b$ 的自然边界条件即 $\left.\frac{dy}{dx}\right|_{x=b}=0$（第二类边界条件）。

考虑依赖于多个函数的泛函

$$
J (y _ {1}, \dots, y _ {n}) = \int_ {a} ^ {b} F (x, y _ {1}, \dots, y _ {n}; y _ {1} ^ {\prime}, \dots, y _ {n} ^ {\prime})\, d x,
$$

或依赖于高阶导数的泛函，例如

$$
J (y) = \int_ {a} ^ {b} F \left(x, y, y ^ {\prime}, y ^ {\prime \prime}\right) d x,
$$

也可以考虑依赖于多元函数的泛函，例如

$$
J (u) = \iint_ {\Omega} F (x, y, u, u _ {x}, u _ {y}) d x d y.
$$

可以分别讨论其极值的必要条件，即 Euler 方程。

## 1.3 $\mathbb{R}^n$ 中的变分问题

为了与无穷维（函数空间上的）泛函极值问题作对照，并为后续有限元等数值离散中常用的代数化思路做铺垫，我们先回顾有限维空间 $\mathbb{R}^n$ 中的变分与极值问题。

设向量 $x=(x_1,\dots,x_n)^T\in\mathbb{R}^n$，矩阵 $A=(a_{ij})\in\mathbb{R}^{n\times n}$。设 $\det A\neq 0$，则方程组 $Ax=b$ 有惟一解，其中已知向量 $b\in\mathbb{R}^n$。

**命题 1.1** 向量 $x \in \mathbb{R}^n$ 满足方程 $Ax - b = 0$。

**命题 1.2** $(Ax-b,y)=0$，$\forall y\in\mathbb{R}^n$。

若命题 1.1 成立，则命题 1.2 成立。若命题 1.2 成立，令 $y=e_i$（$i=1,\dots,n$）为各坐标轴上的单位向量，则可推出 $Ax-b=0$。命题 1.1 与命题 1.2 等价。

$$
J (x) = \frac {1}{2} \sum_ {i, j = 1} ^ {n} a _ {i j} x _ {i} x _ {j} - \sum_ {i = 1} ^ {n} b _ {i} x _ {i}
$$

矩阵形式为 $J(x)=\frac{1}{2}(Ax,x)-(b,x)$。设 $A$ 对称，则 $(Ax,y)=(Ay,x)$，$\forall x,y\in\mathbb{R}^n$。

$$
\begin{array}{l} J (x + \alpha y) = \frac {1}{2} (A (x + \alpha y), x + \alpha y) - (b, x + \alpha y) \\ = J (x) + \alpha (A x - b, y) + \frac {\alpha^ {2}}{2} (A y, y). \\ \end{array}
$$

**命题 1.3** 设 $x\in\mathbb{R}^n$ 使得 $J(x)\le J(y)$，$\forall y\in\mathbb{R}^n$。

若命题 1.3 成立，即 $x$ 使 $J(x)$ 达到极小，则对任意 $y\in\mathbb{R}^n$，$\alpha\mapsto J(x+\alpha y)$ 在 $\alpha=0$ 取极小，故

$$
\left. \frac{d}{d\alpha} J(x+\alpha y) \right|_{\alpha=0} = 0,
$$

即 $(Ax-b,y)=0$，$\forall y\in\mathbb{R}^n$，命题 1.2 成立。

设 $A$ 对称正定，即 $(Ax,x)>0$，$\forall x\neq 0$。若命题 1.2 成立，则

$$
J (x + \alpha y) = J (x) + \frac {\alpha^ {2}}{2} (A y, y), \quad \forall y \in \mathbb{R}^n,\ \alpha \in \mathbb{R},
$$

从而 $J(x)\le J(x+\alpha y)$，命题 1.3 成立。

**小结** 若 $A$ 对称，则命题 1.3 可推出命题 1.2。若 $A$ 对称正定，则命题 1.2 可推出命题 1.3。不论 $A$ 是否对称，命题 1.1 与命题 1.2 总是等价的。若 $A$ 对称正定，则三个命题等价。

关于矩阵的特征值问题：求数 $\lambda$ 与非零向量

$$
x \in \mathbb{R}^n,\quad Ax = \lambda x.
$$

也有对应的变分问题。若 $A$ 对称正定，则 $A$ 的特征值为实数，且 $\lambda_1\ge\lambda_2\ge\cdots\ge\lambda_n>0$，$x_k$ 为对应 $\lambda_k$ 的特征向量。

定义关于 $A$ 的 Rayleigh 商

$$
R (x) = \frac {(A x , x)}{(x , x)}, \quad x \in \mathbb{R}^n,\ x \neq 0.
$$

$R(x)$ 的最大、最小值对应 $A$ 的最大、最小特征值：

$$
\lambda_ {1} = R \left(x _ {1}\right) = \max  _ {x \neq 0} R (x).
$$

其他的特征值（$k>1$）满足极小–极大（在正交补上的极大）原理：

$$
\lambda_k = R(x_k) = \max_{\substack{x \neq 0 \\ (x_1, x) = \dots = (x_{k-1}, x) = 0}} R(x).
$$

又

$$
\lambda_ {n} = R \left(x _ {n}\right) = \min  _ {x \neq 0} R (x).
$$

最大特征值可以看成在归一化条件 $(x,x)=1$ 下求 $(Ax,x)$ 的条件极值问题的解；$\lambda_k$ 可以看成在满足 $(x,x)=1$ 以及**与前 $k-1$ 个特征向量正交**（即 $(x_1,x)=\cdots=(x_{k-1},x)=0$，其中 $x_1,\dots,x_{k-1}$ 依次为对应 $\lambda_1,\dots,\lambda_{k-1}$ 的特征向量）的条件下，求 $(Ax,x)$ 的条件极值问题的解。
