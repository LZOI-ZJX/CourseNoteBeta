# 2 一维数学物理问题的变分问题

弦振动方程为

$$
\rho \frac {\partial ^ {2} u}{\partial t ^ {2}} = T \frac {\partial ^ {2} u}{\partial x ^ {2}} + F (x , t).
$$

其中 $T$ 为弦的张力，$F(x,t)$ 为外力密度。下面考虑**平衡问题**：外力 $F$ 与边界条件均与时间无关，位移 $u(x)$ 满足 $-T u''(x)=f(x)$（记静外力为 $f$）。

分析弦平衡时的位移 $u(x)$ 有两种途径。

**1. 微分方程提法（力的平衡）** $u(x)$ 满足

$$
(P _ {1}) \quad \left\{ \begin{array}{l l} - T \dfrac {d ^ {2} u}{d x ^ {2}} = f (x), & 0 < x < l, \\ u (0) = u (l) = 0. & \end{array} \right.
$$

**2. 能量与虚功（变分提法）** 由力学中的**极小位能原理**或**虚功原理**也可得到同一方程。

### 极小位能原理

弦变形前（零位移近似）可取 $W _ {1} = \displaystyle \int_ {0} ^ {l} T \, d x$；变形后弧长近似给出

$$
W _ {2} = \int_ {0} ^ {l} T \sqrt {1 + u _ {x} ^ {2}}\, d x \approx \int_ {0} ^ {l} T \left(1 + \frac {1}{2} u _ {x} ^ {2}\right) d x.
$$

应变能与外力功

$$
W _ {\text {内}} = W _ {2} - W _ {1} = \frac {1}{2} \int_ {0} ^ {l} T u _ {x} ^ {2} d x ,\quad W _ {\text {外}} = \int_ {0} ^ {l} f u d x.
$$

总位能

$$
J (u) = W _ {\text {内}} - W _ {\text {外}} = \frac {1}{2} \int_ {0} ^ {l} T u _ {x} ^ {2} d x - \int_ {0} ^ {l} f u d x.
$$

弦平衡时 $J(u)$ 最小。取容许集

$$
K = \left\{ v \in C ^ {1} [ 0, l ] \mid v (0) = v (l) = 0 \right\},
$$

得到与极小位能原理对应的变分问题

$$
(P _ {3}) \quad \begin{cases} \text{求 } u \in K \text{，使得 } \displaystyle J (u) = \frac {1}{2} \int_ {0} ^ {l} T u _ {x} ^ {2} d x - \int_ {0} ^ {l} f u d x \text{ 最小}, \\ K = \left\{ u \mid u \in C ^ {1} [ 0, l ], u (0) = u (l) = 0 \right\}. \end{cases}
$$

设 $u$ 为极小解，对任意 $\eta \in K$ 令 $\varphi (\alpha) = J (u + \alpha \eta)$。则

$$
\varphi^ {\prime} (0) = \int_ {0} ^ {l} T u _ {x} \eta_ {x} d x - \int_ {0} ^ {l} f \eta d x = 0.
$$

分部积分并利用 $\eta (0) = \eta (l) = 0$，得

$$
\int_ {0} ^ {l} \left[ - T \frac {d ^ {2} u}{d x ^ {2}} - f (x) \right] \eta d x = 0, \quad \forall \eta \in K.
$$

由此可与 $(P _ {1})$ 的微分方程形式联系起来。

### 虚功原理

虚功原理给出（等价地经分部积分）：

$$
\int_ {0} ^ {l} T u _ {x} v _ {x} d x = \int_ {0} ^ {l} f v d x, \quad \forall v \in K.
$$

记

$$
D (u, v) = \int_ {0} ^ {l} T u _ {x} v _ {x} d x, \quad F (v) = \int_ {0} ^ {l} f v d x,
$$

则**虚功方程**即 Galerkin 变分形式

$$
(P _ {2}) \quad \begin{cases} \text{求 } u \in K \text{ 使得} \\ D (u, v) - F (v) = 0, \quad \forall v \in K. \end{cases}
$$

在一般区间 $[a,b]$ 与更一般的算子 $L$ 下，仍用 $(P _ {1})$、$(P _ {2})$、$(P _ {3})$ 记微分形式与两类变分形式，容许空间常记为 $S_0^1$（见下）。

---

## 2.1 两点边值问题的变分形式

本节在 $[a,b]$ 上考虑含系数 $p(x),q(x)$ 的模型；当 $p\equiv T$, $q\equiv 0$, $[a,b]=[0,l]$ 且两端位移为零时，即回到弦的平衡方程。

考虑一般的二阶自伴微分算子 $L$，

$$
L u = - \frac {d}{d x} \left(p (x) \frac {d u}{d x}\right) + q (x) u \tag {2.1}
$$

其中 $p \in C ^ {1} [ a, b ]$, $p (x) \geq p _ {0} > 0$, $q \in C [ a, b ]$, $q (x) \geq 0$.

设 $f \in C [ a, b ]$。对方程 $L u = f$ 加上两端边界条件，得更一般的两点边值问题

$$
(P _ {1}) \quad \left\{ \begin{array}{l l} - \dfrac {d}{d x} \left(p (x) \dfrac {d u}{d x}\right) + q (x) u = f (x), & a < x < b, \\ u (a) = 0, & \dfrac {d u (b)}{d x} = 0. \end{array} \right.
$$

$(P _ {1})$ 中的方程可写成 $L u - f = 0$. 方程两边乘以函数 $v$ 再在 $(a,b)$ 上积分。记内积

$$
(u, v) = \int_ {a} ^ {b} u (x) v (x) d x,
$$

则 $(L u - f, v) = 0$，即

$$
- \int_ {a} ^ {b} v \frac {d}{d x} \left(p \frac {d u}{d x}\right) d x + \int_ {a} ^ {b} q u v d x - \int_ {a} ^ {b} f v d x = 0.
$$

第一项分部积分，并利用右端边界条件 $\displaystyle \frac {d u (b)}{d x} = 0$，得

$$
\int_ {a} ^ {b} p \frac {d u}{d x} \frac {d v}{d x} d x - \left. \left(p \frac {d u}{d x} v\right) \right| _ {x = a} + \int_ {a} ^ {b} q u v d x - \int_ {a} ^ {b} f v d x = 0.
$$

若规定试探函数满足 $v (a) = 0$，则边界项在 $x = a$ 处消失。令

$$
S _ {0} ^ {1} = \left\{ v \mid v \in C ^ {1} [ a, b ], v (a) = 0 \right\},
$$

并定义

$$
D (u, v) = \int_ {a} ^ {b} \left(p \frac {d u}{d x} \frac {d v}{d x} + q u v\right) d x, \quad F (v) = \int_ {a} ^ {b} f v d x.
$$

则得到与 $(P _ {1})$ 相应的 Galerkin 变分问题

$$
(P _ {2}) \quad \begin{cases} \text{求 } u \in S _ {0} ^ {1} \text{ 使得} \\ D (u, v) - F (v) = 0, \quad \forall v \in S _ {0} ^ {1}. \end{cases}
$$

从 $(P _ {2})$ 可见：给定 $f$ 后，每取定 $v _ {1} \in S _ {0} ^ {1}$，$F (v _ {1})$ 为实数，故 $F$ 是 **$S _ {0} ^ {1}$ 到 $\mathbb{R}$ 的映射（线性泛函）**，满足

$$
F (v _ {1} + v _ {2}) = F (v _ {1}) + F (v _ {2}), \quad F (c v) = c F (v), \quad \forall c \in \mathbb{R}.
$$

同理，$D$ 是 **$S _ {0} ^ {1} \times S _ {0} ^ {1}$ 到 $\mathbb{R}$ 的映射（双线性泛函）**；当 $u$ 与 $v$ 中有一个固定时，$D (u, v)$ 关于另一个是线性的。若 $D (u, v) = D (v, u)$，则称 $D$ 为 $S _ {0} ^ {1} \times S _ {0} ^ {1}$ 上的**对称**双线性泛函。在许多问题中还满足 $D (u, u) \geq 0$，$\forall u \in S _ {0} ^ {1}$.

若 $u$ 是 $(P _ {1})$ 的（充分光滑）解，则 $u \in S _ {0} ^ {1}$（因 $u (a) = 0$），且通过分部积分可知 $u$ 满足 $(P _ {2})$. 反之，由 $(P _ {2})$ 分部积分可写为（对一切 $v \in S _ {0} ^ {1}$）

$$
\left. - \int_ {a} ^ {b} v \frac {d}{d x} \left(p \frac {d u}{d x}\right) d x + \left(p \frac {d u}{d x} v\right) \right| _ {x = b} - \left. \left(p \frac {d u}{d x} v\right) \right| _ {x = a} + \int_ {a} ^ {b} q u v d x - \int_ {a} ^ {b} f v d x = 0. \tag {2.5}
$$

当 $v \in C _ {0} ^ {1} [ a, b ] \subset S _ {0} ^ {1}$ 时，$v (a) = v (b) = 0$，上式给出

$$
\int_ {a} ^ {b} (L u - f) v d x = 0, \quad \forall v \in C _ {0} ^ {1} [ a, b ].
$$

**引理 2.1（变分法基本引理）** 设 $u \in C [ a, b ]$，且 $\displaystyle \int_ {a} ^ {b} u (x) v (x) d x = 0$ 对一切 $v \in C _ {0} ^ {1} [ a, b ]$ 成立，则 $u (x) \equiv 0$（$a \leq x \leq b$）。

由此可得 $L u - f = 0$，再将方程代回 (2.5)，并利用 $v \in S _ {0} ^ {1}$ 任意，得到右端自然边界条件

$$
\left. p (x) \frac {d u}{d x} \right| _ {x = b} = 0 \quad (\text{在 } x = b \text{ 处}).
$$

（若 $p (b) > 0$，则等价于 $\displaystyle \frac {d u (b)}{d x} = 0$。）左端 $u (a) = 0$ 为强加边界条件；$x = b$ 处为**自然边界条件**（力的边界条件），与第一类边界条件相对。

### Ritz 变分问题（极小势能）

考虑

$$
\begin{array}{l} J (u) = \frac {1}{2} D (u, u) - F (u) \tag {2.6} \\ = \dfrac {1}{2} \displaystyle \int_ {a} ^ {b} \left[ p \left(\frac {d u}{d x}\right) ^ {2} + q u ^ {2} - 2 f u \right] d x, \\ \end{array}
$$

它表示在所给左端约束下问题的总势能。由极小势能原理，考虑

$$
(P _ {3}) \quad \begin{cases} \text{求 } u \in S _ {0} ^ {1} \text{ 使得} \\ J (u) \leq J (w), \quad \forall w \in S _ {0} ^ {1}. \end{cases}
$$

设 $u$ 是 $(P _ {3})$ 的解，则对任意 $v \in S _ {0} ^ {1}$，$w = u + \alpha v \in S _ {0} ^ {1}$，有 $J (u) \leq J (u + \alpha v)$. 又

$$
\begin{array}{l} J (u + \alpha v) = \frac {1}{2} D (u + \alpha v, u + \alpha v) - F (u + \alpha v) \\ = J (u) + \alpha \left[ D (u, v) - F (v) \right] + \frac {\alpha^ {2}}{2} D (v, v). \tag {2.7} \\ \end{array}
$$

视 $J (u + \alpha v)$ 为 $\alpha$ 的一元函数，在 $\alpha = 0$ 取极小，故

$$
\left. \frac {d}{d \alpha} J (u + \alpha v) \right| _ {\alpha = 0} = 0 \implies D (u, v) - F (v) = 0, \quad \forall v \in S _ {0} ^ {1}.
$$

即 $(P _ {3})$ 的解必是 $(P _ {2})$ 的解。反过来，若 $u$ 是 $(P _ {2})$ 的解，则由 (2.7) 及 $D (v, v) \geq 0$，有 $J (u + \alpha v) \geq J (u)$，从而 $J (u) \leq J (w)$，$\forall w \in S _ {0} ^ {1}$，即 $(P _ {2})$ 的解也是 $(P _ {3})$ 的解。

**小结** $(P _ {1})$、$(P _ {2})$、$(P _ {3})$ 在上述光滑性与正定性假设下两两等价：$(P _ {1})$ 是微分方程两点边值问题；$(P _ {2})$ 为 Galerkin 变分问题（虚功原理）；$(P _ {3})$ 为 Ritz 变分问题（最小势能原理）。

---

## 2.2 非齐次约束边界条件的处理

考虑第一类非齐次边界条件：

$$
\left\{ \begin{array}{l l} - \dfrac {d}{d x} \left(p (x) \dfrac {d u}{d x}\right) + q (x) u = f (x), & a < x < b, \\ u (a) = u _ {0}, & \dfrac {d u (b)}{d x} = 0. \end{array} \right.
$$

为导出相应的变分问题，可先**齐次化**边界：令 $\bar {u} (x) = u (x) - u _ {0}$，则

$$
\left\{ \begin{array}{l l} - \dfrac {d}{d x} \left(p (x) \dfrac {d \bar {u}}{d x}\right) + q (x) \bar {u} = f (x) - q (x) u _ {0}, & a < x < b, \\ \bar {u} (a) = 0, & \dfrac {d \bar {u} (b)}{d x} = 0. \end{array} \right.
$$

对 $\bar {u}$ 在 $S _ {0} ^ {1}$ 上写出 $(P _ {2})$、$(P _ {3})$（右端项与系数按上式替换），再把 $u = \bar {u} + u _ {0}$ 代回即可得到关于原 $u$ 的变分提法。也可直接取**仿射容许集**

$$
S = \left\{ v \in C ^ {1} [ a, b ] \mid v (a) = u _ {0} \right\},
$$

求 $u \in S$ 使得 Galerkin 方程 $D (u, v) - F (v) = 0$ 对一切满足 $v (a) = 0$ 的试探函数 $v$ 成立（或写出与 (2.6) 同型的 Ritz 泛函在 $S$ 上极小）。其中 $D$、$F$ 的定义与 2.1 节相同。

---

## 2.3 第二、三类边界条件

$$
\left\{ \begin{array}{l} - \dfrac {d}{d x} \left(p (x) \dfrac {d u}{d x}\right) + q (x) u = f (x), \quad a < x < b, \\ u (a) = 0, \quad p (b) \dfrac {d u (b)}{d x} + \alpha u (b) = g. \end{array} \right.
$$

其中 $p (x) \geq p _ {0} > 0$, $q (x) \geq 0$，常数 $\alpha \geq 0$. 在 $x = b$ 的边界条件同时涵盖第二、第三类边界；可仿照 2.1 节推导出相应的变分问题，例如取

$$
D (u, v) = \int_ {a} ^ {b} \left[ p \frac {d u}{d x} \frac {d v}{d x} + q u v \right] d x + \alpha u (b) v (b),
$$

$$
F (v) = \int_ {a} ^ {b} f v d x + g v (b).
$$
