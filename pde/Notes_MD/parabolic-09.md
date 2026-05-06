## 第二节 初边值问题

### 2.1、第一类边界条件

考虑扩散方程的第一边值问题

$$
\left\{ \begin{array}{l l} \frac {\partial u}{\partial t} - a \frac {\partial^ {2} u}{\partial x ^ {2}} = 0 & \quad 0 <   x <   1 \quad t > 0 \\ u (x, 0) = f (x) & \quad 0 <   x <   1 \\ u (0, t) = \varphi (t) & \quad t \geq 0 \\ u (1, t) = \psi (t) & \quad t \geq 0 \end{array} \right.
$$

将变量 $x$ 的计算区域 $[0,1]$ 进行剖分:

$$
0 = x _ {0} <   x _ {1} <   \dots <   x _ {J} = 1
$$

$$
\begin{array}{l} h = \frac {1}{J}, x _ {j} = j h, x _ {J} = J h = 1 \\ \boldsymbol {t} _ {n} = \boldsymbol {\tau n}, \quad \boldsymbol {n} \geq 0 \\ \end{array}
$$

初值离散: $u_{j}^{0} = f(x_{j}) = f_{j} \quad j = 1,2,\dots J - 1$ 

边值离散: $\left\{ \begin{array}{l l} u _ {0} ^ {n} = \varphi \left(t _ {n}\right) & n \geq 0 \\ u _ {J} ^ {n} = \psi \left(t _ {n}\right) & n \geq 0 \end{array} \right.$ 

内点用相应的差分格式即可求解上述方程。

### 2.2、第三类边界条件

$$
\left\{ \begin{array}{l l} \frac {\partial u}{\partial t} - a \frac {\partial^ {2} u}{\partial x ^ {2}} = 0 & 0 <   x <   1 \quad t > 0 \\ u (x, 0) = f (x) & 0 <   x <   1 \\ \frac {\partial u}{\partial x} (0, t) = \alpha u (0, t) + \mu (t) & t > 0 \\ \frac {\partial u}{\partial x} (1, t) = \beta u (1, t) + \gamma (t) & t > 0 \end{array} \right.
$$

网格部分、初值离散和前面一样，边界条件的离散可按如下两种方法处理：

(1)在点 $\left(0, t_{n}\right)$ 和点 $\left(1, t_{n}\right)$ 分别用向前和向后差商

近似 $\frac{\partial u}{\partial x}$ , 有:

$$
\left\{ \begin{array}{l} \frac {u _ {1} ^ {n} - u _ {0} ^ {n}}{h} = \alpha u _ {0} ^ {n} + \mu (t _ {n}) \\ \frac {u _ {J} ^ {n} - u _ {J - 1} ^ {n}}{h} = \beta u _ {J} ^ {n} + \gamma (t _ {n}) \end{array} \right.
$$

与内点差分格式一起，即可求解，但边界条件的离散是一阶的，不适用与二阶格式。

(2)用中心差商近似 $\frac{\partial u}{\partial x}$ , 有：

$$
\left\{ \begin{array}{l} \frac {u _ {1} ^ {n} - u _ {- 1} ^ {n}}{2 h} = \alpha u _ {0} ^ {n} + \mu (t _ {n}) \\ \frac {u _ {J + 1} ^ {n} - u _ {J - 1} ^ {n}}{2 h} = \beta u _ {J} ^ {n} + \gamma (t _ {n}) \end{array} \right.
$$

这种离散是二阶的, $E = O(h^{2})$ , 但式中含有区域之外的点 $\left(x_{-1}, t_{n}\right), \left(x_{J+1}, t_{n}\right)$ , 需要借助内点差分格式消去 $u_{-1}^{n}, u_{J+1}^{n}$ 

例如：内点计算用古典显式格式:

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} - a \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}} = 0
$$

分别取 $j = 0, j = J$ , 得:

$$
\begin{array}{l} u _ {0} ^ {n + 1} = u _ {0} ^ {n} + a \lambda (u _ {1} ^ {n} - 2 u _ {0} ^ {n} + u _ {- 1} ^ {n}) \\ u _ {J} ^ {n + 1} = u _ {J} ^ {n} + a \lambda \left(u _ {J + 1} ^ {n} - 2 u _ {J} ^ {n} + u _ {J - 1} ^ {n}\right) \\ \end{array}
$$

与前面的边界离散格式联立，解得:

$$
\left\{ \begin{array}{l} u _ {0} ^ {n + 1} = \big [ 1 - 2 a \lambda (1 + \alpha h) \big ] u _ {0} ^ {n} + 2 a \lambda u _ {1} ^ {n} - 2 a \lambda h \, \mu (t _ {n}) \\ u _ {J} ^ {n + 1} = \big [ 1 - 2 a \lambda (1 - \beta h) \big ] u _ {J} ^ {n} + 2 a \lambda u _ {J - 1} ^ {n} + 2 a \lambda h \, \gamma (t _ {n}) \end{array} \right.
$$

这样的边界值是二阶精度的, 与内点格式联立可求解方程。

### 2.3、Saul'ev 算法

$$
\frac {u (x _ {j} , t _ {n + 1}) - u (x _ {j} , t _ {n})}{\tau} = \frac {a}{h} \left(\frac {\partial u}{\partial x} \bigg | _ {(x _ {j + \frac {1}{2}},\, t _ {n})} - \frac {\partial u}{\partial x} \bigg | _ {(x _ {j - \frac {1}{2}},\, t _ {n})}\right) + O (\tau + h ^ {2})
$$

$$
\frac {\partial u}{\partial x} \bigg | _ {(x _ {j + \frac {1}{2}},\, t _ {n})} = \frac {\partial u}{\partial x} \bigg | _ {(x _ {j + \frac {1}{2}},\, t _ {n + 1})} - \tau \frac {\partial^ {2} u}{\partial x \partial t} \bigg | _ {(x _ {j + \frac {1}{2}},\, t _ {n} + \theta \tau)}
$$

$$
\frac {\partial u}{\partial x} \bigg | _ {(x _ {j - \frac {1}{2}},\, t _ {n})} = \frac {\partial u}{\partial x} \bigg | _ {(x _ {j - \frac {1}{2}},\, t _ {n + 1})} - \tau \frac {\partial^ {2} u}{\partial x \partial t} \bigg | _ {(x _ {j - \frac {1}{2}},\, t _ {n} + \theta \tau)}
$$

$$
\frac {u (x _ {j} , t _ {n + 1}) - u (x _ {j} , t _ {n})}{\tau} = \frac {a}{h} \left(\frac {\partial u}{\partial x} \bigg | _ {(x _ {j + \frac {1}{2}},\, t _ {n})} - \frac {\partial u}{\partial x} \bigg | _ {(x _ {j - \frac {1}{2}},\, t _ {n + 1})}\right) + O (\tau + h ^ {2})
$$

$$
\frac {u (x _ {j} , t _ {n + 1}) - u (x _ {j} , t _ {n})}{\tau} = \frac {a}{h} \left(\frac {\partial u}{\partial x} \bigg | _ {(x _ {j + \frac {1}{2}},\, t _ {n + 1})} - \frac {\partial u}{\partial x} \bigg | _ {(x _ {j - \frac {1}{2}},\, t _ {n})}\right) + O (\tau + h ^ {2})
$$

进一步采用中心差分:

$$
\frac {\partial u}{\partial x} \bigg | _ {(x _ {j + \frac {1}{2}},\, t _ {n})} = \frac {u (x _ {j + 1} , t _ {n}) - u (x _ {j} , t _ {n})}{h} + O (h ^ {2})
$$

#### 显式计算式（交替方向）

$$
\begin{array}{l} u _ {j} ^ {n + 1} = \frac {1 - a \lambda}{1 + a \lambda} u _ {j} ^ {n} + \frac {a \lambda}{1 + a \lambda} \left(u _ {j - 1} ^ {n + 1} + u _ {j + 1} ^ {n}\right) \quad (1) \\ u _ {j} ^ {n + 1} = \frac {1 - a \lambda}{1 + a \lambda} u _ {j} ^ {n} + \frac {a \lambda}{1 + a \lambda} \left(u _ {j + 1} ^ {n + 1} + u _ {j - 1} ^ {n}\right) \quad (2) \\ n = 1, 3, 5 \dots \text {用 (2) 式 ,} n = 2, 4, 6 \dots \text {用 (1) 式 。} \\ \end{array}
$$

性质: 1、截断误差 $O\left(\frac{\tau}{h} + \tau + h^{2}\right)$ 

2、无条件稳定的。

### 2.4、分组显式方法

令 $a \lambda = r$ , 在网格点 $(x_{j+1}, t_n)$ 采用 (1), 在 $(x_j, t_n)$ 采用 (2)

$$
\left\{ \begin{array}{l} (1 + r) u _ {j} ^ {n + 1} - r u _ {j + 1} ^ {n + 1} = (1 - r) u _ {j} ^ {n} + r u _ {j - 1} ^ {n} \\ (1 + r) u _ {j + 1} ^ {n + 1} - r u _ {j} ^ {n + 1} = (1 - r) u _ {j + 1} ^ {n} + r u _ {j + 2} ^ {n} \end{array} \right.
$$

解得:

$$
\left\{ \begin{array}{l} u _ {j} ^ {n + 1} = \frac {1}{1 + 2 r} [ r (1 + r) u _ {j - 1} ^ {n} + (1 - r ^ {2}) u _ {j} ^ {n} + r (1 - r) u _ {j + 1} ^ {n} + r ^ {2} u _ {j + 2} ^ {n} ] \\ u _ {j + 1} ^ {n + 1} = \frac {1}{1 + 2 r} [ r ^ {2} u _ {j - 1} ^ {n} + r (1 - r) u _ {j} ^ {n} + (1 - r ^ {2}) u _ {j + 1} ^ {n} + r (1 + r) u _ {j + 2} ^ {n} ] \end{array} \right.
$$

此方法易于并行计算

## 第三节 对流扩散方程

对流扩散方程的初值问题:

$$
\left\{ \begin{array}{l l} \frac {\partial u}{\partial t} + a \frac {\partial u}{\partial x} = \nu \frac {\partial^ {2} u}{\partial x ^ {2}} & - \infty <   x <   + \infty , t > 0 \\ u (x, 0) = f (x) & - \infty <   x <   + \infty , \end{array} \right.
$$

其中 $a, \nu$ 为常数，且 $\nu > 0$ .

对流占优问题: $a \gg \nu$ 

### 3.1、中心显式差分格式

$$
\frac {\boldsymbol {u} _ {j} ^ {n + 1} - \boldsymbol {u} _ {j} ^ {n}}{\tau} + a \frac {\boldsymbol {u} _ {j + 1} ^ {n} - \boldsymbol {u} _ {j - 1} ^ {n}}{2 h} = \nu \frac {\boldsymbol {u} _ {j + 1} ^ {n} - 2 \boldsymbol {u} _ {j} ^ {n} + \boldsymbol {u} _ {j - 1} ^ {n}}{h ^ {2}}
$$

其截断误差为: $E = O\left(\tau + h^{2}\right)$ , 而且当 $\nu = 0$ 时,格式为近似对流方程的无条件不稳定格式; 当 $a = 0$ 时，格式退化为近似扩散方程的古典显式格式，只有 $\nu \dfrac{\tau}{h^{2}} \leq \dfrac{1}{2}$ 时，格式才稳定。

下面讨论稳定性：

设 $\lambda = a \tau / h, \mu = \nu \tau / h^{2}$ , 格式改写成:

$$
u _ {j} ^ {n + 1} = u _ {j} ^ {n} - \frac {1}{2} \lambda (u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}) + \mu (u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n})
$$

增长因子为：

$$
G (\lambda , \mu , k) = 1 - 2 \mu (1 - \cos k h) - i \lambda \sin k h
$$

$$
\begin{array}{l} \left| G (\lambda , \mu) \right| ^ {2} = \lambda^ {2} \sin^ {2} k h + 1 + 4 \mu^ {2} (1 - \cos k h) ^ {2} - 4 \mu (1 - \cos k h) \\ = 1 - (1 - \cos k h) \left[ 4 \mu - 4 \mu^ {2} (1 - \cos k h) - \lambda^ {2} (1 + \cos k h) \right] \\ \end{array}
$$

只需验证 $|G| \leq 1$ , 由于

$1 - \cos {kh} \geq 0$ , 条件转化为:

$$
4 \mu - 4 \mu^ {2} (1 - \cos k h) - \lambda^ {2} (1 + \cos k h) \geq 0
$$

即 $4 \mu - 2 \lambda^{2} - (4 \mu^{2} - \lambda^{2})(1 - \cos k h) \geq 0$ 

由于 $\frac{1 - \cos kh}{2} \in [0, 1]$ , 上述不等式转化为

$$
4 \mu - 2 \lambda^ {2} \geq 0, \quad 4 \mu - 2 \lambda^ {2} + 2 (\lambda^ {2} - 4 \mu^ {2}) \geq 0
$$

将 $\lambda, \mu$ 代入，即得条件：

$$
\tau \leq \frac {2 \nu}{a ^ {2}}, \quad \nu \frac {\tau}{h ^ {2}} \leq \frac {1}{2}
$$

此两不等式为中心显式格式稳定的条件。

注：(1) $\nu \dfrac{\tau}{h^{2}} \leq \dfrac{1}{2}$，是相应扩散方程的稳定性条件。

(2) 对于对流占优的问题，稳定性约束过严。

### 3.2、修正中心差分格式

假设对流扩散方程的解充分光滑，对方程两边同时对 $t$ 求导有：

$$
\frac {\partial^ {2} u}{\partial t ^ {2}} = \nu ^ {2} \frac {\partial^ {4} u}{\partial x ^ {4}} - 2 a \nu \frac {\partial^ {3} u}{\partial x ^ {3}} + a ^ {2} \frac {\partial^ {2} u}{\partial x ^ {2}}
$$

中心显格式的截断误差为：

$$
\begin{array}{l} E = \frac {\partial u}{\partial t} + a \frac {\partial u}{\partial x} - \nu \frac {\partial^ {2} u}{\partial x ^ {2}} \\ + \frac {\tau}{2} \frac {\partial^ {2} u}{\partial t ^ {2}} + h ^ {2} \frac {\partial^ {4} u}{\partial x ^ {4}} + O (\tau^ {2} + h ^ {4}) \\ = \frac {\partial u}{\partial t} + a \frac {\partial u}{\partial x} - \nu \frac {\partial^ {2} u}{\partial x ^ {2}} \\ + \frac {\tau}{2} \left(\nu ^ {2} \frac {\partial^ {4} u}{\partial x ^ {4}} - 2 a \nu \frac {\partial^ {3} u}{\partial x ^ {3}} + a ^ {2} \frac {\partial^ {2} u}{\partial x ^ {2}}\right) + h ^ {2} \frac {\partial^ {4} u}{\partial x ^ {4}} + O \left(\tau^ {2} + h ^ {4}\right) \\ \end{array}
$$

$$
\begin{array}{l} = \frac {\partial u}{\partial t} + a \frac {\partial u}{\partial x} - (\nu - \frac {a ^ {2} \tau}{2}) \frac {\partial^ {2} u}{\partial x ^ {2}} \\ + \frac {\tau}{2} \left(\nu ^ {2} \frac {\partial^ {4} u}{\partial x ^ {4}} - 2 a \nu \frac {\partial^ {3} u}{\partial x ^ {3}}\right) + h ^ {2} \frac {\partial^ {4} u}{\partial x ^ {4}} + O (\tau^ {2} + h ^ {4}) \\ \end{array}
$$

由此有当 $\tau$ 不趋近于 0 时, 差分格式与下面的方程相容

$$
\frac {\partial u}{\partial t} + a \frac {\partial u}{\partial x} = (\nu - \frac {\tau}{2} a ^ {2}) \frac {\partial^ {2} u}{\partial x ^ {2}}
$$

然而在数值计算的时候 $\tau$ 是一不为 0 的固定常数,于是导致了扩散效应的损失, 特别在 $\nu - \frac{\tau}{2} a^{2} \rightarrow 0$ 时, 中心差分格式**相容于对流方程**, 而此时中心格式是绝对不稳定的, 为了减少扩散效应的损失,在相应的扩散项增加扩散的系数为 $\nu + \frac{\tau}{2} a^{2}$ 。这样得到如下差分格式:

#### 修正中心显式格式

$$
\boxed { \begin{array}{c} \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}}{2 h} = (\nu + \frac {\tau}{2} a ^ {2}) \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}} \end{array} }
$$

稳定性分析完全类似于中心差分格式,

令 $\bar{\nu} = \nu + \frac{\tau}{2} a^{2}$ , 得

$$
\tau \leq \frac {2 \bar {\nu}}{a ^ {2}}, \quad \bar {\nu} \frac {\tau}{h ^ {2}} \leq \frac {1}{2}
$$

所以修正中心显式格式的稳定性条件是

$$
\nu \frac {\tau}{h ^ {2}} + \frac {1}{2} \left(a \frac {\tau}{h}\right) ^ {2} \leq \frac {1}{2}
$$

### 3.3、迎风差分格式

在中心显式差分格式的稳定性条件中，当 $\dfrac{\nu}{a^{2}}$ 很小时，$\tau$ 只能取得很小, 显式格式计算量很大。

当 $\nu = 0$ (极限情况) 时, 微分方程化为对流方程,中心显式格式转化为绝对不稳定的。

考虑迎风格式:

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j} ^ {n} - u _ {j - 1} ^ {n}}{h} = \nu \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}} \quad (a > 0)
$$

即: $u_{j}^{n+1} = (\lambda + \mu) u_{j-1}^{n} + (1 - \lambda - 2\mu) u_{j}^{n} + \mu u_{j+1}^{n}, (a > 0)$ 

其增长因子

$$
G = 1 - (\lambda + 2 \mu) (1 - \cos k h) - i \lambda \sin k h
$$

$$
\begin{array}{l} \left| G \right| ^ {2} = \lambda^ {2} \sin^ {2} k h + (1 - (\lambda + 2 \mu) (1 - \cos k h)) ^ {2} \\ = 1 - (1 - \cos k h) \left\{2 (\lambda + 2 \mu) - 2 \lambda^ {2} + \left(\lambda^ {2} - (\lambda + 2 \mu) ^ {2}\right) (1 - \cos k h) \right\} \\ \end{array}
$$

类似的迎风格式稳定的充要条件是 $|G| \leq 1$ , 相当于要求:

$$
2 (\lambda + 2 \mu) - 2 \lambda^ {2} + (\lambda^ {2} - (\lambda + 2 \mu) ^ {2}) (1 - \cos k h) \geq 0
$$

由于 $1 - \cos kh \in [0,2]$ , 条件化为:

$$
\begin{array}{l} 2 (\lambda + 2 \mu) - 2 \lambda^ {2} \geq 0 \text {和} \\ 2 (\lambda + 2 \mu) - 2 \lambda^ {2} + 2 \left(\lambda^ {2} - (\lambda + 2 \mu) ^ {2}\right) \geq 0 \\ \end{array}
$$

而 $2(\lambda + 2\mu) - 2\lambda^2 + 2(\lambda^2 - (\lambda + 2\mu)^2)$ 

$$
= 2 (\lambda + 2 \mu) (1 - (\lambda + 2 \mu)) \geq 0 \Leftrightarrow 1 - (\lambda + 2 \mu) \geq 0
$$

即 $\lambda + 2 \mu \leq 1$ , 此时有 $\lambda < 1, (\lambda + 2 \mu - \lambda^{2}) \geq 0$ 自然成立。

所以第二个不等式可以保证第一个不等式。

所以迎风格式的稳定性条件是

$$
\tau \leq \frac {h ^ {2}}{2 \nu + a h}
$$

当 $a < 0$ 时, 情况类似, 稳定性条件是

$$
\tau \leq \frac {h ^ {2}}{2 \nu + | a | h}
$$

也可以利用中心显格式来讨论稳定性，于是将上面格式改为：

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}}{2 h} = \left(\nu + \frac {a h}{2}\right) \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}
$$

取 $\bar{\nu} = \nu + \frac{a h}{2}$（$a>0$ 时与迎风格式通过截断误差引入的**人工粘度**项一致）, 则变为中心格式, 于是

$$
\tau \leq \frac {2 \bar {\nu}}{a ^ {2}} \quad (1), \qquad \bar {\nu} \frac {\tau}{h ^ {2}} \leq \frac {1}{2} \quad (2)
$$

通过简单的推导, 可以发现第一个稳定条件可以由第2个条件推出, 于是迎风格式的稳定条件就是 (2).

### 3.4、Samarskii 格式

Samarskii格式是具有迎风效应的关于空间的二阶格式,为了简单方便,设 $a > 0$ , 先对方程作扰动, 得到另外一对流扩散方程

$$
\frac {\partial u}{\partial t} + a \frac {\partial u}{\partial x} = \frac {1}{1 + R} \nu \frac {\partial^ {2} u}{\partial x ^ {2}}
$$

其中 $R = \dfrac{ah}{2\nu}$。

对上面的方程构造迎风格式

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j} ^ {n} - u _ {j - 1} ^ {n}}{h} = \frac {1}{1 + R} \nu \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}
$$

称为逼近对流扩散方程的Samarskii格式.

为便于截断误差分析，将该差分格式两端重组并记为 $\alpha_j^n + \beta_j^n = 0$，其中：

$$
\alpha_ {j} ^ {n} = \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} - \nu \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}
$$

$$
\beta_ {j} ^ {n} = a \frac {u _ {j} ^ {n} - u _ {j - 1} ^ {n}}{h} - (\frac {1}{1 + R} - 1) \quad \nu \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}
$$

#### 由 Taylor 公式

$$
\alpha_ {j} ^ {n} = (\frac {\partial u}{\partial t}) _ {j} ^ {n} - \nu (\frac {\partial^ {2} u}{\partial x ^ {2}}) _ {j} ^ {n} + O (\tau + h ^ {2})
$$

$$
\beta_ {j} ^ {n} = a \left(\frac {\partial u}{\partial x}\right) _ {j} ^ {n} - \frac {a h}{2} \left(\frac {\partial^ {2} u}{\partial x ^ {2}}\right) _ {j} ^ {n} + \frac {R}{1 + R} \nu \left(\frac {\partial^ {2} u}{\partial x ^ {2}}\right) _ {j} ^ {n} + O \left(h ^ {2}\right)
$$

令 $R = C h$，将 $\beta_j^n$ 中 $\partial^2 u/\partial x^2$ 的系数合并为 $\dfrac{Ch}{1+Ch}\nu - \dfrac{ah}{2}$，通分后得

$$
\begin{array}{l} \beta_ {j} ^ {n} = a \left(\frac {\partial u}{\partial x}\right) _ {j} ^ {n} + \frac {h (2 \nu C - a) - C a h ^ {2}}{2 (1 + C h)} \left(\frac {\partial^ {2} u}{\partial x ^ {2}}\right) _ {j} ^ {n} + O \left(h ^ {2}\right) \\ = a \left(\frac {\partial u}{\partial x}\right) _ {j} ^ {n} - \frac {R ^ {2}}{1 + R} \nu \left(\frac {\partial^ {2} u}{\partial x ^ {2}}\right) _ {j} ^ {n} + O \left(h ^ {2}\right) \quad \left(\therefore C = \frac {a}{2 \nu}\right) \\ \end{array}
$$

于是截断误差有 $O\left(\tau + h^{2}\right)$ 

类似迎风格式的稳定性分析，可以得到此格式的稳定性条件：

$$
\left[ \frac {a h}{2} + \frac {\nu}{1 + \frac {a h}{2 \nu}} \right] \frac {\tau}{h ^ {2}} \leq \frac {1}{2}
$$

### 3.5、指数型差分格式

考虑定常形式的对流扩散方程 $d + a \frac{\partial u}{\partial x} = \nu \frac{\partial^2 u}{\partial x^2}$（其中在局部采用冻结系数，将时间微商视作常数源项 $d=\dfrac{\partial u}{\partial t}$），其解析解为 $u = \beta_{1} e^{\lambda x} + \beta_{2} - \dfrac{d}{a} x$，其中 $\lambda = \dfrac{a}{\nu}$。

对空间进行剖分，有

$$
u _ {j - 1} = \beta_ {1} e ^ {\lambda (j - 1) h} + \beta_ {2} - \frac {d}{a} (j - 1) h
$$

$$
u _ {j + 1} = \beta_ {1} e ^ {\lambda (j + 1) h} + \beta_ {2} - \frac {d}{a} (j + 1) h
$$

可以求解 $\beta_{1}, \beta_{2}$ , 并且代入解析解得:

$$
u _ {j} = \frac {u _ {j + 1} + e ^ {\lambda h} u _ {j - 1}}{1 + e ^ {\lambda h}} - d \frac {h \cdot (e ^ {\lambda h} - 1)}{a \cdot (e ^ {\lambda h} + 1)}
$$

上式通过恒等变形可以改写为

$$
d = - a \frac {u _ {j + 1} - u _ {j - 1}}{2 h} - \frac {a h}{2} \frac {1 + e ^ {\lambda h}}{1 - e ^ {\lambda h}} \frac {u _ {j + 1} - 2 u _ {j} + u _ {j - 1}}{h ^ {2}}
$$

用时间差商代替 $d$ , 得到对流扩散方程的差分格式

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}}{2 h} = - \frac {a h}{2} \frac {1 + e ^ {\lambda h}}{1 - e ^ {\lambda h}} \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}
$$

由前一式，记**辅助因子**（仅用于恒等变形）

$$
\omega := \frac{1 + e^{\lambda h}}{1 - e^{\lambda h}} = -\coth \frac{\lambda h}{2} = -\coth \frac{a h}{2\nu} \quad \left(\lambda = \frac{a}{\nu}\right),
$$

则扩散项可写为 $-\dfrac{ah}{2}\,\omega\,\dfrac{\delta_x^2 u}{h^2}$；若要把方程写成 **$\nu \sigma$** 与 Laplacian 相乘的形式，应令

$$
\nu \sigma = \frac{ah}{2}\coth\frac{ah}{2\nu}, \qquad \sigma = \frac{ah}{2\nu}\coth\frac{ah}{2\nu}.
$$

二者满足 $\nu\sigma = -\dfrac{ah}{2}\omega$（注意：**此处的 $\sigma$ 并不**等于上式中的 $\omega$，切勿混用）。

#### 指数型差分格式（放大因子形式）

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}}{2 h} = \nu \sigma \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}
$$

截断误差为 $O\left(h^{2} + \tau\right)$ 。

与中心显格式比较可得出稳定性条件

稳定性条件为: $\nu \sigma \frac{\tau}{h^{2}} \leq \frac{1}{2}$ 

对流扩散方程的指数型差分格式的两种形式:

$$
\begin{array}{l} (1) \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}}{2 h} = - \frac {a h}{2} \frac {1 + e ^ {\lambda h}}{1 - e ^ {\lambda h}} \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}} \\ a \frac {u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}}{2 h} = a \frac {2 u _ {j} ^ {n} - 2 u _ {j - 1} ^ {n}}{2 h} + a \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{2 h} \\ = a \frac {u _ {j} ^ {n} - u _ {j - 1} ^ {n}}{h} + \frac {a h}{2} \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}} \\ \end{array}
$$

$$
(2) \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j} ^ {n} - u _ {j - 1} ^ {n}}{h} = - \frac {a h}{1 - e ^ {\lambda h}} \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}
$$

考虑迎风格式、指数格式、Samarskii 格式的关系,首先改写指数格式如下:

**注意**：下文及 Taylor 近似中恒有 **$\lambda = a/\nu$（网格 Péclet 数的指数）**，绝非 $\nu/a$。

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + a \frac {u _ {j} ^ {n} - u _ {j - 1} ^ {n}}{h} = - \frac {a h}{1 - e ^ {\lambda h}} \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}
$$

取 $h$ 充分小, $e^{\lambda h} \approx 1 + \lambda h$ , 可以得到迎风格式。

如果 $e^{\lambda h} \approx 1 + \lambda h + \dfrac{(\lambda h)^2}{2}$，则可以得到 Samarskii 格式。于是这两种格式均是指数型格式的特例。

#### 对流占优情形下稳定性条件的比较

迎风格式的稳定性条件是

$$
(\frac {a h}{2} + \nu) \frac {\tau}{h ^ {2}} \leq \frac {1}{2}
$$

Samarskii格式的稳定性条件：

$$
\left[ \frac {a h}{2} + \frac {\nu}{1 + \frac {a h}{2 \nu}} \right] \frac {\tau}{h ^ {2}} \leq \frac {1}{2}
$$

指数型格式稳定性条件为: $\nu \sigma \frac{\tau}{h^{2}} \leq \frac{1}{2}$ 

$$
\sigma = \frac {a h}{2 \nu} \coth \left(\frac {a h}{2 \nu}\right)
$$

### 3.6、隐式格式

为了提高计算精度，放宽稳定性条件的限制，考虑隐式格式。

**Crank–Nicolson 型格式**：

$$
\begin{array}{l} \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} + \frac {a}{2} (\frac {u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}}{2 h} + \frac {u _ {j + 1} ^ {n + 1} - u _ {j - 1} ^ {n + 1}}{2 h}) \\ = \frac {\nu}{2} \left(\frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}} + \frac {u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}}{h ^ {2}}\right) \\ \end{array}
$$

令 $\lambda = \frac{a \tau}{h}, \mu = \frac{\nu \tau}{h^{2}}$ , 格式可改写成:

$$
\begin{array}{l} u _ {j} ^ {n + 1} + \frac {\lambda}{2} \left(\frac {u _ {j + 1} ^ {n + 1} - u _ {j - 1} ^ {n + 1}}{2}\right) - \frac {\mu}{2} \left(u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}\right) \\ = u _ {j} ^ {n} - \frac {\lambda}{2} (\frac {u _ {j + 1} ^ {n} - u _ {j - 1} ^ {n}}{2}) + \frac {\mu}{2} (u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}) \\ \end{array}
$$

将各项在 $\left[\cdot\right]_{j}^{n + \frac{1}{2}}$ 处展开, 知 $E = O\left(\tau^{2} + h^{2}\right)$ , 且

$$
\begin{array}{l} G = \frac {(1 - \mu + \mu \cos k h) - i \frac {\lambda}{2} \sin k h}{(1 + \mu - \mu \cos k h) + i \frac {\lambda}{2} \sin k h} \\ = \frac {\left(1 - 2 \mu \sin^ {2} \frac {k h}{2}\right) - i \frac {\lambda}{2} \sin k h}{\left(1 + 2 \mu \sin^ {2} \frac {k h}{2}\right) + i \frac {\lambda}{2} \sin k h} \\ \end{array}
$$

$$
| G | ^ {2} = \frac {(1 - 2 \mu \sin^ {2} \frac {k h}{2}) ^ {2} + \frac {\lambda^ {2}}{4} \sin^ {2} k h}{(1 + 2 \mu \sin^ {2} \frac {k h}{2}) ^ {2} + \frac {\lambda^ {2}}{4} \sin^ {2} k h} \leq 1
$$

格式是无条件稳定的。