# 第五节 多维问题

# 二维问题

初边值问题：

$$
\left\{ \begin{array}{l l} \frac {\partial u}{\partial t} = a \left(\frac {\partial^ {2} u}{\partial x ^ {2}} + \frac {\partial^ {2} u}{\partial y ^ {2}}\right) & 0 <   x, y <   1 \\ u (x,y,0)=f(x,y),\quad 0<x,y<1, \\ u (0, y, t) = u (1, y, t)=u (x, 0, t) = u (x, 1, t) = 0, & t\ge 0 \end{array} \right.
$$

取 $\Delta t = \tau, \Delta x = \Delta y = h$ , 形成网格:

$$
D _ {h} = \left\{ \begin{array}{l l} & x _ {j} = j h, \quad j = 0, 1 \dots J \quad J h = 1 \\ (x _ {j}, y _ {l}, t _ {n}) & y _ {l} = l h, \quad l = 0, 1 \dots J \quad J h = 1 \\ & t _ {n} = n \tau \quad n \geq 0 \end{array} \right\}
$$

并记: $\delta_{x}^{2} u_{j, l}^{n} = u_{j+1, l}^{n} - 2u_{j, l}^{n} + u_{j-1, l}^{n}$ ,

$$
\delta_ {y} ^ {2} u _ {j, l} ^ {n} = u _ {j, l + 1} ^ {n} - 2 u _ {j, l} ^ {n} + u _ {j, l - 1} ^ {n}
$$

# 5.1、一维格式的直接推广

一维扩散方程古典显式格式推广：

$$
\frac {u _ {j l} ^ {n + 1} - u _ {j l} ^ {n}}{\tau} = a \frac {1}{h ^ {2}} \left(\delta_ {x} ^ {2} u _ {j l} ^ {n} + \delta_ {y} ^ {2} u _ {j l} ^ {n}\right) \tag {5.1}
$$

利用 Taylor 展开式分析截断误差为 $O\left(\tau + h^{2}\right)$ 

计算格式

$$
u _ {j l} ^ {n + 1} = u _ {j l} ^ {n} + a \frac {\tau}{h ^ {2}} \left(\delta_ {x} ^ {2} u _ {j l} ^ {n} + \delta_ {y} ^ {2} u _ {j l} ^ {n}\right) \tag {5.2}
$$

可用Fourier法分析稳定性：

设 $\lambda = \frac{\tau}{h^{2}}, u_{jl}^{n} = v^{n} e^{ik_{1}jh} e^{ik_{2}lh}$ , 代入(5.2), 有

$$
v ^ {n + 1} = (1 + 2 a \lambda (\cos k _ {1} h - 1) + 2 a \lambda (\cos k _ {2} h - 1)) v ^ {n}
$$

传播因子: $G\left(\lambda , k_{1}, k_{2}\right) = 1 - 4 a \lambda \left(\sin^{2} \frac{k_{1} h}{2} + \sin^{2} \frac{k_{2} h}{2}\right)$ 

解得当 $a \lambda \leq \frac{1}{4}$ 时, $|G| \leq 1$ , 格式稳定。

一般 $p$ 维显式格式稳定的条件是 $a \lambda \leq \frac{1}{2 p}$ 。

一维古典隐式格式的推广

$$
\frac {u _ {j l} ^ {n} - u _ {j l} ^ {n - 1}}{\tau} = \frac {a}{h ^ {2}} \left(\delta_ {x} ^ {2} u _ {j l} ^ {n} + \delta_ {y} ^ {2} u _ {j l} ^ {n}\right) \tag {5.3}
$$

格式的截断误差仍为 $O\left(\tau + h^{2}\right)$ , 但增长因子：

$$
G (\lambda , k _ {1}, k _ {2}) = \frac {1}{1 + 4 a \lambda \left(\sin^ {2} \frac {k _ {1} h}{2} + \sin^ {2} \frac {k _ {2} h}{2}\right)}
$$

$\left| \mathrm{G}(\lambda, k) \right| \leq 1.$ 格式无条件稳定。

Crank - Nicolson格式：

$$
\frac {u _ {j l} ^ {n + 1} - u _ {j l} ^ {n}}{\tau} = \frac {a}{2 h ^ {2}} \left[ \delta_ {x} ^ {2} \left(u _ {j l} ^ {n + 1} + u _ {j l} ^ {n}\right) + \delta_ {y} ^ {2} \left(u _ {j l} ^ {n + 1} + u _ {j l} ^ {n}\right) \right] \tag {5.4}
$$

格式具有二阶精度, $E = O\left(\tau^{2} + h^{2}\right)$ 。

$$
G (\lambda , k _ {1}, k _ {2}) = \frac {1 - 2 a \lambda \sin^ {2} \frac {k _ {1} h}{2} - 2 a \lambda \sin^ {2} \frac {k _ {2} h}{2}}{1 + 2 a \lambda \sin^ {2} \frac {k _ {1} h}{2} + 2 a \lambda \sin^ {2} \frac {k _ {2} h}{2}}
$$

对于任何 $\lambda, k_{1,2}, |G(\lambda, k_1, k_2)| \leq 1$ , 格式绝对稳定。

# 5.2、交替方向隐式格式

# 1: PR（Peaceman-Rachford）格式

$$
\left[ \frac {u _ {j l} ^ {n + \frac {1}{2}} - u _ {j l} ^ {n}}{\frac {\tau}{2}} = a \frac {1}{h ^ {2}} \left(\delta_ {x} ^ {2} u _ {j l} ^ {n + \frac {1}{2}} + \delta_ {y} ^ {2} u _ {j l} ^ {n}\right) \right.
$$

$$
\frac {u _ {j l} ^ {n + 1} - u _ {j l} ^ {n + \frac {1}{2}}}{\frac {\tau}{2}} = a \frac {1}{h ^ {2}} \left(\delta_ {x} ^ {2} u _ {j l} ^ {n + \frac {1}{2}} + \delta_ {y} ^ {2} u _ {j l} ^ {n + 1}\right)
$$

注：求解时每层需用两次追赶法求解

(5.5) 

# 格式的特点

# 1、截断误差

消去 $u_{jl}^{n + \frac{1}{2}}, \delta_x^2 u_{jl}^{n + \frac{1}{2}}$ , 可得

$$
\begin{array}{l} (1 + \frac {1}{4} \frac {\tau^ {2} a ^ {2}}{h ^ {4}} \delta_ {x} ^ {2} \delta_ {y} ^ {2}) \frac {u _ {j l} ^ {n + 1} - u _ {j l} ^ {n}}{\tau} \\ = \frac {a}{h ^ {2}} \left(\delta_ {x} ^ {2} + \delta_ {y} ^ {2}\right) \frac {u _ {j l} ^ {n + 1} + u _ {j l} ^ {n}}{2} \tag {5.6} \\ \end{array}
$$

Taylor展开分析截断误差: $O\left(\tau^{2} + h^{2}\right)$ 

所以PR格式是二阶精度的格式

# 2、稳定性

改写（5.6）为

$$
\begin{array}{l} \big( 1 - \frac {a \lambda}{2} \delta_ {x} ^ {2} \big) \big( 1 - \frac {a \lambda}{2} \delta_ {y} ^ {2} \big) \, u _ {j l} ^ {n + 1} \\ = \big( 1 + \frac {a \lambda}{2} \delta_ {x} ^ {2} \big) \big( 1 + \frac {a \lambda}{2} \delta_ {y} ^ {2} \big) \, u _ {j l} ^ {n} \tag {5.7} \\ \end{array}
$$

Fourier方法分析稳定性，得增长因子为

$$
G (\lambda , k _ {1}, k _ {2}) = \frac {(1 - 2 a \lambda \sin^ {2} \frac {k _ {1} h}{2}) (1 - 2 a \lambda \sin^ {2} \frac {k _ {2} h}{2})}{(1 + 2 a \lambda \sin^ {2} \frac {k _ {1} h}{2}) (1 + 2 a \lambda \sin^ {2} \frac {k _ {2} h}{2})}
$$

所以PR格式是无条件稳定的。

# 附、PR格式的另一种推导方法

也可以通过 Crank—Nicolson 格式变形的途径得到 Peaceman-Rachford 格式，首先将 Crank—Nicolson 格式变形为

$$
u _ {j l} ^ {n + 1} - \frac {a \lambda}{2} \Big [ \delta_ {x} ^ {2} u _ {j l} ^ {n + 1} + \delta_ {y} ^ {2} u _ {j l} ^ {n + 1} \Big ] = u _ {j l} ^ {n} + \frac {a \lambda}{2} \Big [ \delta_ {x} ^ {2} u _ {j l} ^ {n} + \delta_ {y} ^ {2} u _ {j l} ^ {n} \Big ]
$$

$$
\frac {a ^ {2} \lambda^ {2}}{4} \delta_ {x} ^ {2} \delta_ {y} ^ {2} (u _ {j l} ^ {n + 1} - u _ {j l} ^ {n}) = \frac {a ^ {2} \lambda^ {2}}{4} \delta_ {x} ^ {2} \delta_ {y} ^ {2} [ \delta_ {t} u _ {j l} ^ {n + \frac {1}{2}} + O (\tau^ {2}) ]
$$

$$
\begin{array}{l} \boldsymbol {u} _ {j l} ^ {n + 1} - \frac {a \lambda}{2} \Big [ \delta_ {x} ^ {2} \boldsymbol {u} _ {j l} ^ {n + 1} + \delta_ {y} ^ {2} \boldsymbol {u} _ {j l} ^ {n + 1} \Big ] + \frac {a ^ {2} \lambda^ {2}}{4} \delta_ {x} ^ {2} \delta_ {y} ^ {2} \boldsymbol {u} _ {j l} ^ {n + 1} \\ = u _ {j l} ^ {n} + \frac {a \lambda}{2} \left[ \delta_ {x} ^ {2} u _ {j l} ^ {n} + \delta_ {y} ^ {2} u _ {j l} ^ {n} \right] + \frac {a ^ {2} \lambda^ {2}}{4} \delta_ {x} ^ {2} \delta_ {y} ^ {2} \left[ u _ {j l} ^ {n} + \delta_ {t} u _ {j l} ^ {n + \frac {1}{2}} + O (\tau^ {2}) \right] \\ \end{array}
$$

忽略高阶无穷小项，并因式分解得到

$$
(1 - \frac {a \lambda}{2} \delta_ {x} ^ {2}) (1 - \frac {a \lambda}{2} \delta_ {y} ^ {2}) u _ {j l} ^ {n + 1} = (1 + \frac {a \lambda}{2} \delta_ {x} ^ {2}) (1 + \frac {a \lambda}{2} \delta_ {y} ^ {2}) u _ {j l} ^ {n}
$$

令 $(1 - \frac{a\lambda}{2}\delta_x^2)u_{jl}^{n + \frac{1}{2}} = (1 + \frac{a\lambda}{2}\delta_y^2)u_{jl}^n$ 

则 $(1 - \frac{a\lambda}{2}\delta_y^2)u_{jl}^{n + 1} = (1 + \frac{a\lambda}{2}\delta_x^2)u_{jl}^{n + \frac{1}{2}}$ 

所以有: $\left\{ \begin{array}{l} (1 - \frac{a \lambda}{2} \delta_{x}^{2}) u_{j l}^{n + \frac{1}{2}} = (1 + \frac{a \lambda}{2} \delta_{y}^{2}) u_{j l}^{n} \\ (1 - \frac{a \lambda}{2} \delta_{y}^{2}) u_{j l}^{n + 1} = (1 + \frac{a \lambda}{2} \delta_{x}^{2}) u_{j l}^{n + \frac{1}{2}} \end{array} \right.$ 

注：PR格式与CN格式有相同的精度，但计算量小很多。

# 实际上还可以构造其他形式的ADI

$$
\left\{ \begin{array}{l} (1 - \frac {a \lambda}{2} \delta_ {x} ^ {2}) u _ {j l} ^ {n + \frac {1}{2}} = (1 + \frac {a \lambda}{2} \delta_ {x} ^ {2}) (1 + \frac {a \lambda}{2} \delta_ {y} ^ {2}) u _ {j l} ^ {n} \\ (1 - \frac {a \lambda}{2} \delta_ {y} ^ {2}) u _ {j l} ^ {n + 1} = u _ {j l} ^ {n + \frac {1}{2}} \end{array} \right.
$$

格式特点：1、二阶精度

2、无条件稳定

# 2: Douglas格式

Peaceman-Rachford 格式的缺点是无法推广到三维的情况。因为它不能象二维格式那样具有对称形式的传播因子, 无条件稳定性不成立。下面的格式可以推广到三维的情况。

改写（5.6）为

$$
\begin{array}{l} (1 - \frac {a \lambda}{2} \delta_ {x} ^ {2}) (1 - \frac {a \lambda}{2} \delta_ {y} ^ {2}) \frac {u _ {j l} ^ {n + 1} - u _ {j l} ^ {n}}{\tau} \\ = \frac {a}{h ^ {2}} \left(\delta_ {x} ^ {2} + \delta_ {y} ^ {2}\right) u _ {j l} ^ {n} \tag {5.8} \\ \end{array}
$$

把(5.8)分裂为

$$
\left\{ \begin{array}{l} (1 - \frac {a \lambda}{2} \delta_ {x} ^ {2}) \frac {\boldsymbol {u} _ {j l} ^ {n + \frac {1}{2}} - \boldsymbol {u} _ {j l} ^ {n}}{\tau} = \frac {a}{h ^ {2}} (\delta_ {x} ^ {2} + \delta_ {y} ^ {2}) \boldsymbol {u} _ {j l} ^ {n} \\ (1 - \frac {a \lambda}{2} \delta_ {y} ^ {2}) \frac {\boldsymbol {u} _ {j l} ^ {n + 1} - \boldsymbol {u} _ {j l} ^ {n}}{\tau} = \frac {\boldsymbol {u} _ {j l} ^ {n + \frac {1}{2}} - \boldsymbol {u} _ {j l} ^ {n}}{\tau} \end{array} \right.
$$

还可以表示为:

$$
\left\{ \begin{array}{l} (1 - \frac {a \lambda}{2} \delta_ {x} ^ {2}) u _ {j l} ^ {n + \frac {1}{2}} - u _ {j l} ^ {n} = a \lambda (\delta_ {x} ^ {2} + \delta_ {y} ^ {2}) u _ {j l} ^ {n} \\ u _ {j l} ^ {n + 1} - u _ {j l} ^ {n + \frac {1}{2}} = \frac {a}{2} \lambda \delta_ {y} ^ {2} (u _ {j l} ^ {n + 1} - u _ {j l} ^ {n}) \end{array} \right.
$$

此格式称为Douglas格式，它的稳定性、局部截断误差和PR格式相同，易于向高维推广。

# 三维的Douglas格式

$$
\left\{ \begin{array}{c} \bigl( 1 - \frac {a \lambda}{2} \delta_ {x} ^ {2} \bigr) \bigl( u _ {i j l} ^ {n + \frac {1}{3}} - u _ {i j l} ^ {n} \bigr) = a \lambda \bigl( \delta_ {x} ^ {2} + \delta_ {y} ^ {2} + \delta_ {z} ^ {2} \bigr) u _ {i j l} ^ {n} \\ u _ {i j l} ^ {n + \frac {2}{3}} - u _ {i j l} ^ {n + \frac {1}{3}} = \frac {a \lambda}{2} \delta_ {y} ^ {2} \bigl( u _ {i j l} ^ {n + \frac {2}{3}} - u _ {i j l} ^ {n} \bigr) \\ u _ {i j l} ^ {n + 1} - u _ {i j l} ^ {n + \frac {2}{3}} = \frac {a \lambda}{2} \delta_ {z} ^ {2} \bigl( u _ {i j l} ^ {n + 1} - u _ {i j l} ^ {n} \bigr) \end{array} \right.
$$

同书（5.35）

# 5.3、局部一维格式

对于方程建立如下的局部一维格式:

$$
\left\{ \begin{array}{l} \frac {1}{2} \frac {u _ {j l} ^ {n + \frac {1}{2}} - u _ {j l} ^ {n}}{\frac {\tau}{2}} = \frac {a}{h ^ {2}} \delta_ {x} ^ {2} \left(\frac {u _ {j l} ^ {n + \frac {1}{2}} + u _ {j l} ^ {n}}{2}\right) \\ \frac {1}{2} \frac {u _ {j l} ^ {n + 1} - u _ {j l} ^ {n + \frac {1}{2}}}{\frac {\tau}{2}} = \frac {a}{h ^ {2}} \delta_ {y} ^ {2} \left(\frac {u _ {j l} ^ {n + 1} + u _ {j l} ^ {n + \frac {1}{2}}}{2}\right) \end{array} \right.
$$

计算时，可化成：

$$
\left\{ \begin{array}{l} \left(1 - \frac {a}{2} \lambda \delta_ {x} ^ {2}\right) u _ {j l} ^ {n + \frac {1}{2}} = \left(1 + \frac {a}{2} \lambda \delta_ {x} ^ {2}\right) u _ {j l} ^ {n} \\ \left(1 - \frac {a}{2} \lambda \delta_ {y} ^ {2}\right) u _ {j l} ^ {n + 1} = \left(1 + \frac {a}{2} \lambda \delta_ {y} ^ {2}\right) u _ {j l} ^ {n + \frac {1}{2}} \end{array} \right. \tag {6.6}
$$

用2次追赶法即可求得 $\left\{\boldsymbol{u}_{jl}^{n+1}\right\}$ 。

此格式与 $P - R$ 格式等价, 是二阶绝对稳定格式。且可推广至三维。

# 5.4、预测-校正格式

$$
\begin{array}{l} \left[ \frac {u _ {j l} ^ {n + \frac {1}{4}} - u _ {j l} ^ {n}}{\frac {\tau}{2}} = \frac {a}{h ^ {2}} \delta_ {x} ^ {2} u _ {j l} ^ {n + \frac {1}{4}} \right. \\ \frac {\boldsymbol {u} _ {j l} ^ {n + \frac {1}{2}} - \boldsymbol {u} _ {j l} ^ {n + \frac {1}{4}}}{\frac {\tau}{2}} = \frac {\boldsymbol {a}}{h ^ {2}} \delta_ {y} ^ {2} \boldsymbol {u} _ {j l} ^ {n + \frac {1}{2}} \\ \frac {\boldsymbol {u} _ {j l} ^ {n + 1} - \boldsymbol {u} _ {j l} ^ {n}}{\tau} = \frac {\boldsymbol {a}}{h ^ {2}} \left(\delta_ {x} ^ {2} + \delta_ {y} ^ {2}\right) \boldsymbol {u} _ {j l} ^ {n + \frac {1}{2}} \tag {6.7} \\ \end{array}
$$

消去 $u^{n + \frac{1}{4}}$ 与 $u^{n + \frac{1}{2}}$ , 可得(5.8)

所以此格式与 $P - R$ 格式等价, 是二阶绝对稳定格式
