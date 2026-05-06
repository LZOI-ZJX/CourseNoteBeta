# 第四节 变系数方程

变系数初值问题（齐次方程）：

$$
\left\{ \begin{array}{l}
\dfrac{\partial u}{\partial t} - a(x)\dfrac{\partial^2 u}{\partial x^2} = 0 \\
u(x,0) = g(x)
\end{array} \right.
$$

若含源项的非齐次方程，则右端写为 $f(x,t)$，即 $\dfrac{\partial u}{\partial t} - a(x)\dfrac{\partial^2 u}{\partial x^2} = f(x,t)$。

$$
a(x) \geq a_0 > 0,\quad a(x)\ \text{连续}
$$

方程两端在节点 $(j,n)$ 取值，即

$$
\left\{ \begin{array}{l}
\dfrac{\partial u}{\partial t}\Big|_{j}^{n} - a(x)\big|_{j}\,\dfrac{\partial^2 u}{\partial x^2}\Big|_{j}^{n} = 0 \\
u(x_j,0) = g(x_j)
\end{array} \right.
$$

（非齐次时在右端代之以 $f(x_j,t_n)$。）

$$
a(x) \geq a_0 > 0
$$

# 4.1、Taylor级数展开法

由于

$$
\begin{array}{l}
\left[\dfrac{\partial^2 u}{\partial x^2}\right]_{j}^{n}
= \dfrac{u(x_{j+1},t_n)-2u(x_j,t_n)+u(x_{j-1},t_n)}{h^2}+O(h^2) \\
= \left[\dfrac{1}{a(x)}\dfrac{\partial u}{\partial t}\right]_{j}^{n}+O(h^2)
\end{array}
$$

引入时间变量差商，得近似方程的差分格式（以下记网格比 $\lambda := \tau/h^2$，即时间步长 $\tau$ 与空间步长平方 $h^2$ 之比）：

$$
\frac{u_j^{n+1}-u_j^n}{\tau}=a_j\frac{u_{j+1}^n-2u_j^n+u_{j-1}^n}{h^2}
$$

由此得到古典显式格式，其截断误差为 $E = O(\tau + h^2)$。稳定性条件（冻结系数法）为 $\max_x |a(x)|\,\lambda \leq \dfrac{1}{2}$。

$$
\boxed{\frac{u_j^{n+1}-u_j^n}{\tau} = a_j\frac{u_{j+1}^n-2u_j^n+u_{j-1}^n}{h^2}}
$$

古典隐式格式，其截断误差为 $E = O(\tau + h^{2})$，稳定性条件（冻结系数法）：无条件稳定。

$$
\boxed{\frac{u_j^{n+1}-u_j^n}{\tau} = a_j\frac{u_{j+1}^{n+1}-2u_j^{n+1}+u_{j-1}^{n+1}}{h^2}}
$$

若 $a(x)$ 高阶连续可导，可构造更高精度的差分格式，设 $\delta_{x}^{2} u_{j}^{n} = u(x_{j+1}, t_{n}) - 2u(x_{j}, t_{n}) + u(x_{j-1}, t_{n})$

那么：

$$
\begin{array}{l}
\frac{\delta_x^2 u_j^n}{h^2} = \left[\frac{\partial^2 u}{\partial x^2}\right]_{j}^{n} + \frac{h^2}{12}\left[\frac{\partial^4 u}{\partial x^4}\right]_{j}^{n} + O(h^4) \\
= \left(\frac{1}{a(x)}\frac{\partial u}{\partial t}\right)_{j}^{n} + \frac{h^2}{12}\left(\frac{\partial^2}{\partial x^2}\left(\frac{1}{a(x)}\frac{\partial u}{\partial t}\right)\right)_{j}^{n} + O(h^4)
\end{array}
$$

$$
\begin{array}{l}
= \left(\frac{1}{a(x)}\frac{\partial u}{\partial t}\right)_{j}^{n} + \frac{h^2}{12}\left\{\frac{\delta_x^2\left(\frac{1}{a(x)}\frac{\partial u}{\partial t}\right)}{h^2}\Big|_{j}^{n} + O(h^2)\right\} + O(h^4) \\
= \frac{5}{6}\left(\frac{1}{a(x)}\frac{\partial u}{\partial t}\right)_{j}^{n} + \frac{1}{12}\left(\frac{1}{a(x)}\frac{\partial u}{\partial t}\right)_{j+1}^{n} + \frac{1}{12}\left(\frac{1}{a(x)}\frac{\partial u}{\partial t}\right)_{j-1}^{n} + O(h^4)
\end{array}
$$

至此只是消去 $u_{xxxx}$ 后在**固定时刻** $t_n$ 上得到的**空间紧逼近**（右端为时间导数项的组合）。要在时间上达到 $O(\tau^2)$，需用梯形公式把方程在 $[t_n,t_{n+1}]$ 上积分，或对时间**中心差分**，即把时间导数 $u_t$ 在 $t_{n+\frac{1}{2}}$ 近似成 $\dfrac{u^{n+1}-u^n}{\tau}$，同时把上述空间算子取为 $t_n$ 与 $t_{n+1}$ 上的**算术平均**（与 Crank–Nicolson 同一思想）。

引入上述时间离散，得差分格式：

$$
\begin{array}{l}
\frac{1}{12}\frac{u_{j+1}^{n+1}-u_{j+1}^{n}}{a_{j+1}\tau} + \frac{5}{6}\frac{u_j^{n+1}-u_j^n}{a_j\tau} + \frac{1}{12}\frac{u_{j-1}^{n+1}-u_{j-1}^{n}}{a_{j-1}\tau} \\
= \frac{\delta_x^2 u_j^{n+1}+\delta_x^2 u_j^n}{2h^2}
\end{array}
$$

截断误差 $E = O(\tau^{2} + h^{4})$，即为紧格式。

稳定性条件（冻结系数法）：无条件稳定

# 4.2、Keller盒式格式

1970年，H.B Keller提出了解抛物型初边值问题的一个盒式格式，这一格式具有二阶精度，且无条件稳定。

$$
\begin{array}{r}
\text{初边值问题:}\left\{\begin{array}{ll}
\dfrac{\partial u}{\partial t}=\dfrac{\partial}{\partial x}\!\left(a(x)\dfrac{\partial u}{\partial x}\right) & a(x)\geq a_0>0 \\
u(x,0)=g(x) & 0<x<1 \\
u(0,t)=g_0(t) & t>0 \\
u(1,t)=g_1(t) &
\end{array}\right.
\end{array}
$$

首先将方程化为一阶微分方程组：

$$
\left\{\begin{array}{l}
a(x)\dfrac{\partial u}{\partial x}=v(x,t) \qquad (x,t)\in\{(x,t)\mid 0\leq x\leq 1,\ 0\leq t\leq T\} \\
\dfrac{\partial v}{\partial x}=\dfrac{\partial u}{\partial t}
\end{array}\right.
$$

在定解区域上取网格：

$$
R_h(T)=\left\{(x_j,t_n)\,\middle|\,\begin{array}{l} x_j=jh,\ j=0,1,\ldots,J,\ Jh=1 \\ t_n=n\tau,\ n=0,1,\ldots,N,\ N\tau=T \end{array}\right\}
$$

引入记号：

(1) $x_{j\pm\frac{1}{2}}=\frac{1}{2}(x_j+x_{j\pm 1}),\quad t_{n\pm\frac{1}{2}}=\frac{1}{2}(t_n+t_{n\pm 1})$。

(2) $\phi_j^{n\pm\frac{1}{2}}=\frac{1}{2}(\phi_j^n+\phi_j^{n\pm 1}),\quad \phi_{j-\frac{1}{2}}^n=\frac{1}{2}(\phi_j^n+\phi_{j-1}^n)$（**注意**：$u_{j-\frac{1}{2}}^n$ 表示 $x_{j-\frac{1}{2}}$ 处由左右整格点**算术平均**，不是 $u$ 存储在半格点上）。

**Keller 盒式格式（Box）** 的要点：在每一矩形盒 $[x_{j-1},x_j]\times[t_{n-1},t_n]$ 上，把所有待离散的量都看成盒的**四角** $(j{-}1,n{-}1),(j,n{-}1),(j{-}1,n),(j,n)$ 上的取值，对方程在盒的**几何中心** $(x_{j-\frac{1}{2}},t_{n-\frac{1}{2}})$ 离散——$\partial/\partial x$ 用该盒上**左右**两层值的差分，$\partial/\partial t$ 用**上下**两层值的差分，$a u_x=v$ 中的 $v$ 在中心取为**两竖边**上 $v$ 的算术平均。这样两条方程的局部离散算子都落在**同一点** $(j-\frac{1}{2},n-\frac{1}{2})$，与单纯的“交错存数 MAC 网格”或蛙跃等做法在概念上不同。

记 $u_j^n\approx u(x_j,t_n)$，$v_j^{n-\frac{1}{2}}\approx v(x_j,t_{n-\frac{1}{2}})$（辅助变量在**时间半步**、**整空间格点**取值）。盒心处的离散方程组为

$$
\left\{\begin{array}{l}
a_{j-\frac{1}{2}}\dfrac{u_j^{n-\frac{1}{2}}-u_{j-1}^{n-\frac{1}{2}}}{h}
=\dfrac{v_j^{n-\frac{1}{2}}+v_{j-1}^{n-\frac{1}{2}}}{2}, \\[8pt]
\dfrac{v_j^{n-\frac{1}{2}}-v_{j-1}^{n-\frac{1}{2}}}{h}
=\dfrac{u_{j-\frac{1}{2}}^n-u_{j-\frac{1}{2}}^{n-1}}{\tau},
\end{array}\right.
$$

其中  
$u_j^{n-\frac{1}{2}}=\frac{1}{2}(u_j^n+u_j^{n-1})$，  
$u_{j-\frac{1}{2}}^n=\frac{1}{2}(u_j^n+u_{j-1}^n)$（均由记号 (2)）。

第二式正是在 $(j-\frac{1}{2},n-\frac{1}{2})$ 上逼近 $\partial v/\partial x=\partial u/\partial t$：空间步长对应 $x_{j-1}\!\to\!x_j$，时间步长对应 $t_{n-1}\!\to\!t_n$，**右端所用的是 $u$ 在半格点 $x_{j-\frac{1}{2}}$ 上的平均值，而不是 $x_j$ 上的单层值**。

$j=1,2,\ldots,J$；$n=1,2,\ldots,N$。

初始条件：$u_j^0=g(x_j)$。$n\ge1$ 时上式与 $u_j^{n-\frac{1}{2}}=\frac{1}{2}(u_j^n+u_j^{n-1})$ 联立；$n=1$ 时 $u_j^1$ 与 $v_j^{1/2}$ 由该方程组与边界条件一起确定（亦可用与 $au_x=v$ 在 $t=0$ 相容的 $v_j^{1/2}$ 启动，见教材）。

边界条件离散：$\left\{\begin{array}{l} u_0^n=g_0(t_n) \\ u_J^n=g_1(t_n) \end{array}\right.$，$1 \leq n \leq N$

格式特点：

1. 截断误差：二阶精度。
2. 稳定性条件：无条件稳定。
3. 容易编制程序。
4. 可以应用到非线性问题及间断系数问题。

# 4.3、有限体积法

考虑方程：

$$
\frac{\partial u}{\partial t}=\frac{\partial}{\partial x}\!\left(a(x)\frac{\partial u}{\partial x}\right),\qquad a(x)\geq a_0>0
$$

令 $w(x,t)=a(x)\dfrac{\partial u}{\partial x}$，那么方程变成 $\dfrac{\partial u}{\partial t}=\dfrac{\partial w}{\partial x}$。

对上面的方程在 $[x_{j-\frac{1}{2}}, x_{j+\frac{1}{2}}]\times[t_n,t_{n+1}]$ 上积分有

$$
\begin{array}{l}
\displaystyle\int_{x_{j-\frac{1}{2}}}^{x_{j+\frac{1}{2}}}\bigl[u(x,t_{n+1})-u(x,t_n)\bigr]\,dx \\
\displaystyle =\int_{t_n}^{t_{n+1}}\bigl[w(x_{j+\frac{1}{2}},t)-w(x_{j-\frac{1}{2}},t)\bigr]\,dt
\end{array}
$$

据 $w(x,t)$ 的定义有 $\dfrac{\partial u}{\partial x}=\dfrac{w(x,t)}{a(x)}$，在 $[x_j,x_{j+1}]$ 上积分有

$$
u(x_{j+1},t)-u(x_j,t)=\int_{x_j}^{x_{j+1}}\frac{w(x,t)}{a(x)}\,dx
$$

$$
\approx w\!\left(x_{j+\frac{1}{2}},t\right)\int_{x_j}^{x_{j+1}}\frac{1}{a(x)}\,dx
$$

$$
w\!\left(x_{j+\frac{1}{2}},t\right)\approx \frac{h}{\displaystyle\int_{x_j}^{x_{j+1}}\frac{1}{a(x)}\,dx}\cdot\frac{u(x_{j+1},t)-u(x_j,t)}{h}
$$

令

$$
A_{j+1}=h\left(\int_{x_j}^{x_{j+1}}\frac{1}{a(x)}\,dx\right)^{-1}
$$

则有

$$
w\!\left(x_{j+\frac{1}{2}},t\right)\approx A_{j+1}\cdot\frac{u(x_{j+1},t)-u(x_j,t)}{h}
$$

积分近似：

$$
\int_{x_{j-\frac{1}{2}}}^{x_{j+\frac{1}{2}}}u(x,t_n)\,dx\approx h\,u(x_j,t_n),
$$

$$
\int_{t_n}^{t_{n+1}}w(x_{j+\frac{1}{2}},t)\,dt\approx \tau\bigl[\theta w(x_{j+\frac{1}{2}},t_{n+1})+(1-\theta)w(x_{j+\frac{1}{2}},t_n)\bigr]
$$

**于是有**（时间层与区间 $[t_n,t_{n+1}]$ 上的积分一致：$n\to n{+}1$）

$$
\begin{array}{l}
\displaystyle\frac{u(x_j,t_{n+1})-u(x_j,t_n)}{\tau} \\
\displaystyle =\frac{1}{h^2}\bigl[\theta\,\Delta(A_j\nabla u(x,t_{n+1}))+(1-\theta)\,\Delta(A_j\nabla u(x,t_n))\bigr]+R_j^n
\end{array}
$$

$\Delta$、$\nabla$ 分别表示向前和向后差分，即

$$
\begin{array}{l}
\Delta u(x_j,t_n)=u(x_{j+1},t_n)-u(x_j,t_n) \\
\nabla u(x_j,t_n)=u(x_j,t_n)-u(x_{j-1},t_n)
\end{array}
$$

**忽略误差项，得差分格式**

$$
\frac{u_j^{n+1}-u_j^n}{\tau}-\frac{1}{h^2}\bigl[\theta\,\Delta(A_j\nabla u_j^{n+1})+(1-\theta)\,\Delta(A_j\nabla u_j^n)\bigr]=0
$$

$A_j$ 的计算：

1. 数值积分的误差阶不低于差分格式的截断误差。
2. 由定义，$A_j=h\bigl(\int_{x_{j-1}}^{x_j}\frac{1}{a(x)}\,dx\bigr)^{-1}$ 是相邻控制体界面上的**调和平均**型传导系数（对 $\frac{1}{a}$ 积分后取倒数），**不能**用算术平均 $\frac{1}{2}\bigl[a(x_{j-1})+a(x_j)\bigr]$ 代替——在变系数或间断附近会严重扭曲通量。对 $\frac{1}{a(x)}$ 在 $[x_{j-1},x_j]$ 上用梯形公式，即得常用的两点调和平均近似（见下式）。另一种常见选择是 $A_j\approx a\!\left(x_{j-\frac{1}{2}}\right)$（对 $\frac{1}{a}$ 的中点法则）；它与上式都是界面系数的积分逼近，但与**算术平均**不是同一回事。

$$
A_j\approx\frac{2a(x_{j-1})a(x_j)}{a(x_{j-1})+a(x_j)}.
$$

**格式的特点**

1. 截断误差：$\theta=\dfrac{1}{2}$ 时 $R_j^n=O(\tau^2+h^2)$；$\theta\neq\dfrac{1}{2}$ 时 $R_j^n=O(\tau+h^2)$。
2. 稳定性：$\dfrac{1}{2}\leq\theta\leq 1$ 时无条件稳定；$0\leq\theta<\dfrac{1}{2}$ 时条件稳定。

实际上：$\theta=0$ 为显式格式，$\theta=1$ 为隐式格式；$\theta=\dfrac{1}{2}$ 为 Crank–Nicolson 格式。

# 4.4、间断系数问题

考虑方程

$$
\frac{\partial u}{\partial t}-\frac{\partial}{\partial x}\!\left(a(x)\frac{\partial u}{\partial x}\right)+b(x)\,u=0,\qquad a(x)\geq a_0>0
$$

假定在 $x=\xi$ 处，$a(x)$、$b(x)$ 有第一类间断；在间断线上解满足如下联结条件

$$
\left\{\begin{array}{l}
u(\xi-0,t)=u(\xi+0,t) \\
a(\xi-0)\dfrac{\partial u}{\partial x}(\xi-0,t)=a(\xi+0)\dfrac{\partial u}{\partial x}(\xi+0,t)
\end{array}\right.
$$

为了处理问题的方便，通常使间断线与网格线重合。

即 $x_j=\xi$。

采用有限体积方法：

积分区域：$D=\{(x,t)\mid x_{j-\frac{1}{2}}\leq x\leq x_{j+\frac{1}{2}},\ t_n\leq t\leq t_{n+1}\}$

对上面的方程在 $D$ 上积分有

$$
\begin{array}{l}
\displaystyle\int_{x_{j-\frac{1}{2}}}^{x_{j+\frac{1}{2}}}\int_{t_n}^{t_{n+1}}\frac{\partial u(x,t)}{\partial t}\,dt\,dx \\
\displaystyle\quad=\int_{t_n}^{t_{n+1}}\int_{x_{j-\frac{1}{2}}}^{x_{j+\frac{1}{2}}}\frac{\partial w(x,t)}{\partial x}\,dx\,dt
-\int_{t_n}^{t_{n+1}}\int_{x_{j-\frac{1}{2}}}^{x_{j+\frac{1}{2}}}b(x)\,u(x,t)\,dx\,dt .
\end{array}
$$

左端对时间的积分化为 $\displaystyle\int_{x_{j-\frac{1}{2}}}^{x_{j+\frac{1}{2}}}\bigl[u(x,t_{n+1})-u(x,t_n)\bigr]\,dx$。对通量项，因 $x_j$ 为间断点，将 $\displaystyle\int_{x_{j-\frac{1}{2}}}^{x_{j+\frac{1}{2}}}\frac{\partial w}{\partial x}\,dx$ 写成相邻两个子区间积分之和，即左半单元积分 **+** 右半单元积分，有

$$
\begin{array}{l}
\displaystyle\int_{t_n}^{t_{n+1}}\int_{x_{j-\frac{1}{2}}}^{x_{j+\frac{1}{2}}}\frac{\partial w(x,t)}{\partial x}\,dx\,dt \\
\displaystyle\quad=\int_{t_n}^{t_{n+1}}\left(\int_{x_{j-\frac{1}{2}}}^{x_j}\frac{\partial w(x,t)}{\partial x}\,dx+\int_{x_j}^{x_{j+\frac{1}{2}}}\frac{\partial w(x,t)}{\partial x}\,dx\right)dt \\[0.5em]
\displaystyle\quad=\int_{t_n}^{t_{n+1}}\bigl[w(x_{j+\frac{1}{2}},t)-w(x_j+0,t)+w(x_j-0,t)-w(x_{j-\frac{1}{2}},t)\bigr]\,dt \\[0.5em]
\displaystyle\quad=\int_{t_n}^{t_{n+1}}\bigl[w(x_{j+\frac{1}{2}},t)-w(x_{j-\frac{1}{2}},t)\bigr]\,dt ,
\end{array}
$$

其中最后一步用到联结条件所蕴含的通量连续：$w(x_j-0,t)=w(x_j+0,t)$。以下记单侧极限的简写 $b_j^-:=b(x_j-0)$，$b_j^+:=b(x_j+0)$。对含 $b(x)\,u$ 的源项积分，同样在 $x_j$ 处拆分区间，并写清相加关系：

$$
\begin{array}{l}
\displaystyle\int_{t_n}^{t_{n+1}}\int_{x_{j-\frac{1}{2}}}^{x_{j+\frac{1}{2}}}b(x)\,u(x,t)\,dx\,dt
=\int_{t_n}^{t_{n+1}}\left(\int_{x_{j-\frac{1}{2}}}^{x_j}b(x)\,u(x,t)\,dx+\int_{x_j}^{x_{j+\frac{1}{2}}}b(x)\,u(x,t)\,dx\right)dt \\[0.6em]
\displaystyle\approx\int_{t_n}^{t_{n+1}}\frac{h}{2}(b_j^-+b_j^+)\,u(x_j,t)\,dt \\[0.6em]
\displaystyle\approx\frac{\tau h}{2}(b_j^-+b_j^+)\,u(x_j,t_n) .
\end{array}
$$

**可得差分格式**

$$
\frac{u_j^{n+1}-u_j^n}{\tau}-\frac{1}{h^2}\Delta(A_j\nabla u_j^n)+\frac{1}{2}(b_j^-+b_j^+)u_j^n=0
$$

$$
A_j=h\left(\int_{x_{j-1}}^{x_j}\frac{1}{a(x)}\,dx\right)^{-1}
$$

注：由此看出有限体积法在处理间断系数问题的有效性。对于右端项有间断的情形同样适用。
