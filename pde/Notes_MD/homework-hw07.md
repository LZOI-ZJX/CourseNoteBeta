## 题目 1：扩散方程差分格式的截断误差与精度调整

**【题目描述】**
试讨论扩散方程 $\frac{\partial u}{\partial t} = \frac{\partial^2 u}{\partial x^2}$ 的差分格式：
$$\theta \frac{u_j^{n+1} - u_j^{n-1}}{2\tau} + (1 - \theta) \frac{u_j^n - u_j^{n-1}}{\tau} = \frac{1}{h^2}\delta_x^2 u_j^n, \quad 0 \leqslant \theta \leqslant 1$$
的截断误差，并调整 $\theta$ 使差分格式为二阶精度的。

**【解答过程】**
设 $u(x, t)$ 为原扩散方程的充分光滑的精确解。将其代入差分格式，定义截断误差 $T_j^n$ 为：
$$T_j^n = \theta \frac{u(x_j, t_{n+1}) - u(x_j, t_{n-1})}{2\tau} + (1 - \theta) \frac{u(x_j, t_n) - u(x_j, t_{n-1})}{\tau} - \frac{1}{h^2} \delta_x^2 u(x_j, t_n)$$

在节点 $(x_j, t_n)$ 处对各项进行泰勒（Taylor）展开：

**1. 时间项展开：**
$$u(x_j, t_{n+1}) = u + \tau u_t + \frac{\tau^2}{2} u_{tt} + \frac{\tau^3}{6} u_{ttt} + O(\tau^4)$$
$$u(x_j, t_{n-1}) = u - \tau u_t + \frac{\tau^2}{2} u_{tt} - \frac{\tau^3}{6} u_{ttt} + O(\tau^4)$$

分别代入中心差商和向后差商：
中心差商：
$$\frac{u(x_j, t_{n+1}) - u(x_j, t_{n-1})}{2\tau} = u_t + \frac{\tau^2}{6} u_{ttt} + O(\tau^4)$$
向后差商：
$$\frac{u(x_j, t_n) - u(x_j, t_{n-1})}{\tau} = u_t - \frac{\tau}{2} u_{tt} + \frac{\tau^2}{6} u_{ttt} + O(\tau^3)$$

**2. 空间项展开：**
中心二阶差分算子 $\delta_x^2 u(x_j, t_n) = u(x_{j+1}, t_n) - 2u(x_j, t_n) + u(x_{j-1}, t_n)$，展开得：
$$\frac{1}{h^2} \delta_x^2 u(x_j, t_n) = u_{xx} + \frac{h^2}{12} u_{xxxx} + O(h^4)$$

**3. 合并截断误差：**
代回 $T_j^n$ 的表达式：
$$T_j^n = \theta \left( u_t + \frac{\tau^2}{6} u_{ttt} + \dots \right) + (1 - \theta) \left( u_t - \frac{\tau}{2} u_{tt} + \dots \right) - \left( u_{xx} + \frac{h^2}{12} u_{xxxx} + \dots \right)$$

化简并合并时间导数部分：
$$T_j^n = u_t - \frac{\tau}{2}(1 - \theta) u_{tt} + \frac{\tau^2}{6} u_{ttt} - u_{xx} - \frac{h^2}{12} u_{xxxx} + O(\tau^3 + h^4)$$

由于 $u(x,t)$ 是精确解，满足 $u_t = u_{xx}$，两者抵消。最终截断误差为：
$$T_j^n = - \frac{\tau}{2}(1 - \theta) u_{tt} - \frac{h^2}{12} u_{xxxx} + O(\tau^2 + h^4)$$

**4. 调整 $\theta$ 以达到二阶精度：**
为了使格式在时间上达到二阶精度（即截断误差为 $O(\tau^2 + h^2)$），必须消除 $O(\tau)$ 的一阶项。
令一阶项系数为零：
$$- \frac{1}{2}(1 - \theta) u_{tt} = 0$$
解得：**$\theta = 1$**。

当 $\theta = 1$ 时，差分格式变为 Leapfrog 中心差分格式，截断误差为 $O(\tau^2 + h^2)$。

---

## 题目 2：变系数抛物型方程的二阶精度差分格式

**【题目描述】**
试构造变系数抛物型方程：
$$\frac{\partial u}{\partial t} = \frac{\partial}{\partial x} \left[ (0.1 + \sin^2 x) \frac{\partial u}{\partial x} \right]$$
的一个二阶精度差分格式。

**【解答过程】**
为了构造一个时空均具有二阶精度（即截断误差为 $O(\tau^2 + h^2)$）的差分格式，采用 **Crank-Nicolson 格式** 结合 **空间半步长中心差分**。

记变系数 $a(x) = 0.1 + \sin^2 x$。设空间步长为 $h$，时间步长为 $\tau$。

**1. 空间导数的二阶中心差分近似：**
在固定时间层上，对空间项 $\frac{\partial}{\partial x} \left( a(x) \frac{\partial u}{\partial x} \right)$ 使用半节点 $(x_{j+1/2})$ 和 $(x_{j-1/2})$ 进行两次中心差分。

$$\left. \frac{\partial}{\partial x} \left( a(x) \frac{\partial u}{\partial x} \right) \right|_{x_j} \approx \frac{1}{h} \left[ a_{j+1/2} \left( \frac{\partial u}{\partial x} \right)_{j+1/2} - a_{j-1/2} \left( \frac{\partial u}{\partial x} \right)_{j-1/2} \right]$$

展开一阶导数得空间差分算子 $L_h u_j$：
$$L_h u_j = \frac{1}{h^2} \left[ a_{j+1/2} (u_{j+1} - u_j) - a_{j-1/2} (u_j - u_{j-1}) \right]$$
即：
$$L_h u_j = \frac{1}{h^2} \left[ a_{j+1/2} u_{j+1} - (a_{j+1/2} + a_{j-1/2}) u_j + a_{j-1/2} u_{j-1} \right]$$
其中，$a_{j\pm1/2} = 0.1 + \sin^2 \left( (j \pm \frac{1}{2})h \right)$。

**2. 时间导数的 Crank-Nicolson 离散：**
为了在时间上达到二阶精度 $O(\tau^2)$，在半时间层 $t_{n+1/2}$ 处展开。
左端时间导数使用中心差分：
$$\left. \frac{\partial u}{\partial t} \right|_{t_{n+1/2}} \approx \frac{u_j^{n+1} - u_j^n}{\tau}$$

右端空间导数取第 $n$ 层和第 $n+1$ 层的算术平均值：
$$\left. \frac{\partial}{\partial x} \left( a(x) \frac{\partial u}{\partial x} \right) \right|_{t_{n+1/2}} \approx \frac{1}{2} (L_h u_j^{n+1} + L_h u_j^n)$$

**3. 最终格式组装：**
将两者结合，得到完整的 Crank-Nicolson 差分格式：
$$\frac{u_j^{n+1} - u_j^n}{\tau} = \frac{1}{2h^2} \left[ a_{j+1/2} (u_{j+1}^{n+1} + u_{j+1}^n) - (a_{j+1/2} + a_{j-1/2}) (u_j^{n+1} + u_j^n) + a_{j-1/2} (u_{j-1}^{n+1} + u_{j-1}^n) \right]$$

该格式具有 $O(\tau^2 + h^2)$ 的截断误差，达到了题目要求的二阶精度，且是无条件稳定的。