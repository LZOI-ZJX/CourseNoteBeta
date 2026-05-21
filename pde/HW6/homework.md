# 题目与解答：Strang 分裂法二阶精度推导

## 一、 题目背景与问题提出

**考虑最简单的一阶常系数方程组 (6.9)：**
$$\frac{\partial \boldsymbol{u}}{\partial t} + \boldsymbol{A} \frac{\partial \boldsymbol{u}}{\partial x} + \boldsymbol{B} \frac{\partial \boldsymbol{u}}{\partial y} = \mathbf{0}$$
其中 $\boldsymbol{u} = [u_1, \cdots, u_p]^\mathrm{T}$，$\boldsymbol{A}, \boldsymbol{B}$ 为实的 $p$ 阶方阵。

**已知一维 Lax-Wendroff 格式的差分算子逼近相应的微分算子：**
* $x$ 方向：$L^x(\tau) = \boldsymbol{I} - \tau \boldsymbol{A} \frac{\partial}{\partial x} + \frac{1}{2}\tau^2 \boldsymbol{A}^2 \frac{\partial^2}{\partial x^2}$
* $y$ 方向：$L^y(\tau) = \boldsymbol{I} - \tau \boldsymbol{B} \frac{\partial}{\partial y} + \frac{1}{2}\tau^2 \boldsymbol{B}^2 \frac{\partial^2}{\partial y^2}$

**引入 Strang 提出的分裂算子：**
$$\tilde{L}_h(\tau) = L_h^x\left(\frac{\tau}{2}\right) L_h^y(\tau) L_h^x\left(\frac{\tau}{2}\right)$$

**构造差分格式 (6.15)：**
$$\boldsymbol{u}_{jm}^{n+1} = \tilde{L}_h(\tau) \boldsymbol{u}_{jm}^n$$

**【问题】：** 证明差分格式 (6.15) 是以二阶精度逼近方程组 (6.9)。

---

## 二、 详细解答过程

为证明该格式具有二阶精度，需要证明**真实解的时间泰勒展开**与**Strang分裂算子的展开式**在 $O(\tau^2)$ 的量级上完全一致。

### 1. 真实解的泰勒展开 (连续偏微分方程层面)

对真实解 $\boldsymbol{u}(t+\tau)$ 在 $t$ 处进行泰勒展开：
$$\boldsymbol{u}(t+\tau) = \boldsymbol{u}(t) + \tau \boldsymbol{u}_t + \frac{1}{2}\tau^2 \boldsymbol{u}_{tt} + O(\tau^3)$$

利用原始方程组 $\boldsymbol{u}_t = -\boldsymbol{A}\boldsymbol{u}_x - \boldsymbol{B}\boldsymbol{u}_y$，将其对时间求导，把对时间的导数转换为对空间的导数：
$$\boldsymbol{u}_{tt} = \frac{\partial}{\partial t}(-\boldsymbol{A}\boldsymbol{u}_x - \boldsymbol{B}\boldsymbol{u}_y) = -\boldsymbol{A}(\boldsymbol{u}_t)_x - \boldsymbol{B}(\boldsymbol{u}_t)_y$$

将 $\boldsymbol{u}_t$ 再次代入：
$$\boldsymbol{u}_{tt} = -\boldsymbol{A}(-\boldsymbol{A}\boldsymbol{u}_x - \boldsymbol{B}\boldsymbol{u}_y)_x - \boldsymbol{B}(-\boldsymbol{A}\boldsymbol{u}_x - \boldsymbol{B}\boldsymbol{u}_y)_y$$
$$\boldsymbol{u}_{tt} = \boldsymbol{A}^2 \boldsymbol{u}_{xx} + \boldsymbol{A}\boldsymbol{B}\boldsymbol{u}_{yx} + \boldsymbol{B}\boldsymbol{A}\boldsymbol{u}_{xy} + \boldsymbol{B}^2 \boldsymbol{u}_{yy}$$

假设解足够光滑，混合偏导数可交换（即 $\boldsymbol{u}_{xy} = \boldsymbol{u}_{yx}$），则有：
$$\boldsymbol{u}_{tt} = \boldsymbol{A}^2 \boldsymbol{u}_{xx} + (\boldsymbol{A}\boldsymbol{B} + \boldsymbol{B}\boldsymbol{A})\boldsymbol{u}_{xy} + \boldsymbol{B}^2 \boldsymbol{u}_{yy}$$

将其代回泰勒展开式，提取出连续微分算子：
$$\boldsymbol{u}(t+\tau) = \left[ \boldsymbol{I} - \tau\left(\boldsymbol{A}\frac{\partial}{\partial x} + \boldsymbol{B}\frac{\partial}{\partial y}\right) + \frac{1}{2}\tau^2\left(\boldsymbol{A}^2\frac{\partial^2}{\partial x^2} + (\boldsymbol{A}\boldsymbol{B} + \boldsymbol{B}\boldsymbol{A})\frac{\partial^2}{\partial x \partial y} + \boldsymbol{B}^2\frac{\partial^2}{\partial y^2}\right) \right] \boldsymbol{u}(t) + O(\tau^3)$$

这说明，任何一个二阶精度的时间格式，其等效微分算子展开后必须与中括号内的各项完全一致。

### 2. Strang 分裂算子的展开与对比

为了书写简洁，记 $\partial_x = \frac{\partial}{\partial x}$, $\partial_y = \frac{\partial}{\partial y}$。将一维连续微分算子代入 Strang 分裂算子 $\tilde{L}_h(\tau)$，并展开到 $O(\tau^2)$ 级。

首先，时间步长减半的 $x$ 方向算子为：
$$L^x\left(\frac{\tau}{2}\right) = \boldsymbol{I} - \frac{\tau}{2} \boldsymbol{A} \partial_x + \frac{\tau^2}{8} \boldsymbol{A}^2 \partial_x^2 + O(\tau^3)$$

计算后两项的乘积 $L^y(\tau) L^x\left(\frac{\tau}{2}\right)$，忽略高阶项：
$$\left( \boldsymbol{I} - \tau \boldsymbol{B} \partial_y + \frac{\tau^2}{2} \boldsymbol{B}^2 \partial_y^2 \right) \left( \boldsymbol{I} - \frac{\tau}{2} \boldsymbol{A} \partial_x + \frac{\tau^2}{8} \boldsymbol{A}^2 \partial_x^2 \right)$$
$$= \boldsymbol{I} - \frac{\tau}{2} \boldsymbol{A} \partial_x + \frac{\tau^2}{8} \boldsymbol{A}^2 \partial_x^2 - \tau \boldsymbol{B} \partial_y + \frac{\tau^2}{2} \boldsymbol{B}\boldsymbol{A} \partial_y \partial_x + \frac{\tau^2}{2} \boldsymbol{B}^2 \partial_y^2$$

接着，左乘第一个算子 $L^x\left(\frac{\tau}{2}\right)$ 得到完整的 $\tilde{L}(\tau)$：
$$\tilde{L}(\tau) = \left( \boldsymbol{I} - \frac{\tau}{2} \boldsymbol{A} \partial_x + \frac{\tau^2}{8} \boldsymbol{A}^2 \partial_x^2 \right) \left[ \boldsymbol{I} - \frac{\tau}{2} \boldsymbol{A} \partial_x - \tau \boldsymbol{B} \partial_y + \frac{\tau^2}{8} \boldsymbol{A}^2 \partial_x^2 + \frac{\tau^2}{2} \boldsymbol{B}\boldsymbol{A} \partial_y \partial_x + \frac{\tau^2}{2} \boldsymbol{B}^2 \partial_y^2 \right]$$

按 $\tau$ 的幂次对乘积进行合并提取：
* **零阶项 $O(1)$:**
    $$\boldsymbol{I}$$
* **一阶项 $O(\tau)$:**
    $$\boldsymbol{I} \cdot \left(-\frac{\tau}{2} \boldsymbol{A} \partial_x - \tau \boldsymbol{B} \partial_y\right) + \left(-\frac{\tau}{2} \boldsymbol{A} \partial_x\right) \cdot \boldsymbol{I} = -\tau(\boldsymbol{A} \partial_x + \boldsymbol{B} \partial_y)$$
* **二阶项 $O(\tau^2)$:**
    由各项相乘组合得到：
    $$\left(\frac{1}{8} + \frac{1}{4} + \frac{1}{8}\right)\tau^2 \boldsymbol{A}^2 \partial_x^2 + \frac{1}{2}\tau^2 (\boldsymbol{A}\boldsymbol{B} + \boldsymbol{B}\boldsymbol{A}) \partial_x \partial_y + \frac{1}{2}\tau^2 \boldsymbol{B}^2 \partial_y^2$$
    化简后为：
    $$\frac{1}{2}\tau^2 \boldsymbol{A}^2 \partial_x^2 + \frac{1}{2}\tau^2 (\boldsymbol{A}\boldsymbol{B} + \boldsymbol{B}\boldsymbol{A}) \partial_x \partial_y + \frac{1}{2}\tau^2 \boldsymbol{B}^2 \partial_y^2$$

### 3. 结论归纳

将推导出的 Strang 分裂算子 $\tilde{L}(\tau)$ 重组，得到：
$$\tilde{L}(\tau) = \boldsymbol{I} - \tau(\boldsymbol{A} \partial_x + \boldsymbol{B} \partial_y) + \frac{1}{2}\tau^2\left[ \boldsymbol{A}^2 \partial_x^2 + (\boldsymbol{A}\boldsymbol{B} + \boldsymbol{B}\boldsymbol{A}) \partial_x \partial_y + \boldsymbol{B}^2 \partial_y^2 \right] + O(\tau^3)$$

对比第 1 步中真实解的泰勒展开式，可以看出两者在 $O(1)$、$O(\tau)$ 和 $O(\tau^2)$ 的系数上**完全相等**。这严格证明了差分格式 (6.15) 是以二阶精度逼近原方程组 (6.9) 的。

*(注：该方法的巧妙之处在于通过 $x \to y \to x$ 的对称分裂，自然地构造出了交叉偏导数的对称项 $(\boldsymbol{A}\boldsymbol{B} + \boldsymbol{B}\boldsymbol{A})$，从而在矩阵 $\boldsymbol{A}$ 和 $\boldsymbol{B}$ 不可交换的情况下，依然消除了时间分裂带来的误差。)*


# 题目与解答：扩散方程差分格式的稳定性分析

## 一、 题目提出

**试讨论扩散方程 $\frac{\partial u}{\partial t} = \frac{\partial^2 u}{\partial x^2}$ 的差分格式**

$$\frac{u_j^{n+1} - \frac{1}{2}(u_{j+1}^n + u_{j-1}^n)}{\tau} = \frac{1}{h^2}\delta_x^2 u_j^n$$

**的稳定性。**

---

## 二、 详细解答过程

为了分析该差分格式的稳定性，我们采用 **von Neumann（冯·诺依曼）稳定性分析方法**。

### 1. 整理差分格式的显式表达

已知中心二阶差分算子展开为：
$$\delta_x^2 u_j^n = u_{j+1}^n - 2u_j^n + u_{j-1}^n$$

定义网格比（Courant数）为 $r = \frac{\tau}{h^2}$。将算子展开式代入原差分格式中：
$$\frac{u_j^{n+1} - \frac{1}{2}(u_{j+1}^n + u_{j-1}^n)}{\tau} = \frac{1}{h^2}(u_{j+1}^n - 2u_j^n + u_{j-1}^n)$$

两边同乘 $\tau$，并移项整理出 $u_j^{n+1}$ 的显式表达式：
$$u_j^{n+1} - \frac{1}{2}(u_{j+1}^n + u_{j-1}^n) = r(u_{j+1}^n - 2u_j^n + u_{j-1}^n)$$
$$u_j^{n+1} = \left(r + \frac{1}{2}\right)u_{j+1}^n - 2ru_j^n + \left(r + \frac{1}{2}\right)u_{j-1}^n$$

### 2. 求增长因子 (放大因子)

令差分方程的解具有傅里叶模的形式：$u_j^n = V^n e^{ikjh}$，其中 $i = \sqrt{-1}$ 为虚数单位，$k$ 为波数，$h$ 为空间步长，$V^n$ 为第 $n$ 个时间步的振幅。

将该解代入显式格式中：
$$V^{n+1}e^{ikjh} = \left(r + \frac{1}{2}\right)V^ne^{ik(j+1)h} - 2rV^ne^{ikjh} + \left(r + \frac{1}{2}\right)V^ne^{ik(j-1)h}$$

等式两边同时除以 $V^ne^{ikjh}$，求得增长因子 $G = \frac{V^{n+1}}{V^n}$：
$$G = \left(r + \frac{1}{2}\right)e^{ikh} - 2r + \left(r + \frac{1}{2}\right)e^{-ikh}$$
$$G = \left(r + \frac{1}{2}\right)(e^{ikh} + e^{-ikh}) - 2r$$

利用欧拉公式 $e^{ikh} + e^{-ikh} = 2\cos(kh)$ 进行化简：
$$G = 2\left(r + \frac{1}{2}\right)\cos(kh) - 2r$$
$$G = 2r\cos(kh) + \cos(kh) - 2r$$
$$G = \cos(kh) - 2r(1 - \cos(kh))$$

利用半角公式的变形 $1 - \cos(kh) = 2\sin^2\left(\frac{kh}{2}\right)$，以及 $\cos(kh) = 1 - 2\sin^2\left(\frac{kh}{2}\right)$，将其完全转化为正弦平方项：
$$G = 1 - 2\sin^2\left(\frac{kh}{2}\right) - 4r\sin^2\left(\frac{kh}{2}\right)$$
$$G = 1 - 2(1 + 2r)\sin^2\left(\frac{kh}{2}\right)$$

### 3. 稳定性条件讨论

根据 von Neumann 稳定性准则，差分格式稳定的必要条件是，对于任意的相位角 $kh$，增长因子必须满足 $|G| \leq 1$。即：
$$-1 \leq 1 - 2(1 + 2r)\sin^2\left(\frac{kh}{2}\right) \leq 1$$

我们拆解为两个不等式进行讨论：

**（1）右侧不等式：**
$$1 - 2(1 + 2r)\sin^2\left(\frac{kh}{2}\right) \leq 1$$
$$\implies -2(1 + 2r)\sin^2\left(\frac{kh}{2}\right) \leq 0$$
由于物理上网格比 $r = \frac{\tau}{h^2} > 0$，且平方项恒大于等于 0，此不等式恒成立。

**（2）左侧不等式：**
$$-1 \leq 1 - 2(1 + 2r)\sin^2\left(\frac{kh}{2}\right)$$
化简得：
$$(1 + 2r)\sin^2\left(\frac{kh}{2}\right) \leq 1$$
为了保证格式稳定，该不等式必须对所有的 $kh$ 成立。当考虑最极端的高频误差情况（即 $\sin^2\left(\frac{kh}{2}\right) = 1$ 时），必须有：
$$1 + 2r \leq 1$$
$$2r \leq 0$$
$$r \leq 0$$

---

## 三、 结论

在实际的数值计算中，时间步长 $\tau$ 和空间步长 $h$ 必须为正数，因此网格比 $r = \frac{\tau}{h^2}$ 必须严格大于 0。

条件 $r \leq 0$ 在计算中是无法满足的。这意味着无论时间步长取多小，某些高频误差模态总是会满足 $|G| > 1$，从而随着时间步的推进而呈指数级放大。

**因此，该差分格式是绝对不稳定（无条件不稳定）的，不能用于实际求解。**