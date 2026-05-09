## 题目5.4

在区域 $D=\{(x,y)|0 \leqslant x,y \leqslant 1\}$ 上给出边值问题：

$$
\begin{cases}
-\left(\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}\right) = 16, & 0 < x,y < 1, \\
u \big|_{x=1} = 0, \\
\left.\frac{\partial u}{\partial y}\right|_{y=1} = -u, \\
\left.\frac{\partial u}{\partial x}\right|_{x=0} = 0, \\
\left.\frac{\partial u}{\partial y}\right|_{y=0} = 0.
\end{cases}
$$

取 $h=\frac{1}{4}$，试用五点差分格式求此问题的数值解。

---

**差分解法与格式推导**

为了保持五点差分格式 $O(h^2)$ 的截断误差精度，对于含有导数的边界条件（Neumann 边界和 Robin 边界），采用**中心差分法**并引入**虚节点（Ghost points）**进行离散。

已知步长 $h = \frac{1}{4}$，定义网格节点 $(x_i, y_j) = (ih, jh)$，其中 $i, j = 0, 1, 2, 3, 4$。记 $u_{i,j}$ 为精确解 $u(x_i, y_j)$ 的近似值。

内部节点的标准五点差分格式为：
$$4u_{i,j} - u_{i-1,j} - u_{i+1,j} - u_{i,j-1} - u_{i,j+1} = 16h^2$$
代入 $h = \frac{1}{4}$，右端项 $16h^2 = 1$。

**边界条件离散**

1. **右边界 $x=1$ （Dirichlet 边界）**
   边界上函数值已知，直接赋值：
   $$u_{4,j} = 0 \quad (j = 0, 1, 2, 3, 4)$$

2. **左边界 $x=0$ （Neumann 边界）**
   条件为 $\left.\frac{\partial u}{\partial x}\right|_{x=0} = 0$。引入 $x$ 方向的虚节点 $x_{-1} = -h$：
   $$\frac{u_{1,j} - u_{-1,j}}{2h} = 0 \implies u_{-1,j} = u_{1,j}$$
   代入 $i=0$ 处的内部五点差分方程，消去虚节点得：
   $$4u_{0,j} - 2u_{1,j} - u_{0,j-1} - u_{0,j+1} = 1 \quad (j = 1, 2, 3)$$

3. **下边界 $y=0$ （Neumann 边界）**
   条件为 $\left.\frac{\partial u}{\partial y}\right|_{y=0} = 0$。引入 $y$ 方向的虚节点 $y_{-1} = -h$：
   $$\frac{u_{i,1} - u_{i,-1}}{2h} = 0 \implies u_{i,-1} = u_{i,1}$$
   代入 $j=0$ 处的方程，消去虚节点得：
   $$4u_{i,0} - u_{i-1,0} - u_{i+1,0} - 2u_{i,1} = 1 \quad (i = 1, 2, 3)$$

4. **上边界 $y=1$ （Robin 边界）**
   条件为 $\left.\frac{\partial u}{\partial y}\right|_{y=1} = -u$。引入 $y$ 方向的虚节点 $y_5 = 1+h$：
   $$\frac{u_{i,5} - u_{i,3}}{2h} = -u_{i,4} \implies u_{i,5} = u_{i,3} - 2hu_{i,4}$$
   代入 $h = \frac{1}{4}$，即 $u_{i,5} = u_{i,3} - \frac{1}{2}u_{i,4}$。将该虚节点代入 $j=4$ 处的方程中：
   $$4u_{i,4} - u_{i-1,4} - u_{i+1,4} - u_{i,3} - \left(u_{i,3} - \frac{1}{2}u_{i,4}\right) = 16h^2$$
   整理得到：
   $$4.5u_{i,4} - u_{i-1,4} - u_{i+1,4} - 2u_{i,3} = 1 \quad (i = 1, 2, 3)$$

*(注：在两导数边界相交的角点处，需同时应用两个方向的虚节点进行消元。)*

---

**矩阵形式与求解**

将未知数按 $y$ 优先、$x$ 其次的顺序排列为列向量 $U = [U_0^T, U_1^T, U_2^T, U_3^T, U_4^T]^T$，其中 $U_j = [u_{0,j}, u_{1,j}, u_{2,j}, u_{3,j}]^T$。此时原问题转化为求解大型稀疏线性方程组 $AU = b$：

$$
\begin{bmatrix}
T & -2I & 0 & 0 & 0 \\
-I & T & -I & 0 & 0 \\
0 & -I & T & -I & 0 \\
0 & 0 & -I & T & -I \\
0 & 0 & 0 & -2I & T_4
\end{bmatrix}
\begin{bmatrix}
U_0 \\ U_1 \\ U_2 \\ U_3 \\ U_4
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{1} \\ \mathbf{1} \\ \mathbf{1} \\ \mathbf{1} \\ \mathbf{1}
\end{bmatrix}
$$

其中 $I$ 是 $4 \times 4$ 单位矩阵，$\mathbf{1} = [1, 1, 1, 1]^T$。块矩阵定义为：

$$T = \begin{bmatrix} 4 & -2 & 0 & 0 \\ -1 & 4 & -1 & 0 \\ 0 & -1 & 4 & -1 \\ 0 & 0 & -1 & 4 \end{bmatrix}, \quad T_4 = \begin{bmatrix} 4.5 & -2 & 0 & 0 \\ -1 & 4.5 & -1 & 0 \\ 0 & -1 & 4.5 & -1 \\ 0 & 0 & -1 & 4.5 \end{bmatrix}$$

**MATLAB 求解代码**
```matlab
%| eval: false
%========================================================
% 五点差分法求解泊松方程边值问题
% 区域：0 <= x, y <= 1, 步长 h = 1/4
%========================================================
clear; clc;

% 未知数总数: i=0,1,2,3 (4个), j=0,1,2,3,4 (5个), 共计20个
N = 20; 

% 初始化系数矩阵 A 和右端项 b
A = zeros(N, N);
b = ones(N, 1); % 右端项 16 * h^2 = 16 * (1/4)^2 = 1

% 定义 4x4 的基础块矩阵
I = eye(4);

T = [ 4, -2,  0,  0;
     -1,  4, -1,  0;
      0, -1,  4, -1;
      0,  0, -1,  4];

T4 = [4.5, -2,  0,  0;
     -1, 4.5, -1,  0;
      0, -1, 4.5, -1;
      0,  0, -1, 4.5];

% 组装 20x20 的大型分块矩阵 A
% j = 0
A(1:4, 1:4) = T;       A(1:4, 5:8) = -2*I;
% j = 1
A(5:8, 1:4) = -I;      A(5:8, 5:8) = T;       A(5:8, 9:12) = -I;
% j = 2
A(9:12, 5:8) = -I;     A(9:12, 9:12) = T;     A(9:12, 13:16) = -I;
% j = 3
A(13:16, 9:12) = -I;   A(13:16, 13:16) = T;   A(13:16, 17:20) = -I;
% j = 4
A(17:20, 13:16) = -2*I; A(17:20, 17:20) = T4;

% 求解线性方程组 A * U = b
U = A \ b;

% 将一维解向量重构成 5x5 的二维网格矩阵以便观察
u_grid = zeros(5, 5); % 初始化含边界的全网格

for j = 0:4
    for i = 0:3
        idx = j * 4 + i + 1;
        u_grid(j+1, i+1) = U(idx);
    end
end
% 右边界 u=0 已由 zeros 初始化，无需修改。

% 翻转矩阵：使得输出符合直角坐标系视觉习惯 (y=1 在上方)
u_display = flipud(u_grid);

% 打印结果
disp('【网格点数值解 (按笛卡尔坐标系排列)】:');
disp('            x=0       x=1/4      x=1/2      x=3/4      x=1');
disp('        ---------------------------------------------------');
y_labels = {'y=1  |', 'y=3/4|', 'y=1/2|', 'y=1/4|', 'y=0  |'};
for row = 1:5
    fprintf('%s %10.4f %10.4f %10.4f %10.4f %10.4f\n', y_labels{row}, u_display(row, :));
end
```

## 题目5.5

用柱坐标表示的 Poisson 方程的形式为：

$$
\frac{1}{r} \frac{\partial}{\partial r}\left(r \frac{\partial u}{\partial r}\right) + \frac{1}{r^2} \frac{\partial^2 u}{\partial \varphi^2} + \frac{\partial^2 u}{\partial z^2} = -f(r, \varphi, z),
$$

其中 $r = \sqrt{x^2+y^2}, \tan\varphi = \frac{y}{x}$。试写出逼近上述方程的一个差分格式。

---

**一、 通用网格划分**

设柱坐标系下三个方向的步长分别为 $\Delta r$、$\Delta \varphi$ 和 $\Delta z$。取空间网格节点 $(r_i, \varphi_j, z_k)$，记 $u_{i,j,k}$ 为精确解在该点的数值近似，同理右端项记为 $f_{i,j,k} = f(r_i, \varphi_j, z_k)$。

**二、 解法一：基于泰勒展开的有限差分法（FDM）**

**1. 方程展开**
利用乘积求导法则，将第一项关于 $r$ 的导数项展开：
$$
 \frac{1}{r} \frac{\partial}{\partial r}\left(r \frac{\partial u}{\partial r}\right) = \frac{1}{r} \left( \frac{\partial u}{\partial r} + r \frac{\partial^2 u}{\partial r^2} \right) = \frac{\partial^2 u}{\partial r^2} + \frac{1}{r} \frac{\partial u}{\partial r} 
$$

将其代回原方程，得到等价形式：
$$ 
\frac{\partial^2 u}{\partial r^2} + \frac{1}{r} \frac{\partial u}{\partial r} + \frac{1}{r^2} \frac{\partial^2 u}{\partial \varphi^2} + \frac{\partial^2 u}{\partial z^2} = -f(r, \varphi, z) 
$$

**2. 偏导数中心差分离散**

在内部节点 $(r_i, \varphi_j, z_k)$ 处，采用二阶精度的中心差分商来近似各个偏导数项：

* **$r$ 方向二阶导数：** $\frac{\partial^2 u}{\partial r^2} \approx \frac{u_{i+1,j,k} - 2u_{i,j,k} + u_{i-1,j,k}}{(\Delta r)^2}$
* **$r$ 方向一阶导数：** $\frac{\partial u}{\partial r} \approx \frac{u_{i+1,j,k} - u_{i-1,j,k}}{2\Delta r}$
* **$\varphi$ 方向二阶导数：** $\frac{\partial^2 u}{\partial \varphi^2} \approx \frac{u_{i,j+1,k} - 2u_{i,j,k} + u_{i,j-1,k}}{(\Delta \varphi)^2}$
* **$z$ 方向二阶导数：** $\frac{\partial^2 u}{\partial z^2} \approx \frac{u_{i,j,k+1} - 2u_{i,j,k} + u_{i,j,k-1}}{(\Delta z)^2}$

**3. 最终差分格式组合**

将上述离散式代入展开后的方程中，并提取合并 $u_{i+1,j,k}$ 和 $u_{i-1,j,k}$ 的公因子，得到标准二阶中心差分格式：

$$
 \left( \frac{1}{(\Delta r)^2} + \frac{1}{2 r_i \Delta r} \right) u_{i+1,j,k} + \left( \frac{1}{(\Delta r)^2} - \frac{1}{2 r_i \Delta r} \right) u_{i-1,j,k} + \frac{u_{i,j+1,k} - 2u_{i,j,k} + u_{i,j-1,k}}{r_i^2 (\Delta \varphi)^2} + \frac{u_{i,j,k+1} - 2u_{i,j,k} + u_{i,j,k-1}}{(\Delta z)^2} - \frac{2}{(\Delta r)^2} u_{i,j,k} = -f_{i,j,k} 
 $$

**三、 解法二：基于积分守恒的有限体积法（FVM）**

**1. 定义控制体积**

以节点 $(r_i, \varphi_j, z_k)$ 为中心，建立控制体积 $V_{cv}$，边界定义为：$r \in [r_{i-1/2}, r_{i+1/2}]$，$\varphi \in [\varphi_{j-1/2}, \varphi_{j+1/2}]$，$z \in [z_{k-1/2}, z_{k+1/2}]$。其中半整数角标代表节点间的中点（如 $r_{i+1/2} = r_i + \frac{\Delta r}{2}$）。该控制体积的体积近似为 $V_{cv} = r_i \Delta r \Delta\varphi \Delta z$。

**2. 方程积分**

将原方程两边乘以微元 $r \, dr \, d\varphi \, dz$，并在 $V_{cv}$ 上进行三重积分：

$$ 
\iiint_{V_{cv}} \left[ \frac{\partial}{\partial r}\left(r \frac{\partial u}{\partial r}\right) + \frac{1}{r} \frac{\partial^2 u}{\partial \varphi^2} + r \frac{\partial^2 u}{\partial z^2} \right] dr \, d\varphi \, dz = \iiint_{V_{cv}} -f(r, \varphi, z) \, r \, dr \, d\varphi \, dz 
$$

**3. 近似各界面通量与源项**

利用散度定理和中心差分，将积分转化为界面上的净流入通量：

* **$r$ 方向界面通量：** $\Delta\varphi \Delta z \left[ r_{i+1/2} \frac{u_{i+1,j,k} - u_{i,j,k}}{\Delta r} - r_{i-1/2} \frac{u_{i,j,k} - u_{i-1,j,k}}{\Delta r} \right]$
* **$\varphi$ 方向界面通量：** 取 $\int \frac{dr}{r} \approx \frac{\Delta r}{r_i}$，近似为 $\frac{\Delta r \Delta z}{r_i \Delta\varphi} (u_{i,j+1,k} - 2u_{i,j,k} + u_{i,j-1,k})$
* **$z$ 方向界面通量：** $\frac{r_i \Delta r \Delta\varphi}{\Delta z} (u_{i,j,k+1} - 2u_{i,j,k} + u_{i,j,k-1})$
* **源项积分：** $\approx -f_{i,j,k} (r_i \Delta r \Delta\varphi \Delta z)$

**4. 归一化得出守恒格式**

将上述三项通量相加，等于源项积分。为了转化为标准的离散方程形式，将方程两边同时除以控制体积 $V_{cv} = r_i \Delta r \Delta\varphi \Delta z$，即可得到守恒型的有限体积差分格式：

$$
\frac{1}{r_i \Delta r} \left[ r_{i+1/2} \frac{u_{i+1,j,k} - u_{i,j,k}}{\Delta r} - r_{i-1/2} \frac{u_{i,j,k} - u_{i-1,j,k}}{\Delta r} \right] + \frac{u_{i,j+1,k} - 2u_{i,j,k} + u_{i,j-1,k}}{r_i^2 (\Delta \varphi)^2} + \frac{u_{i,j,k+1} - 2u_{i,j,k} + u_{i,j,k-1}}{(\Delta z)^2} = -f_{i,j,k} 
$$

---

有限体积法的核心思想是**积分守恒**。我们将微分方程在控制体积（Control Volume）上进行积分，利用高斯散度定理将体积积分转化为通过控制体积各个界面的通量积分，从而自然得出守恒型的差分格式。

**1. 定义控制体积 (Control Volume)**

在柱坐标 $(r, \varphi, z)$ 中，体积微元为 $dV = r \, dr \, d\varphi \, dz$。
以空间网格节点 $P(r_i, \varphi_j, z_k)$ 为中心，建立一个控制体积 $V_{cv}$，其界面边界划分如下：

*   **$r$ 方向**：$[r_{i-1/2}, r_{i+1/2}]$，其中 $r_{i\pm1/2} = r_i \pm \frac{\Delta r}{2}$
*   **$\varphi$ 方向**：$[\varphi_{j-1/2}, \varphi_{j+1/2}]$，其中 $\varphi_{j\pm1/2} = \varphi_j \pm \frac{\Delta \varphi}{2}$
*   **$z$ 方向**：$[z_{k-1/2}, z_{k+1/2}]$，其中 $z_{k\pm1/2} = z_k \pm \frac{\Delta z}{2}$

该控制体积的精确体积积分值为：
$$V_{cv} = \int_{z_{k-1/2}}^{z_{k+1/2}} \int_{\varphi_{j-1/2}}^{\varphi_{j+1/2}} \int_{r_{i-1/2}}^{r_{i+1/2}} r \, dr \, d\varphi \, dz = \frac{1}{2}(r_{i+1/2}^2 - r_{i-1/2}^2) \Delta\varphi \Delta z = r_i \Delta r \Delta\varphi \Delta z$$

**2. 对原方程在控制体积上积分**

将原方程两边同时乘以体积微元 $r \, dr \, d\varphi \, dz$，并在控制体积 $V_{cv}$ 上作三重积分：

$$
\iiint_{V_{cv}} \left[ \frac{\partial}{\partial r}\left(r \frac{\partial u}{\partial r}\right) + \frac{1}{r} \frac{\partial^2 u}{\partial \varphi^2} + r \frac{\partial^2 u}{\partial z^2} \right] dr \, d\varphi \, dz = \iiint_{V_{cv}} -f(r, \varphi, z) \, r \, dr \, d\varphi \, dz
$$

**3. 各项积分与界面通量近似计算**

我们将等式左边分为三项，分别计算通过界面的通量：

**第一项（$r$ 方向界面通量）：**

$$
\int_{z}\int_{\varphi}\int_{r} \frac{\partial}{\partial r}\left(r \frac{\partial u}{\partial r}\right) dr \, d\varphi \, dz = \Delta\varphi \Delta z \left[ \left( r \frac{\partial u}{\partial r} \right)_{r_{i+1/2}} - \left( r \frac{\partial u}{\partial r} \right)_{r_{i-1/2}} \right]
$$

利用相邻节点的中心差分近似界面上的导数：

$$
\approx \Delta\varphi \Delta z \left[ r_{i+1/2} \frac{u_{i+1,j,k} - u_{i,j,k}}{\Delta r} - r_{i-1/2} \frac{u_{i,j,k} - u_{i-1,j,k}}{\Delta r} \right]
$$

**第二项（$\varphi$ 方向界面通量）：**

$$
\int_{z}\int_{\varphi}\int_{r} \frac{1}{r} \frac{\partial}{\partial \varphi} \left(\frac{\partial u}{\partial \varphi}\right) dr \, d\varphi \, dz = \Delta z \left(\int_{r_{i-1/2}}^{r_{i+1/2}} \frac{dr}{r} \right) \left[ \left( \frac{\partial u}{\partial \varphi} \right)_{\varphi_{j+1/2}} - \left( \frac{\partial u}{\partial \varphi} \right)_{\varphi_{j-1/2}} \right]
$$

取节点中心 $r_i$ 近似积分 $\int \frac{dr}{r} \approx \frac{\Delta r}{r_i}$，并使用中心差分近似界面导数：

$$
\approx \Delta z \frac{\Delta r}{r_i} \left[ \frac{u_{i,j+1,k} - u_{i,j,k}}{\Delta\varphi} - \frac{u_{i,j,k} - u_{i,j-1,k}}{\Delta\varphi} \right] = \frac{\Delta r \Delta z}{r_i \Delta\varphi} (u_{i,j+1,k} - 2u_{i,j,k} + u_{i,j-1,k})
$$

**第三项（$z$ 方向界面通量）：**
$$
\int_{z}\int_{\varphi}\int_{r} r \frac{\partial}{\partial z} \left(\frac{\partial u}{\partial z}\right) dr \, d\varphi \, dz = \Delta\varphi \left(\int_{r_{i-1/2}}^{r_{i+1/2}} r dr \right) \left[ \left( \frac{\partial u}{\partial z} \right)_{z_{k+1/2}} - \left( \frac{\partial u}{\partial z} \right)_{z_{k-1/2}} \right]
$$

由于 $\int r dr = r_i \Delta r$，代入中心差分：

$$
\approx r_i \Delta r \Delta\varphi \left[ \frac{u_{i,j,k+1} - u_{i,j,k}}{\Delta z} - \frac{u_{i,j,k} - u_{i,j,k-1}}{\Delta z} \right] = \frac{r_i \Delta r \Delta\varphi}{\Delta z} (u_{i,j,k+1} - 2u_{i,j,k} + u_{i,j,k-1})
$$

**源项积分（等式右侧）：**

假定源项 $f$ 在整个控制体积内近似等于中心点的值 $f_{i,j,k}$：
$$
\iiint_{V_{cv}} -f \, r \, dr \, d\varphi \, dz \approx -f_{i,j,k} (r_i \Delta r \Delta\varphi \Delta z)
$$

**4. 组装与归一化**

将上述近似项相加，并将方程两边同时除以控制体积 $V_{cv} = r_i \Delta r \Delta\varphi \Delta z$，得到各个方向的离散形式：

1.  **$r$ 方向项**：
    $$ \frac{1}{r_i \Delta r} \left[ r_{i+1/2} \frac{u_{i+1,j,k} - u_{i,j,k}}{\Delta r} - r_{i-1/2} \frac{u_{i,j,k} - u_{i-1,j,k}}{\Delta r} \right] $$
2.  **$\varphi$ 方向项**：
    $$ \frac{1}{r_i^2 (\Delta\varphi)^2} (u_{i,j+1,k} - 2u_{i,j,k} + u_{i,j-1,k}) $$
3.  **$z$ 方向项**：
    $$ \frac{1}{(\Delta z)^2} (u_{i,j,k+1} - 2u_{i,j,k} + u_{i,j,k-1}) $$

## 题目 6.1（1）

试推导下列泛函在集合 $K$ 中达到极值的必要条件，列出对应的 Euler 方程或方程组：

$$
J(y_1, \cdots, y_n) = \int_a^b F(x; y_1, \cdots, y_n; y_1', \cdots, y_n') \text{d}x,
$$

$$
K = \{(y_1, \cdots, y_n) | y_i \in C^2[a,b], y_i(a)=y_{ia}, y_i(b)=y_{ib}, i=1,\cdots,n\}.
$$

---

**一、 引入变分函数**

假设向量函数 $\mathbf{y}(x) = (y_1(x), y_2(x), \cdots, y_n(x))$ 是使泛函 $J$ 在集合 $K$ 中达到极值的解。

引入 $n$ 个任意的、在 $[a, b]$ 上具有二阶连续导数的变分函数 $\eta_i(x)$（$i=1, 2, \cdots, n$），且满足边界条件：
$$ \eta_i(a) = 0, \quad \eta_i(b) = 0 \quad (i = 1, 2, \cdots, n) $$

构造一族依赖于参数 $\epsilon$ 的函数：
$$ \tilde{y}_i(x) = y_i(x) + \epsilon \eta_i(x) \quad (i = 1, 2, \cdots, n) $$

显然，当 $\epsilon = 0$ 时，$\tilde{y}_i(x) = y_i(x)$。由于 $\eta_i$ 在端点处为零，$\tilde{y}_i(x)$ 满足集合 $K$ 中的边界条件 $\tilde{y}_i(a) = y_{ia}, \tilde{y}_i(b) = y_{ib}$，因此构造的函数族始终有 $(\tilde{y}_1, \cdots, \tilde{y}_n) \in K$。

**二、 利用极值必要条件**

将构造的函数族代入原泛函，泛函 $J$ 转化为关于参数 $\epsilon$ 的单变量函数：
$$ \Phi(\epsilon) = J(y_1 + \epsilon \eta_1, \cdots, y_n + \epsilon \eta_n) = \int_a^b F(x; y_1 + \epsilon \eta_1, \cdots, y_n + \epsilon \eta_n; y_1' + \epsilon \eta_1', \cdots, y_n' + \epsilon \eta_n') \text{d}x $$

由于函数族在 $\epsilon = 0$ 时达到极值，极值的必要条件是泛函的变分为零，即 $\Phi(\epsilon)$ 在 $\epsilon=0$ 处的一阶导数必须为零：
$$ \left. \frac{\text{d}\Phi}{\text{d}\epsilon} \right|_{\epsilon=0} = 0 $$

对 $\Phi(\epsilon)$ 积分号内利用多元复合函数求导法则，并代入 $\epsilon = 0$，得到泛函的一阶变分等式：
$$ \delta J = \int_a^b \sum_{i=1}^n \left( \frac{\partial F}{\partial y_i} \eta_i(x) + \frac{\partial F}{\partial y_i'} \eta_i'(x) \right) \text{d}x = 0 $$

**三、 分部积分化简**

为了提取公因子 $\eta_i(x)$，对上式中的第二部分 $\frac{\partial F}{\partial y_i'} \eta_i'(x)$ 利用分部积分法：
$$ \int_a^b \frac{\partial F}{\partial y_i'} \eta_i'(x) \text{d}x = \left[ \frac{\partial F}{\partial y_i'} \eta_i(x) \right]_a^b - \int_a^b \frac{\text{d}}{\text{d}x} \left( \frac{\partial F}{\partial y_i'} \right) \eta_i(x) \text{d}x $$

根据设定的边界条件 $\eta_i(a) = \eta_i(b) = 0$，上述公式的边界代入项计算结果为 $0$，即 $\left[ \frac{\partial F}{\partial y_i'} \eta_i(x) \right]_a^b = 0$。

将分部积分的结果代回变分等式中，提取公因子 $\eta_i(x)$ 并合并同类项，得到：
$$ \int_a^b \sum_{i=1}^n \left[ \frac{\partial F}{\partial y_i} - \frac{\text{d}}{\text{d}x} \left( \frac{\partial F}{\partial y_i'} \right) \right] \eta_i(x) \text{d}x = 0 $$

**四、 得出 Euler 方程组**

上述积分等式必须对所有满足边界条件的、**任意且相互独立**的变分函数 $\eta_1(x), \cdots, \eta_n(x)$ 都成立。根据变分法基本引理（拉格朗日引理），被积函数中每一个 $\eta_i(x)$ 前面的系数在积分区间 $[a,b]$ 内必须恒等于零。

由此，我们得到了使得多未知函数泛函达到极值的必要条件，即对应的 **Euler 方程组**：
$$ \frac{\partial F}{\partial y_i} - \frac{\text{d}}{\text{d}x} \left( \frac{\partial F}{\partial y_i'} \right) = 0 \quad (i = 1, 2, \cdots, n) $$

## 题目 6.1(2)

试推导下列泛函在集合 $K$ 中达到极值的必要条件：

$$
J(y) = \int_a^b F(x; y, y', y'') \text{d}x,
$$

$$
K = \{y | y \in C^4[a,b], y(a)=y_a, y(b)=y_b, y'(a)=y'_a, y'(b)=y'_b\}.
$$

---

**一、 引入变分函数**

假设函数 $y(x) \in K$ 是使泛函 $J(y)$ 达到极值的解。

引入一个任意的变分函数 $\eta(x)$，要求它在区间 $[a,b]$ 上具有足够高阶的连续导数，并且满足端点条件：
$$ \eta(a) = 0, \quad \eta(b) = 0 $$
$$ \eta'(a) = 0, \quad \eta'(b) = 0 $$

构造依赖于参数 $\epsilon$ 的函数族：
$$ \tilde{y}(x) = y(x) + \epsilon \eta(x) $$

由于 $\eta$ 及其一阶导数在边界处为零，所以对于任意的 $\epsilon$，函数 $\tilde{y}(x)$ 都严格满足集合 $K$ 规定的边界条件，即 $\tilde{y} \in K$。

**二、 利用极值的必要条件**

将构造的函数代入泛函中，得到一个关于 $\epsilon$ 的单变量函数：
$$ \Phi(\epsilon) = J(y + \epsilon \eta) = \int_a^b F(x; y + \epsilon \eta, y' + \epsilon \eta', y'' + \epsilon \eta'') \text{d}x $$

由于 $y(x)$ 是极值函数（对应 $\epsilon = 0$ 的情况），所以 $\Phi(\epsilon)$ 在 $\epsilon = 0$ 处取得极值。其必要条件是 $\Phi(\epsilon)$ 在 $\epsilon = 0$ 处的一阶导数（即一阶变分）为零：
$$ \delta J = \left. \frac{\text{d}\Phi}{\text{d}\epsilon} \right|_{\epsilon=0} = 0 $$

利用积分号下求导法则和多元复合函数求导的链式法则，对上式求导并代入 $\epsilon = 0$，得到：
$$ \int_a^b \left( \frac{\partial F}{\partial y} \eta(x) + \frac{\partial F}{\partial y'} \eta'(x) + \frac{\partial F}{\partial y''} \eta''(x) \right) \text{d}x = 0 $$

**三、 分部积分化简**

为了提取出公因子 $\eta(x)$，我们需要对含有 $\eta'(x)$ 和 $\eta''(x)$ 的项进行分部积分。

首先处理 $\eta'(x)$ 项，进行一次分部积分：
$$ \int_a^b \frac{\partial F}{\partial y'} \eta'(x) \text{d}x = \left[ \frac{\partial F}{\partial y'} \eta(x) \right]_a^b - \int_a^b \frac{\text{d}}{\text{d}x} \left( \frac{\partial F}{\partial y'} \right) \eta(x) \text{d}x $$
因为 $\eta(a) = \eta(b) = 0$，第一项边界值代入后为 $0$，所以：
$$ \int_a^b \frac{\partial F}{\partial y'} \eta'(x) \text{d}x = - \int_a^b \frac{\text{d}}{\text{d}x} \left( \frac{\partial F}{\partial y'} \right) \eta(x) \text{d}x $$

接着处理 $\eta''(x)$ 项，需要进行**两次**分部积分：
第一次分部积分：
$$ \int_a^b \frac{\partial F}{\partial y''} \eta''(x) \text{d}x = \left[ \frac{\partial F}{\partial y''} \eta'(x) \right]_a^b - \int_a^b \frac{\text{d}}{\text{d}x} \left( \frac{\partial F}{\partial y''} \right) \eta'(x) \text{d}x $$
因为 $\eta'(a) = \eta'(b) = 0$，边界项为 $0$。对剩下的积分再次进行分部积分：
$$ - \int_a^b \frac{\text{d}}{\text{d}x} \left( \frac{\partial F}{\partial y''} \right) \eta'(x) \text{d}x = - \left[ \frac{\text{d}}{\text{d}x} \left( \frac{\partial F}{\partial y''} \right) \eta(x) \right]_a^b + \int_a^b \frac{\text{d}^2}{\text{d}x^2} \left( \frac{\partial F}{\partial y''} \right) \eta(x) \text{d}x $$
同样由于 $\eta(a) = \eta(b) = 0$，边界项依然为 $0$。因此：
$$ \int_a^b \frac{\partial F}{\partial y''} \eta''(x) \text{d}x = \int_a^b \frac{\text{d}^2}{\text{d}x^2} \left( \frac{\partial F}{\partial y''} \right) \eta(x) \text{d}x $$

**四、 得出 Euler-Poisson 方程**

将上述两次分部积分的结果代回第一变分的表达式中，合并同类项，得到：
$$ \int_a^b \left[ \frac{\partial F}{\partial y} - \frac{\text{d}}{\text{d}x} \left( \frac{\partial F}{\partial y'} \right) + \frac{\text{d}^2}{\text{d}x^2} \left( \frac{\partial F}{\partial y''} \right) \right] \eta(x) \text{d}x = 0 $$

由于变分函数 $\eta(x)$ 在区间 $(a,b)$ 内是完全任意的，根据变分法基本引理（拉格朗日引理），被积函数中 $\eta(x)$ 前面的系数必须在整个区间内恒为零。

由此得出该泛函达到极值的必要条件，即高阶导数情形下的 **Euler-Poisson 方程**：
$$ \frac{\partial F}{\partial y} - \frac{\text{d}}{\text{d}x} \left( \frac{\partial F}{\partial y'} \right) + \frac{\text{d}^2}{\text{d}x^2} \left( \frac{\partial F}{\partial y''} \right) = 0 $$
