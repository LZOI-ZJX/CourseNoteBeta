# 第五章 椭圆型方程的差分方法（续二）

## 第四节 变系数方程

考虑变系数椭圆方程:

$$
\frac {\partial}{\partial x} (a (x, y) \frac {\partial u}{\partial x}) + \frac {\partial}{\partial y} (a (x, y) \frac {\partial u}{\partial y}) - c (x, y) u = f (x, y)
$$

其中 $a(x, y)>0,\ c(x,y)\geqslant 0$（原稿中多出未出现的 $b(x,y)>0$，已删）。

### 直接差分方法

$$
\frac {1}{h ^ {2}} \delta_ {x} (a _ {i j} \delta_ {x} u _ {i j}) + \frac {1}{k ^ {2}} \delta_ {y} (a _ {i j} \delta_ {y} u _ {i j}) - c _ {i j} u _ {i j} = f _ {i j}
$$

其中 $h$、$k$ 分别为 $x$、$y$ 方向的网格步长（与下节有限体积离散一致）。

注：为书写简便，式中 $\delta_x(a_{ij}\delta_x u_{ij})$、$\delta_y(a_{ij}\delta_y u_{ij})$ 为简记；严格的半节点（$i\pm\frac12,\ j\pm\frac12$）展开与下文有限体积离散一致。

$$
\begin{array}{l} \iint_ {D _ {i j}} \frac {\partial}{\partial x} (a (x, y) \frac {\partial u}{\partial x}) d x d y \\ = \int_ {y _ {j - \frac {1}{2}}} ^ {y _ {j + \frac {1}{2}}} [ a (x _ {i + \frac {1}{2}}, y) \frac {\partial u}{\partial x} (x _ {i + \frac {1}{2}}, y) - a (x _ {i - \frac {1}{2}}, y) \frac {\partial u}{\partial x} (x _ {i - \frac {1}{2}}, y) ] d y \\ \end{array}
$$

对上面定积分利用中点矩形公式有

$$
\begin{array}{l} \iint_ {D _ {i j}} \frac {\partial}{\partial x} (a (x, y) \frac {\partial u}{\partial x}) d x d y \approx [ a (x _ {i + \frac {1}{2}}, y _ {j}) \frac {\partial u}{\partial x} (x _ {i + \frac {1}{2}}, y _ {j}) \\ - a (x _ {i - \frac {1}{2}}, y _ {j}) \frac {\partial u}{\partial x} (x _ {i - \frac {1}{2}}, y _ {j}) ] k \\ \end{array}
$$

$$
\begin{array}{l} \iint_ {D _ {i j}} \frac {\partial}{\partial y} \left(a (x, y) \frac {\partial u}{\partial y}\right) d x d y \approx \\ \quad \left[ a \left(x _ {i}, y _ {j + \frac {1}{2}}\right) \frac {\partial u}{\partial y} \left(x _ {i}, y _ {j + \frac {1}{2}}\right) - a \left(x _ {i}, y _ {j - \frac {1}{2}}\right) \frac {\partial u}{\partial y} \left(x _ {i}, y _ {j - \frac {1}{2}}\right) \right] h \\ \end{array}
$$

$$
\int_ {D _ {i j}} c (x, y) u d x d y \approx c _ {i j} u _ {i j} h k, \int_ {D _ {i j}} f (x, y) d x d y \approx f _ {i j} h k
$$

于是有

$$
\begin{array}{l} \frac {1}{h ^ {2}} [ a _ {i + \frac {1}{2}, j} (u _ {i + 1, j} - u _ {i j}) - a _ {i - \frac {1}{2}, j} (u _ {i j} - u _ {i - 1, j}) ] \\ + \frac {1}{k ^ {2}} [ a _ {i, j + \frac {1}{2}} (u _ {i, j + 1} - u _ {i j}) - a _ {i , j - \frac {1}{2}} (u _ {i j} - u _ {i, j - 1}) ] - c _ {i j} u _ {i j} = f _ {i j} \\ \end{array}
$$

## 第六节 特征值问题

### 一维模型问题

#### 问题提法（Statement）

Find nontrivial $(u, \lambda)$ such that 

$$
- u _ {x x} = \lambda u, \quad x \in (0, 1)
$$

$$
u (0) = u (1) = 0;
$$

denote solutions $(u^{k},\lambda^{k})$ $k = 1,2,\ldots$ ，with 

$$
0 \leq \lambda^ {1} \leq \lambda^ {2} \leq \dots
$$

#### 精确解（Exact Solution）

$$
- u _ {x x} - \lambda u = 0
$$

$$
u = A \sin \sqrt {\lambda} x + B \cos \sqrt {\lambda} x
$$

$$
u (0) = 0 \Rightarrow B = 0
$$

$$
u (1) = 0 \Rightarrow A = 0 \text{ or } \lambda = k ^ {2} \pi^ {2}, k = 1, 2, \dots
$$

于是（取 $A = 1$） 

$$
\left. \begin{array}{c} u ^ {k} = \sin k \pi x \\ \lambda^ {k} = k ^ {2} \pi^ {2} \end{array} \right\} k = 1, 2 \dots
$$

Larger $k \Rightarrow$ more oscillatory $u^{k} \Rightarrow$ larger $\lambda$ . 

#### 前四阶特征函数示意

**第 1 特征函数**

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/f1511e5124d8cdbc70dfcbe1ccf9a936edc4bf773f081d2318efa9f46b05f9a5.jpg)

**第 2 特征函数**

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/1924dac3714ab98fd1d93fb817dc985346084ddc5717eefa246fc5f2a0a45132.jpg)

**第 3 特征函数**

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/2bb649f88cb5a3a5853bc1d06046f7bdee68f80e27dc599c9692830dce8646b0.jpg)

**第 4 特征函数**

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/4d11d3c553e23b10efdc4138e933e706e1a1e00764070e30945634d41f3727dd.jpg)

### 一维离散方程与差分格式

$$
\begin{array}{l} - u _ {x x} = \lambda u, \quad u (0) = u (1) = 0 \\ \Downarrow \\ \end{array}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/3795aae7500477dcccf46eb9ea34bbbe3b864e88bb7a122a24c37f896e89e2bc.jpg)


$$
\Delta x = \frac {1}{n + 1}
$$

$$
\frac {- 1}{\Delta x ^ {2}} \left(\hat {u} _ {j - 1} - 2 \hat {u} _ {j} + \hat {u} _ {j + 1}\right) = \hat {\lambda} \hat {u} _ {j}, \quad j = 1, \dots , n
$$

$$
\hat {u} _ {0} = \hat {u} _ {n + 1} = 0
$$

### 一维误差分析

#### 离散特征向量与解析解

记离散特征向量为 $\boldsymbol{\hat{u}}^{k} = (\hat{u}_1^k,\ldots,\hat{u}_n^k)^{\mathsf T}$。在网点上与连续本征函数 $u^k(x)=\sin(k\pi x)$ 分量一致，可简记为

$$
\boldsymbol{\hat{u}}^{k} \equiv \boldsymbol{u}^{k} \quad \text{（即 } \hat{u}_j^k = u^k(x_j),\ j=1,\ldots,n\text{）}
$$

$$
\begin{array}{l} \hat {u} _ {j} ^ {k} = u ^ {k} \left(x _ {j}\right) = \sin (k \pi x _ {j}) \\ = \sin (k \pi j \Delta x) = \sin (\frac {k \pi j}{n + 1}), \quad j = 1, \dots , n \\ \end{array}
$$

Note $\hat{u}_0^k = \hat{u}_{n + 1}^k = 0$ since $\sin (0) = \sin (k\pi) = 0.$ 

#### $\hat{\lambda}^k$ 的推导

求 $\hat{\lambda}^k$（下式左端即离散 Laplace 算子作用在 $\hat{u}_j^k$ 上的结果，应等于 $\hat{\lambda}^k \hat{u}_j^k$）：

$$
- \frac {1}{\Delta x ^ {2}} \left\{\hat {u} _ {j - 1} ^ {k} - 2 \hat {u} _ {j} ^ {k} + \hat {u} _ {j + 1} ^ {k} \right\}
$$

$$
\begin{array}{l} = - \frac {1}{\Delta x ^ {2}} \left\{\sin \left(k \pi \left(x _ {j} - \Delta x\right)\right) - 2 \sin \left(k \pi x _ {j}\right) + \sin \left(k \pi \left(x _ {j} + \Delta x\right)\right) \right\} \\ = - \frac {1}{\Delta x ^ {2}} \underbrace {\left. \left. \sin (k \pi x _ {j} - k \pi \Delta x) + \sin (k \pi x _ {j} + k \pi \Delta x) \right. \right.} _ {2 \cos (k \pi \Delta x) \sin (k \pi x _ {j})} - 2 \sin (k \pi x _ {j}) \} \\ \end{array}
$$

#### 特征值恒等式与离散谱

于是

$$
- \frac {1}{\Delta x ^ {2}} \left\{\hat {u} _ {j - 1} ^ {k} - 2 \hat {u} _ {j} ^ {k} + \hat {u} _ {j + 1} ^ {k} \right\}
$$

$$
\begin{array}{l} = - \frac {1}{\Delta x ^ {2}} \left\{2 \cos (k \pi \Delta x) \sin (k \pi x _ {j}) - 2 \sin (k \pi x _ {j}) \right\} \\ = \frac {2}{\Delta x ^ {2}} \{1 - \cos (k \pi \Delta x) \} \hat {u} _ {j} ^ {k}. \\ \end{array}
$$

对比差分方程 $\displaystyle -\frac{1}{\Delta x^2}\bigl(\hat{u}_{j-1}^k - 2\hat{u}_j^k + \hat{u}_{j+1}^k\bigr) = \hat{\lambda}^k \hat{u}_j^k$（即 $-\delta_x^2 \hat{u}_j^k / \Delta x^2 = \hat{\lambda}^k \hat{u}_j^k$），两边约去非零因子 $\hat{u}_j^k$（对内点 $j$ 成立），得离散特征值

$$
\hat {\lambda} ^ {k} = \frac {2}{\Delta x ^ {2}} \{1 - \cos (k \pi \Delta x) \}.
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/767db047d547777a6d8638847233aa0a703de989c7d34519bb0c816f4bc511ae.jpg)


$$
A \boldsymbol{\hat{u}} ^ {k} = \hat {\lambda} ^ {k} \boldsymbol{\hat{u}} ^ {k}
$$

#### 图示补充

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/c0ebd445934ca36d7b9442d461809718c3beef9d1f8f51dddc1b3226e4b6ef08.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/aa4be45d2215558d82bb69d67864b63b988a2a68bcb351d9df525cd2f092e1b2.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/c042306e4f4ade72ec09970751e90d17b81e92b5271ae8fd855bcfb619c58fb3.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/f40c3add0b1e6bad8a7be9baa55c77febac5563ac43f2d4582cced20a6392dd9.jpg)


#### 结论：低模态（Low modes）

For fixed $k$ , $\Delta x \to 0$ : 

$$
\begin{array}{l} \hat {\lambda} ^ {k} = \frac {2}{\Delta x ^ {2}} \{1 - \cos (k \pi \Delta x) \} \\ = \frac {2}{\Delta x ^ {2}} \left\{1 - \left(1 - \frac {1}{2} k ^ {2} \pi^ {2} \Delta x ^ {2} + O \left(\Delta x ^ {4}\right)\right) \right\} \\ = k ^ {2} \pi^ {2} + \mathcal {O} \left(\Delta x ^ {2}\right) \\ \end{array}
$$

$$
\text{second-order convergence},\quad \hat {\lambda} ^ {k} \rightarrow \lambda^ {k}.
$$

#### 结论：高模态（High modes）

For $k = n$ 

$$
\Delta x = \frac {1}{n + 1}
$$

$$
\begin{array}{l} \hat {\lambda} ^ {n} = \frac {2}{\Delta x ^ {2}} \{1 - \cos (\frac {n \pi}{n + 1}) \} \\ = 4 (n + 1) ^ {2} \quad \text{as} \quad \Delta x \rightarrow 0 \\ \neq n ^ {2} \pi^ {2} = \lambda^ {n}. \\ \end{array}
$$

High modes $(k\approx n)$ are not accurate. 

#### 低模态与高模态对比（例）

Low modes vs. high modes 

Example: $n = 19, \Delta x = 1 / 20$ 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/7a638da5254f67211422e7031407d1a81af83baab6059bce7fcc9a336c56dedc.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/1ae3c0a9dd7b208c6ba81f6bbc58a310d44194c28c1143edc8d1b706aaab8e52.jpg)


#### 模态分辨与误差阶

Low modes vs. high modes 

$$
k \ll n
$$

$$
k \approx n
$$

$\boldsymbol{\hat{u}}^{k}$ resolved 

$\hat{\lambda}^{k}$ accurate 

$$
\hat {\lambda} ^ {k} - \lambda^ {k} \sim \mathcal {O} (\Delta x ^ {2})
$$

$\boldsymbol{\hat{u}}^{k}$ not resolved 

$\hat{\lambda}^k$ not accurate 

$\hat{\lambda}^{k}-\lambda^{k}$ is $\mathcal{O}(1)$ 

BUT: as $\Delta x \to 0$ , $n \to \infty$ , so any fixed mode $k$ converges. 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/68e37ae3c0b757e7bde927163653b056e7855a25013c13c21120310bc3b83ec3.jpg)

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/15dfa48f36154ff4c3ad81aee3ac2922846b84f095ccbb7d10949d743cb84cfc.jpg)

### 矩阵 $A$ 的条件数（一维）

For a SPD matrix $M$ , the condition number $\kappa_{M}$ is given by 

$$
\kappa_ {M} = \frac {\text{maximum eigenvalue of } M}{\text{minimum eigenvalue of } M}.
$$

Thus, for our $A$ matrix, 

$$
\kappa_ {A} \rightarrow \frac {4 n ^ {2}}{\pi^ {2}} \quad \text{as} \quad \Delta x \rightarrow 0
$$

grows (in $\mathbb{R}^{1}$) as the number of grid points squared.

Importance: understanding solution procedures. 

### 二维特征值问题

#### 问题提法（Statement）

$$
- \nabla^ {2} u = \lambda u \quad \text{in } \Omega
$$

$$
u = 0 \quad \text{on } \Gamma
$$

Assume (for simplicity) 

$$
L _ {x} = L _ {y} = 1
$$

Solutions $(u(x,y),\lambda)$ 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/c2fcb9094e2d2bf7022e97d20fdd0aa7fb14e47c7521d36ad036585ca8dd31ab.jpg)


#### 连续谱：本征函数与本征值

$$
\begin{array}{l} u ^ {k, l} (x, y) = \sin (k \pi x) \sin (l \pi y) \\ - \nabla^ {2} u ^ {k, l} = \left(k ^ {2} \pi^ {2} + l ^ {2} \pi^ {2}\right) u ^ {k, l} \\ \end{array}
$$

$$
\lambda^ {k, l} = k ^ {2} \pi^ {2} + l ^ {2} \pi^ {2}, \quad k, l = 1, \dots
$$

#### 离散问题与离散特征向量

$$
A \boldsymbol{\hat{u}} ^ {k, l} = \hat {\lambda} ^ {k, l} \boldsymbol{\hat{u}} ^ {k, l} \quad \Rightarrow \quad \left(\boldsymbol{\hat{u}} ^ {k, l}, \hat {\lambda} ^ {k, l}\right)
$$

$$
\begin{array}{l} \hat{u} _ {i, j} ^ {k, l} = \sin (k \pi x _ {i}) \sin (l \pi y _ {j}) \\ = \sin (k \pi i \Delta x) \sin (l \pi j \Delta y) \\ = \sin (\frac {k \pi i}{n + 1}) \sin (\frac {l \pi j}{m + 1}) \\ \end{array}
$$

$$
k, i = 1, \dots , n \quad l, j = 1, \dots , m
$$

#### 离散特征值与高低模态

$$
\hat {\lambda} ^ {k, l} = \frac {2}{\Delta x ^ {2}} \{1 - \cos (k \pi \Delta x) \} + \frac {2}{\Delta y ^ {2}} \{1 - \cos (l \pi \Delta y) \}
$$

Low Modes 

$$
\Delta x, \Delta y \rightarrow 0 \quad (k, l \ \text{fixed})
$$

$$
\begin{array}{l} \hat {\lambda} ^ {k, l} = k ^ {2} \pi^ {2} + l ^ {2} \pi^ {2} \\ + \mathcal {O} (\Delta x ^ {2}, \Delta y ^ {2}) \\ \end{array}
$$

High Modes 

$$
k \approx n, l \approx m
$$

$$
\begin{array}{l} \hat {\lambda} ^ {k, l} = 4 (n + 1) ^ {2} + 4 (m + 1) ^ {2} \\ \text{as} \quad \Delta x, \Delta y \rightarrow 0 \\ \end{array}
$$

#### 矩阵 $A$ 的条件数（二维）

$$
\kappa_ {A} \rightarrow \frac {4 n ^ {2} + 4 m ^ {2}}{2 \pi^ {2}} \quad \text{as} \quad \Delta x, \Delta y \rightarrow 0
$$

$$
\text{If } m \approx n
$$

$$
\kappa_ {A} \rightarrow \frac {4 n ^ {2}}{\pi^ {2}}
$$

grows (in $\mathbb{R}^2$) with the number of grid points in each coordinate direction.

**注：** 若记一维未知量总数为 $N=n$，则 $\kappa_A=O(n^2)=O(N^2)$。若二维 $m\approx n$，未知量总数 $N\approx nm\approx n^2$，且 $\kappa_A=O(n^2)=O(N)$。相对**总自由度** $N$ 而言，二维条件数随 $N$ 的增长为 $O(N)$，慢于一维的 $O(N^2)$，故常说相对一维情形“更温和”（二者比较的是**相同意义下的问题规模 $N$**，而非单边的分段子数 $n$）。