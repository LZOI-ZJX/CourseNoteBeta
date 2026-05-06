# 第四章 抛物型方程的有限差分法

## 主要内容

- 常系数扩散方程  
- 初边值问题  
- 对流扩散方程  
- 变系数方程  
- 多维问题  

## 第一节 常系数扩散方程

考虑一维扩散方程：

$$
\frac {\partial u}{\partial t} = a \frac {\partial^ {2} u}{\partial x ^ {2}}, 0 \leq t \leq T, (1. 1)
$$

其中 $a$ 是正常数。

可将(1.1)的定解问题分为两类：

第一、初值问题（也称Cauchy问题）：求具有所需次数偏微商的函数 $u(x,t)$ ，满足方程(1.1)和初始条件：

$$
\left\{ \begin{array}{c} \frac {\partial u}{\partial t} = a \frac {\partial^ {2} u}{\partial x ^ {2}}, 0 <   t \leq T \\ u (x, 0) = \varphi (x), - \infty <   x <   + \infty \end{array} \right. \tag {1.2}
$$

第二、初边值问题（也称混合问题）：求具有所需次数偏微商的函数 $u(x,t)$ ，满足方程(1.1)、初始条件和边值条件：

$$
\left\{ \begin{array}{c} \frac {\partial u}{\partial t} = a \frac {\partial^ {2} u}{\partial x ^ {2}}, 0 <   t \leq T \\ u (x, 0) = \varphi (x), 0 <   x <   l (1. 3) _ {1} \\ u (0, t) = u (l, t) = 0, 0 \leq t \leq T (1. 3) _ {2} \end{array} \right.
$$

假定 $\varphi(x)$ 在相应区域光滑，并且在 $x = 0, l$ 满足相容条件，使上述问题有唯一充分光滑的解。

### 1、划分网格

现在考虑边值问题(1.1)，(1.3)的差分逼近。取空间

步长 $h = \frac{l}{N}$ 和时间步长 $\tau = \frac{T}{M}$ , 其中 $N, M$ 都是自然数。用

两族平行直线 $x = x_{j} = j h (j = 0,1,\dots ,N)$ 和 $t = t_{n} = n \tau (n =$ 

$0, 1, \dots, M)$ 将矩形域 $\overline{G} = \left\{0 \leq x \leq l ; 0 \leq t \leq T \right\}$ 分割成

矩形网格,网格节点为 $(x_{j}, t_{n})$。

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/5a47ec83-83d1-451d-9b1a-5309ec4fa7fb/04f203e2812b29060f18a33870a5edb5707ca11d7b94a97d869fbcaaa6f977d1.jpg)


### 2、离散微分方程建立差分方程

#### (1) 向前差分格式与向后差分格式

##### 向前差分格式

$$
\begin{array}{l} \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} = a \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}} (1. 4) _ {1} \\ u _ {j} ^ {0} = \varphi_ {j} = \varphi (x _ {j}), u _ {0} ^ {n} = u _ {N} ^ {n} = 0, (1. 4) _ {2} \\ \end{array}
$$

其中 $j = 1,2,\dots ,N - 1,n = 1,2,\dots ,M - 1.$ 

计算的公式:

以 $\lambda = \frac{\tau}{h^{2}}$ 表示网比。

$$
u _ {j} ^ {n + 1} = a \lambda u _ {j + 1} ^ {n} + (1 - 2 a \lambda) u _ {j} ^ {n} + a \lambda u _ {j - 1} ^ {n} \quad (1. 4) _ {1} ^ {\prime}
$$

记 $Lu = \frac{\partial u}{\partial t} - a \frac{\partial^2 u}{\partial x^2}$ ,

$$
L _ {h} ^ {(1)} u _ {j} ^ {n} = \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} - a \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}
$$

用Taylor展开分析截断误差

$$
\begin{array}{l} R (u (x _ {j}, t _ {n})) = L _ {h} ^ {(1)} u (x _ {j}, t _ {n}) - L u (x _ {j}, t _ {n}) \\ = a \tau \left( \frac {a}{2} - \frac {1}{1 2 \lambda} \right) \left( \frac {\partial^ {4} \tilde {u}}{\partial x ^ {4}} \right) _ {j} ^ {n} + O (\tau^ {2} + h ^ {2}) \\ = O (\tau + h ^ {2}) \\ \end{array}
$$

（首项在形式上写成 $u_{xxxx}$：时间离散贡献 $\frac{\tau}{2} u_{tt}$，空间 Laplacian 离散贡献 $-\frac{a h^2}{12} u_{xxxx}$；又用 $u_{tt}=a^2 u_{xxxx}$ 合并，$\lambda=\tau/h^2$。）


##### 向前差分格式的性质

截断误差 $O\left(\tau + h^{2}\right)$ 

利用Fourier方法分析稳定性：

增长因子为: $G(\lambda, k) = 1 - 4a\lambda \sin^2{\frac{k h}{2}}$ 

稳定性条件为: $a \lambda \leq \frac{1}{2}$ 

##### 向后差分格式

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} = a \frac {u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}}{h ^ {2}} \tag {1.6}
$$

$$
u _ {j} ^ {0} = \varphi_ {j} = \varphi (x _ {j}), u _ {0} ^ {n} = u _ {N} ^ {n} = 0, (1. 6) _ {2}
$$

其中 $j = 1,2,\dots ,N - 1,n = 1,2,\dots ,M - 1$ 

将 $(1.6)_{1}$ 改写成便于计算的形式,

$$
- a \lambda u _ {j + 1} ^ {n + 1} + (1 - 2 a \lambda) u _ {j} ^ {n + 1} - a \lambda u _ {j - 1} ^ {n + 1} = u _ {j} ^ {n} (1. 6) _ {1} ^ {\prime}
$$

记 $L_{h}^{(2)}u_{j}^{n} = \frac{u_{j}^{n + 1} - u_{j}^{n}}{\tau} -a\frac{u_{j + 1}^{n + 1} - 2u_{j}^{n + 1} + u_{j - 1}^{n + 1}}{h^{2}}$ 

显然截断误差

$$
\begin{array}{l} R (u (x _ {j}, t _ {n})) = L _ {h} ^ {(2)} u (x _ {j}, t _ {n}) - L u (x _ {j}, t _ {n}) \\ = - a \tau \left( \frac {a}{2} + \frac {1}{1 2 \lambda} \right) \left( \frac {\partial^ {4} \tilde {u}}{\partial x ^ {4}} \right) _ {j} ^ {n} + O (\tau^ {2} + h ^ {2}) \\ = O \left(\tau + h ^ {2}\right) \tag {1.7} \\ \end{array}
$$

（向后差分对时间的展开与向前差分相差一个符号，故出现 $+\frac{a}{2}$；合并方式同前。）


##### 向后差分格式的性质

截断误差 $O\left(\tau + h^{2}\right)$ 

利用Fourier方法分析稳定性：

增长因子为: $G(\lambda, k) = \frac{1}{1 + 4a\lambda\sin^2{\frac{k h}{2}}}$ 

稳定性为：无条件稳定

#### (2) 加权隐式格式

$$
\begin{array}{l} \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} = a \left[ (1 - \theta) \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}} \right. \\ \left. + \theta \frac {u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}}{h ^ {2}} \right] \\ \end{array}
$$

$$
\boldsymbol {u} _ {j} ^ {0} = \varphi_ {j} = \varphi (\boldsymbol {x} _ {j}), \boldsymbol {u} _ {0} ^ {n} = \boldsymbol {u} _ {N} ^ {n} = 0,
$$

其中 $j = 1,2,\dots ,N - 1,n = 1,2,\dots ,M - 1.$ 此处 **$\theta$ 为隐式层（$t_{n+1}$）权重**，$1-\theta$ 为显式层权重，与下文 Von Neumann 增长因子 $G$ 中分母、分子的写法一致。

$$
0 \leq \theta \leq 1
$$

计算公式:

$$
\begin{array}{l} - a \lambda \theta u _ {j + 1} ^ {n + 1} + (1 + 2 a \lambda \theta) u _ {j} ^ {n + 1} - a \lambda \theta u _ {j - 1} ^ {n + 1} = \\ a \lambda (1 - \theta) u _ {j + 1} ^ {n} + (1 - 2 a \lambda (1 - \theta)) u _ {j} ^ {n} + a \lambda (1 - \theta) u _ {j - 1} ^ {n} \\ \end{array}
$$

截断误差: $E = a\left(\theta - \frac{1}{2}\right) \tau \left[\frac{\partial^{3} u}{\partial x^{2}\,\partial t}\right]_{j}^{n} + O\left(\tau^{2} + h^{2}\right)$ 

$$
E = \left\{ \begin{array}{l l} O (\tau + h ^ {2}) & \theta \neq \frac {1}{2} \\ O (\tau^ {2} + h ^ {2}) & \theta = \frac {1}{2} \end{array} \right.
$$

##### Fourier 方法与过渡因子（记号：以下 $G(\lambda,k)$ 中的 $\lambda=\tau/h^2$ 恒为**网比**，勿与矩阵特征根混淆）

$$
G (\lambda , k) = \frac {1 - 4 (1 - \theta) a \lambda \sin^ {2} \frac {k h}{2}}{1 + 4 \theta a \lambda \sin^ {2} \frac {k h}{2}}
$$

稳定性条件为：

$$
2 a \lambda \leq \frac {1}{1 - 2 \theta}, \quad \text {当} 0 \leq \theta <   \frac {1}{2}
$$

$$
\text {无 条 件 稳 定 ， 当} \frac {1}{2} \leq \theta \leq 1
$$

$$
\text {Crank-Nicolson 格式 } \left(\theta = \frac {1}{2}\right)\text{：}
$$

$$
\begin{array}{l} \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} = \frac {a}{2} [ \frac {u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}}{h ^ {2}} + \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}} ] \\ u _ {j} ^ {0} = \varphi_ {j} = \varphi (x _ {j}), u _ {0} ^ {n} = u _ {N} ^ {n} = 0 \\ \end{array}
$$

计算公式:

$$
\begin{array}{l} - \frac {a \lambda}{2} u _ {j + 1} ^ {n + 1} + (1 + a \lambda) u _ {j} ^ {n + 1} - \frac {a \lambda}{2} u _ {j - 1} ^ {n + 1} \\ = \frac {a \lambda}{2} u _ {j + 1} ^ {n} + (1 - a \lambda) u _ {j} ^ {n} + \frac {a \lambda}{2} u _ {j - 1} ^ {n} \\ \end{array}
$$

##### Crank-Nicolson 格式的性质

Crank-Nicolson格式的截断误差:

$$
O (\tau^ {2} + h ^ {2})
$$

Crank-Nicolson格式的稳定性:

无条件稳定

#### (3) 三层显式格式

##### Richardson 格式

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n - 1}}{2 \tau} = a \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}}
$$

计算公式:

$$
u _ {j} ^ {n + 1} = 2 a \lambda (u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}) + u _ {j} ^ {n - 1}
$$

注：Richardson格式是显格式，截断误差的阶为 $O(\tau^{2} + h^{2})$ ，但从稳定性方面来看，它是不可用的。

用 $e_{j}^{n}$ 表示 $\pmb{u}_{j}^{n}$ 的误差

则 $e_{j}^{n}$ 满足与 $(1.10)'$ 相应的齐次方程:

$$
e _ {j} ^ {n + 1} = 2 a \lambda \left(e _ {j + 1} ^ {n} - 2 e _ {j} ^ {n} + e _ {j - 1} ^ {n}\right) + e _ {j} ^ {n - 1}. (1. 1 1)
$$

设误差只在初始层的原点 $(j = 0)$ 发生, 即 $e_{j}^{0} = \delta_{j0} \varepsilon$ 

$(\varepsilon > 0; \delta_{00} = 1, \delta_{j0} = 0$ 当 $j \neq 0), e_j^{-1} = 0,$ 而在以后计算

中都是精确的，则初始误差的传递情况如表1：


表1 $a \lambda = \frac{1}{2}$ 时 Richardson 格式的误差传播


<table><tr><td>\</td><td>-4</td><td>-3</td><td>-2</td><td>-1</td><td>0</td><td>1</td><td>2</td><td>3</td><td>4</td></tr><tr><td>0</td><td>0</td><td>0</td><td>0</td><td>0</td><td>ε</td><td>0</td><td>0</td><td>0</td><td>0</td></tr><tr><td>1</td><td>0</td><td>0</td><td>0</td><td>ε</td><td>-2ε</td><td>ε</td><td>0</td><td>0</td><td>0</td></tr><tr><td>2</td><td>0</td><td>0</td><td>ε</td><td>-4ε</td><td>7ε</td><td>-4ε</td><td>ε</td><td>0</td><td>0</td></tr><tr><td>3</td><td>0</td><td>ε</td><td>-6ε</td><td>17ε</td><td>-24ε</td><td>17ε</td><td>-6ε</td><td>ε</td><td>0</td></tr><tr><td>4</td><td>ε</td><td>-8ε</td><td>31ε</td><td>-68ε</td><td>89ε</td><td>-68ε</td><td>31ε</td><td>-8ε</td><td>ε</td></tr><tr><td>5</td><td>-10ε</td><td>49ε</td><td>-144ε</td><td>273ε</td><td>-388ε</td><td>273ε</td><td>-144ε</td><td>49ε</td><td>-10ε</td></tr><tr><td>6</td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td><td></td></tr></table>

##### Fourier 方法分析 Richardson 格式的稳定性

一般将上述三层格式化为等价的二层方程组（截断误差仍为 $O(\tau^{2}+h^{2})$），引入 $v_{j}^{n+1}=u_{j}^{n}$，写为

$$
\left\{ \begin{array}{l} u _ {j} ^ {n + 1} = u _ {j} ^ {n - 1} + 2 a \lambda (u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}) \\ v _ {j} ^ {n + 1} = u _ {j} ^ {n} \end{array} \right.
$$
（与上文计算公式一致，网比 $\lambda=\tau/h^2$。）

令 $U = (u, v)^T$ , 则有:

$$
U _ {j} ^ {n + 1} = \left( \begin{array}{c c} 2 a \lambda & 0 \\ 0 & 0 \end{array} \right) U _ {j + 1} ^ {n} + \left( \begin{array}{c c} - 4 a \lambda & 1 \\ 1 & 0 \end{array} \right) U _ {j} ^ {n} + \left( \begin{array}{c c} 2 a \lambda & 0 \\ 0 & 0 \end{array} \right) U _ {j - 1} ^ {n}
$$

令 $U_{j}^{n} = V^{n} e^{ikjh}$ 有:

$$
V ^ {n + 1} e ^ {i k j h}
$$

$$
\begin{array}{l} = \left\{\left( \begin{array}{c c} 2 a \lambda & 0 \\ 0 & 0 \end{array} \right) V ^ {n} e ^ {- i k h} + \left( \begin{array}{c c} - 4 a \lambda & 1 \\ 1 & 0 \end{array} \right) V ^ {n} + \left( \begin{array}{c c} 2 a \lambda & 0 \\ 0 & 0 \end{array} \right) V ^ {n} e ^ {i k h} \right\} e ^ {i k j h} \\ = \left( \begin{array}{c c} - 4 a \lambda + 2 a \lambda (e ^ {i k h} + e ^ {- i k h}) & 1 \\ 1 & 0 \end{array} \right) V ^ {n} e ^ {i k j h} \\ \end{array}
$$

传播矩阵 $G(\lambda) = \left( \begin{array}{ccc} - 8a\lambda \sin^2{\frac{k h}{2}} & 1\\ 1 & 0 \end{array} \right)$ 

其特征值 $\mu$ 满足方程: $\mu^{2} + 8a\lambda \mu \sin^{2}\frac{k h}{2} -1 = 0$ 

$$
\mu_ {1, 2} = - 4 a \lambda \sin^ {2} \frac {k h}{2} \pm (1 + 1 6 a ^ {2} \lambda^ {2} \sin^ {4} \frac {k h}{2}) ^ {\frac {1}{2}}
$$

显然对于任何网格比 $\lambda$ , 有

$$
\rho (G) > 1 + 4 a \lambda \sin^ {2} \frac {k h}{2} > 1
$$

破坏了 Von Neumann 条件，格式不稳定。此格式是无条件不稳定的。

##### Du Fort-Frankel 格式

1953年 Du Fort-Frankel 格式

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n - 1}}{2 \tau} - a \frac {u _ {j + 1} ^ {n} - (u _ {j} ^ {n + 1} + u _ {j} ^ {n - 1}) + u _ {j - 1} ^ {n}}{h ^ {2}} = 0
$$

它的节点结构图为

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/5a47ec83-83d1-451d-9b1a-5309ec4fa7fb/c393e905d3edb8306c0eb820f263222c2bce752926e600f56b48d6443eccbf21.jpg)


这是一个三层的格式，其截断误差为：

$$
E = a \left(\frac {\tau}{h}\right) ^ {2} \left[ \frac {\partial^ {2} u}{\partial t ^ {2}} \right] _ {j} ^ {n} + o \left(\tau^ {2} + h ^ {2}\right) + o \left(\frac {\tau^ {4}}{h ^ {2}}\right)
$$

由此可见，若要使差分格式与微分方程相容，要求：当 $\tau \rightarrow 0$ 时，有 $\frac{\tau}{h} \rightarrow 0$ ，即差分方程与微分方程相容的充分条件是： $\tau \rightarrow 0$ 时， $\tau$ 是 $h$ 的高阶无穷小（这当 $\lambda = \frac{\tau}{h^2} \rightarrow \text{const}$ 时是能满足的）。

反之, 如果 $\frac{\tau}{h} \rightarrow \beta (\tau \rightarrow 0)$ , 则差分格式不与扩散方程 $\frac{\partial u}{\partial t} - a \frac{\partial^2 u}{\partial x^2} = 0$ 相容, 而是下面双曲线方程相容

$$
\frac {\partial u}{\partial t} - a \frac {\partial^ {2} u}{\partial x ^ {2}} + a \beta^ {2} \frac {\partial^ {2} u}{\partial t ^ {2}} = 0
$$

下面研究 Du Fort-Frankel 格式的稳定性：先将格式化成等价的矩阵形式：

$$
\left\{ \begin{array}{l} (1 + 2 a \lambda) u _ {j} ^ {n + 1} = (1 - 2 a \lambda) v _ {j} ^ {n} + 2 a \lambda (u _ {j + 1} ^ {n} + u _ {j - 1} ^ {n}) \\ v _ {j} ^ {n + 1} = u _ {j} ^ {n} \end{array} \right.
$$

$$
U _ {j} ^ {n + 1} = \left(u _ {j} ^ {n + 1}, v _ {j} ^ {n + 1}\right) ^ {T} U _ {j} ^ {n} = V ^ {n} e ^ {i k j h}
$$

$$
\begin{array}{l} \left( \begin{array}{c c} 1 + 2 a \lambda & 0 \\ 0 & 1 \end{array} \right) V ^ {n + 1} e ^ {i k j h} = \left\{\left( \begin{array}{c c} 2 a \lambda & 0 \\ 0 & 0 \end{array} \right) V ^ {n} e ^ {- i k h} + \left( \begin{array}{c c} 0 & 1 - 2 a \lambda \\ 1 & 0 \end{array} \right) V ^ {n} + \left( \begin{array}{c c} 2 a \lambda & 0 \\ 0 & 0 \end{array} \right) V ^ {n} e ^ {i k h} \right\} e ^ {i k j h} \\ = \left( \begin{array}{c c} 2 a \lambda (e ^ {i k h} + e ^ {- i k h}) & 1 - 2 a \lambda \\ 1 & 0 \end{array} \right) V ^ {n} e ^ {i k j h} \\ \end{array}
$$

$$
\left( \begin{array}{c c} 1 + 2 a \lambda & 0 \\ 0 & 1 \end{array} \right) G (\lambda , k) = \left( \begin{array}{c c} 0 & 1 - 2 a \lambda \\ 1 & 0 \end{array} \right) + \left( \begin{array}{c c} 4 a \lambda \cos k h & 0 \\ 0 & 0 \end{array} \right)
$$

$G(\lambda, k)$ 的特征方程是 $(\alpha = 2a\lambda)$ 

$$
\mu^ {2} - (\frac {2 \alpha}{1 + \alpha} \cos k h) \mu - \frac {1 - \alpha}{1 + \alpha} = 0
$$

可以证明：增长矩阵 $G(\lambda,k)$（此处 $\lambda=\tau/h^2$ 仍表**网比**）的特征根 $\mu_{1,2}$ 满足 $|\mu_{1,2}| \leq 1$。当两根互异且为**共轭复根**时，由 $\mu^2 - b\mu - c = 0$（本格式中 $c=\frac{1-\alpha}{1+\alpha}$）知 $\mu_1\mu_2=-c$，故 $|\mu_{1,2}| = \sqrt{|\mu_1\mu_2|} = \sqrt{\left|\frac{1-\alpha}{1+\alpha}\right|}$（**不能**误用为 $\left|\frac{\alpha\cos kh}{1+\alpha}\right|$，后者仅为实根情形根的实部相关量）；若 $\alpha>0$ 仍有 $\left|\frac{1-\alpha}{1+\alpha}\right|<1$，故格式无条件稳定。分析基于如下的引理:

引理1.1

引理：实系数二次方程 $\mu^2 - b\mu - c = 0$ 的根满足 $|\mu_{1,2}| \leq 1$ 的充分必要条件是

$$
| b | \leq 1 - c \quad {\text {且}} \quad | c | \leq 1
$$

虽然DFF格式是无条件稳定的，但其相容性是有条件的，即 $\tau / h \rightarrow 0 (\tau \rightarrow 0)$ ，故对于具体的 $\tau, h$ , 不见得一定相容。事实上，无法构造出无条件相容和无条件稳定的显式格式。

#### (4) 三层隐式格式

##### 三层隐式格式 I

$$
\boxed {\frac {3}{2} \frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} - \frac {1}{2} \frac {u _ {j} ^ {n} - u _ {j} ^ {n - 1}}{\tau} - a \frac {u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}}{h ^ {2}} = 0}
$$

或令 $\lambda = \tau / h^{2}$ , 将其改写成:

$$
(3 + 4 a \lambda) u _ {j} ^ {n + 1} - 2 a \lambda (u _ {j + 1} ^ {n + 1} + u _ {j - 1} ^ {n + 1}) = 4 u _ {j} ^ {n} - u _ {j} ^ {n - 1}
$$

此格式的截断误差为: $E = o\left(\tau^{2} + h^{2}\right)$ . 节点图:

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/5a47ec83-83d1-451d-9b1a-5309ec4fa7fb/ca0e7aa24db7ead0eea092a7612cbffc85be575078b4cdc6f5f4be80e3c03a04.jpg)


引入变量 $\boldsymbol{\nu}_{j}^{n+1} = \boldsymbol{u}_{j}^{n}$ , 将格式写成矩阵形式, 从而得到增长矩阵：

$$
G (\lambda , k) = \left( \begin{array}{c c} \frac {4}{3 + 8 \alpha \sin^ {2} \frac {k h}{2}} & \frac {- 1}{3 + 8 \alpha \sin^ {2} \frac {k h}{2}} \\ 1 & 0 \end{array} \right) \qquad \alpha = a \lambda
$$

特征方程：

$$
\mu^ {2} - \frac {4}{3 + 8 \alpha \sin^ {2} \frac {k h}{2}} \cdot \mu + \frac {1}{3 + 8 \alpha \sin^ {2} \frac {k h}{2}} = 0
$$

由于 $|b| < 1 - |c|$ , $|c| < 1$ , 由引理1.1知 $\left|\mu_{1,2}\right| \leq 1$ 

且 $\mu_{1,2} = \frac{2 \pm \sqrt{1 - 8\alpha\sin^2{\frac{k h}{2}}}}{3 + 8\alpha\sin^2{\frac{k h}{2}}}$ 

若重根, $\mu_{1,2} = \frac{2}{3 + 8\alpha \sin^2{\frac{k h}{2}}}< 1$ 

故格式I是无条件稳定的。

##### 三层隐式格式 II

设 $\delta_{x}^{2} u_{j}^{n} = u_{j+1}^{n} - 2u_{j}^{n} + u_{j-1}^{n}$ , 格式:

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n - 1}}{2 \tau} - \frac {a}{3 h ^ {2}} \left(\delta_ {x} ^ {2} u _ {j} ^ {n + 1} + \delta_ {x} ^ {2} u _ {j} ^ {n} + \delta_ {x} ^ {2} u _ {j} ^ {n - 1}\right) = 0
$$

这是用 $\frac{1}{3} \left(\delta_{x}^{2} u_{j}^{n+1} + \delta_{x}^{2} u_{j}^{n} + \delta_{x}^{2} u_{j}^{n-1}\right)$ 代替了Richardson格式中的 $\delta_{x}^{2} u_{j}^{n}$ 而已，也可理解为将二层的Crank-Nicolson格式向三层格式的推广。这也是一个二阶格式: $E = O\left(\tau^{2} + h^{2}\right)$ 节点图：

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/5a47ec83-83d1-451d-9b1a-5309ec4fa7fb/d5c76efc6cfecb85e15c3ae352378d7f6b1d95ac2b5d904cf99b79b08513b106.jpg)


将格式改写成：

$$
(1 - \frac {2}{3} a \lambda \delta_ {x} ^ {2}) u _ {j} ^ {n + 1} = \frac {2}{3} a \lambda \delta_ {x} ^ {2} u _ {j} ^ {n} + (1 + \frac {2}{3} a \lambda \delta_ {x} ^ {2}) u _ {j} ^ {n - 1}
$$

引入变量 $v_{j}^{n + 1} = u_{j}^{n}$ , 差分方程可写成矩阵形式

$$
\left\{ \begin{array}{l} (1 - \frac {2}{3} a \lambda \delta_ {x} ^ {2}) u _ {j} ^ {n + 1} = \frac {2}{3} a \lambda \delta_ {x} ^ {2} u _ {j} ^ {n} + (1 + \frac {2}{3} a \lambda \delta_ {x} ^ {2}) v _ {j} ^ {n} \\ v _ {j} ^ {n + 1} = u _ {j} ^ {n} \end{array} \right.
$$

$$
\left( \begin{array}{c c} 1 - \frac {2}{3} a \lambda \delta_ {x} ^ {2} & 0 \\ 0 & 1 \end{array} \right) U _ {j} ^ {n + 1} = \left( \begin{array}{c c} \frac {2}{3} a \lambda \delta_ {x} ^ {2} & (1 + \frac {2}{3} a \lambda \delta_ {x} ^ {2}) \\ 1 & 0 \end{array} \right) U _ {j} ^ {n}
$$

令 $U_{j}^{n + 1} = V^{n + 1}e^{ikjh}$ 

并注意到: $\delta_{x}^{2} u_{j}^{n} = -4 \sin^{2} \frac{k h}{2} v^{n}$ 

可得增长矩阵（式中 $\lambda=\tau/h^2$ 为网比）：

$$
G (k, \lambda) = \left( \begin{array}{c c} 1 + \frac {8}{3} a \lambda \sin^ {2} \frac {k h}{2} & 0 \\ 0 & 1 \end{array} \right) ^ {- 1} \left( \begin{array}{c c} - \frac {8}{3} a \lambda \sin^ {2} \frac {k h}{2} & 1 - \frac {8}{3} a \lambda \sin^ {2} \frac {k h}{2} \\ 1 & 0 \end{array} \right)
$$

$$
= \left( \begin{array}{c c} \frac {- \alpha}{1 + \alpha} & \frac {1 - \alpha}{1 + \alpha} \\ 1 & 0 \end{array} \right)
$$

$$
\alpha = \frac {8}{3} a \lambda \sin^ {2} \frac {k h}{2}
$$

$G$ 的特征方程（关于特征根 $\mu$）为：

$$
\mu^ {2} + \frac {\alpha}{1 + \alpha} \mu - \frac {1 - \alpha}{1 + \alpha} = 0
$$

由引理1.1 有 $|\mu_{1,2}| \leq 1$ , 且 $\mu_{1,2} = \frac{-\alpha \pm \sqrt{4 - 3\alpha^2}}{2(1 + \alpha)}$ 

当有重根时, 有 $|\mu| < 1$ , 故格式II是无条件稳定的.

#### (5) 跳点格式

所谓跳点格式, 是将网格点 $\left(x_{j}, t_{n}\right)$ 按 $n + j =$ 奇偶数分成两组, 按奇偶网格点的不同, 用不同的格式计算。

从 $t_{n}$ 推进到 $t_{n+1}$ 时, 先在偶数网格上用古典显格式

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} - a \frac {u _ {j + 1} ^ {n} - 2 u _ {j} ^ {n} + u _ {j - 1} ^ {n}}{h ^ {2}} = 0, n + 1 + j = \text {偶 数}
$$

然后在奇数格点上，用隐式格式：

$$
\frac {u _ {j} ^ {n + 1} - u _ {j} ^ {n}}{\tau} - a \frac {u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}}{h ^ {2}} = 0 \quad n + j + 1 = \text {奇 数}
$$

计算奇数格点上的值 $\boldsymbol{u}_{j}^{n + 1}$ , 当 $n + j + 1 =$ 奇数时, $\boldsymbol{u}_{j \pm 1}^{n + 1}$ 是偶数格点上的值, 隐式格式可整理成:

$$
u _ {j} ^ {n + 1} = \frac {u _ {j} ^ {n} + a \lambda \left(u _ {j + 1} ^ {n + 1} + u _ {j - 1} ^ {n + 1}\right)} {1 + 2 a \lambda}
$$

$n + j + 2$ 是偶数，采用如下格式

$$
\frac {u _ {j} ^ {n + 2} - u _ {j} ^ {n + 1}}{\tau} - a \frac {u _ {j + 1} ^ {n + 1} - 2 u _ {j} ^ {n + 1} + u _ {j - 1} ^ {n + 1}}{h ^ {2}} = 0
$$

比较 $n + 1$ 层的隐格式可以得到 $u_{j}^{n + 2} = 2 u_{j}^{n + 1} - u_{j}^{n}$ 

代入奇点的隐格式 $u_{j}^{n + 1} = \dfrac{u_{j}^{n} + a\lambda\left(u_{j + 1}^{n + 1} + u_{j - 1}^{n + 1}\right)}{1 + 2a\lambda}$ 

$$
\frac {u _ {j} ^ {n + 2} - u _ {j} ^ {n}}{2 \tau} - a \frac {u _ {j + 1} ^ {n + 1} - (u _ {j} ^ {n + 2} + u _ {j} ^ {n}) + u _ {j - 1} ^ {n + 1}}{h ^ {2}} = 0
$$

$n + j + 2 =$ 偶数

故奇数格点上的隐式格式实质上是一个显式格式。

可以说明跳点格式与 Du Fort-Frankel 格式等价，故精度和稳定性相同，但跳点格式计算方便，占用内存少。