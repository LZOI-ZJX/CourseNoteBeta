# 第四节 变系数方程

考虑变系数椭圆方程:

$$
\frac {\partial}{\partial x} (a (x, y) \frac {\partial u}{\partial x}) + \frac {\partial}{\partial y} (a (x, y) \frac {\partial u}{\partial y}) - c (x, y) u = f (x, y)
$$

其中 $a(x, y)>0,\ c(x,y)\geqslant 0$（原稿中多出未出现的 $b(x,y)>0$，已删）。

1: 直接差分方法

$$
\frac {1}{h ^ {2}} \delta_ {x} (a _ {i j} \delta_ {x} u _ {i j}) + \frac {1}{h ^ {2}} \delta_ {y} (a _ {i j} \delta_ {y} u _ {i j}) - c _ {i j} u _ {i j} = f _ {i j}
$$

# 2：有限体积法（积分差分方法）

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

# 第六节 特征值问题

# The Eigenvalue Problem

# Model Problem

# Statement

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




# Exact Solution

$$
- u _ {x x} - \lambda u = 0
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/13b6ad3ef885c132a0e28b07e9a5a7a62ad14f75ac182b415320c57154448827.jpg)


$$
u = A \sin \sqrt {\lambda} x + B \cos \sqrt {\lambda} x
$$

$$
u (0) = 0 \Rightarrow B = 0
$$

$$
u (1) = 0 \Rightarrow A = 0 \text {o r} \lambda = k ^ {2} \pi^ {2}, k = 1, 2, \dots
$$

# The Eigenvalue Problem

# Exact Solution

Thus (choose $A = 1$ ) 

$$
\left. \begin{array}{c} u ^ {k} = \sin k \pi x \\ \lambda^ {k} = k ^ {2} \pi^ {2} \end{array} \right\} k = 1, 2 \dots
$$

Larger $k \Rightarrow$ more oscillatory $u^{k} \Rightarrow$ larger $\lambda$ . 

# The Eigenvalue Problem

# Exact Solution


1st Eigenfunction


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/f1511e5124d8cdbc70dfcbe1ccf9a936edc4bf773f081d2318efa9f46b05f9a5.jpg)



2nd Eigenfunction


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/1924dac3714ab98fd1d93fb817dc985346084ddc5717eefa246fc5f2a0a45132.jpg)



3rd Eigenfunction


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/2bb649f88cb5a3a5853bc1d06046f7bdee68f80e27dc599c9692830dce8646b0.jpg)



4th Eigenfunction


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/4d11d3c553e23b10efdc4138e933e706e1a1e00764070e30945634d41f3727dd.jpg)


# The Eigenvalue Problem

# Discrete Equations

# Difference Formulas

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

# The Eigenvalue Problem

# Error Analysis

Analytical Solution: $\hat{\boldsymbol{u}}^k,\hat{\lambda}^k\dots$ 

Claim that 

$$
\underline {{\hat {u}}} ^ {k} \equiv \underline {{u}} ^ {k}
$$

$$
\begin{array}{l} \hat {u} _ {j} ^ {k} = u ^ {k} \left(x _ {j}\right) = \sin (k \pi x _ {j}) \\ = \sin (k \pi j \Delta x) = \sin (\frac {k \pi j}{n + 1}), \quad j = 1, \dots , n \\ \end{array}
$$

Note $\hat{u}_0^k = \hat{u}_{n + 1}^k = 0$ since $\sin (0) = \sin (k\pi) = 0.$ 

# The Eigenvalue Problem

# Error Analysis

...Analytical Solution: $\hat{\boldsymbol{u}}^k,\hat{\lambda}^k\dots$ 

What are $\hat{\lambda}^k$ ? 

$$
- \frac {1}{\Delta x ^ {2}} \left\{\hat {u} _ {j - 1} ^ {k} - 2 \hat {u} _ {j} ^ {k} + \hat {u} _ {j + 1} ^ {k} \right\}
$$

$$
\begin{array}{l} = - \frac {1}{\Delta x ^ {2}} \left\{\sin \left(k \pi \left(x _ {j} - \Delta x\right)\right) - 2 \sin \left(k \pi x _ {j}\right) + \sin \left(k \pi \left(x _ {j} + \Delta x\right)\right) \right\} \\ = - \frac {1}{\Delta x ^ {2}} \underbrace {\left. \left. \sin (k \pi x _ {j} - k \pi \Delta x) + \sin (k \pi x _ {j} + k \pi \Delta x) \right. \right.} _ {2 \cos (k \pi \Delta x) \sin (k \pi x _ {j})} - 2 \sin (k \pi x _ {j}) \} \\ \end{array}
$$

# The Eigenvalue Problem

# Error Analysis

...Analytical Solution: $\hat{\underline{u}}^k, \hat{\lambda}^k$ 

Thus: 

$$
- \frac {1}{\Delta x ^ {2}} \left\{\hat {u} _ {j - 1} ^ {k} - 2 \hat {u} _ {j} ^ {k} + \hat {u} _ {j + 1} ^ {k} \right\}
$$

$$
\begin{array}{l} = - \frac {1}{\Delta x ^ {2}} \left\{2 \cos (k \pi \Delta x) \sin (k \pi x _ {j}) - 2 \sin (k \pi x _ {j}) \right\} \\ = \frac {2}{\Delta x ^ {2}} \{1 - \cos (k \pi \Delta x) \} \sin (k \pi x _ {j}). \\ \end{array}
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/767db047d547777a6d8638847233aa0a703de989c7d34519bb0c816f4bc511ae.jpg)


$$
A \hat {u} ^ {k} = \hat {\lambda} ^ {k} \hat {u} ^ {k}
$$

# The Eigenvalue Problem

# Error Analysis

...Analytical Solution: $\hat{\boldsymbol{u}}^k,\hat{\lambda}^k\dots$ 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/c0ebd445934ca36d7b9442d461809718c3beef9d1f8f51dddc1b3226e4b6ef08.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/aa4be45d2215558d82bb69d67864b63b988a2a68bcb351d9df525cd2f092e1b2.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/c042306e4f4ade72ec09970751e90d17b81e92b5271ae8fd855bcfb619c58fb3.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/f40c3add0b1e6bad8a7be9baa55c77febac5563ac43f2d4582cced20a6392dd9.jpg)


# The Eigenvalue Problem

# Error Analysis

# Conclusions...

# Low modes

For fixed $k$ , $\Delta x \to 0$ : 

$$
\begin{array}{l} \hat {\lambda} ^ {k} = \frac {2}{\Delta x ^ {2}} \{1 - \cos (k \pi \Delta x) \} \\ = \frac {2}{\Delta x ^ {2}} \left\{1 - \left(1 - \frac {1}{2} k ^ {2} \pi^ {2} \Delta x ^ {2} + O \left(\Delta x ^ {4}\right)\right) \right\} \\ = k ^ {2} \pi^ {2} + \mathcal {O} \left(\Delta x ^ {2}\right) \\ \end{array}
$$

$$
\text {s e c o n d - o r d e r c o n v e r g e n c e}, \hat {\lambda} ^ {k} \rightarrow \lambda^ {k}.
$$

# The Eigenvalue Problem

# Error Analysis

# ...Conclusions...

# High modes:

For $k = n$ 

$$
\Delta x = \frac {1}{n + 1}
$$

$$
\begin{array}{l} \hat {\lambda} ^ {n} = \frac {2}{\Delta x ^ {2}} \{1 - \cos (\frac {n \pi}{n + 1}) \} \\ = 4 (n + 1) ^ {2} \quad \text {a s} \quad \Delta x \rightarrow 0 \\ \neq n ^ {2} \pi^ {2} = \lambda^ {n}. \\ \end{array}
$$

High modes $(k\approx n)$ are not accurate. 

# The Eigenvalue Problem

# Error Analysis

...Conclusions... 

Low modes vs. high modes 

Example: $n = 19, \Delta x = 1 / 20$ 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/7a638da5254f67211422e7031407d1a81af83baab6059bce7fcc9a336c56dedc.jpg)


![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/1ae3c0a9dd7b208c6ba81f6bbc58a310d44194c28c1143edc8d1b706aaab8e52.jpg)


# The Eigenvalue Problem

# Error Analysis

...Conclusions... 

Low modes vs. high modes 

$$
k \ll n
$$

$$
k \approx n
$$

N19 

$\hat{u}^{k}$ resolved 

$\hat{\lambda}^{k}$ accurate 

$$
\hat {\lambda} ^ {k} - \lambda^ {k} \sim \mathcal {O} (\Delta x ^ {2})
$$

$\hat{\boldsymbol{u}}^k$ not resolved 

$\hat{\lambda}^k$ not accurate 

$\hat{\lambda}^{k}-\lambda^{k}$ is $\mathcal{O}(1)$ 

BUT: as $\Delta x \to 0$ , $n \to \infty$ , so any fixed mode $k$ converges. 

# The Eigenvalue Problem

# Error Analysis

...Conclusions... 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/68e37ae3c0b757e7bde927163653b056e7855a25013c13c21120310bc3b83ec3.jpg)


# The Eigenvalue Problem

# Error Analysis

...Conclusions 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/15dfa48f36154ff4c3ad81aee3ac2922846b84f095ccbb7d10949d743cb84cfc.jpg)


# The Eigenvalue Problem

# Condition Number of $A$

For a SPD matrix $M$ , the condition number $\kappa_{M}$ is given by 

$$
\kappa_ {M} = \frac {\text {m a x i m u m e i g e n v a l u e o f M}}{\text {m i n i m u m e i g e n v a l u e o f M}}.
$$

Thus, for our $A$ matrix, 

$$
\kappa_ {A} \rightarrow \frac {4 n ^ {2}}{\pi^ {2}} \mathrm {a s} \Delta x \rightarrow 0
$$

grows (in $\mathbb{R}^{1}$ ) as number of grid points squared. N20 

Importance: understanding solution procedures. 

# Eigenvalue Problem in 2D

# Statement

$$
- \nabla^ {2} u = \lambda u \quad \text {i n} \Omega
$$

$$
u = 0 \quad \text {o n} \Gamma
$$

Assume (for simplicity) 

$$
L _ {x} = L _ {y} = 1
$$

Solutions $(u(x,y),\lambda)$ 

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/a958efed-5598-47bd-8351-732d0ca6a589/c2fcb9094e2d2bf7022e97d20fdd0aa7fb14e47c7521d36ad036585ca8dd31ab.jpg)


# Eigenvalue Problem in 2D

# Exact Solution

# Eigenvalues

$$
\begin{array}{l} u ^ {k, l} (x, y) = \sin (k \pi x) \sin (l \pi y) \\ - \nabla^ {2} u ^ {k, l} = \left(k ^ {2} \pi^ {2} + l ^ {2} \pi^ {2}\right) u ^ {k, l} \\ \end{array}
$$

# Eigenvectors

$$
\lambda^ {k, l} = k ^ {2} \pi^ {2} + l ^ {2} \pi^ {2}, \quad k, l = 1, \dots
$$

# Eigenvalue Problem in 2D

# Discrete Problem

# Eigenvectors

$$
A \hat {u} = \hat {\lambda} \hat {u} \quad \Rightarrow \quad \left(\underline {{u}} ^ {k, l}, \hat {\lambda} ^ {k, l}\right)
$$

$$
\begin{array}{l} u _ {i, j} ^ {k, l} = \sin (k \pi x _ {i}) \sin (l \pi y _ {j}) \\ = \sin (k \pi i \Delta x) \sin (l \pi j \Delta y) \\ = \sin (\frac {k \pi i}{n + 1}) \sin (\frac {l \pi j}{m + 1}) \\ \end{array}
$$

$$
k, i = 1, \dots , n \quad l, j = 1, \dots , m
$$

# Eigenvalue Problem in 2D

# Discrete Problem

# Eigenvalues

$$
\hat {\lambda} ^ {k, l} = \frac {2}{\Delta x} \{1 - \cos (k \pi \Delta x) \} + \frac {2}{\Delta y} \{1 - \cos (l \pi \Delta y) \}
$$

Low Modes 

$$
\Delta x, \Delta y \rightarrow 0 (k, l \text {f i x e d})
$$

$$
\begin{array}{l} \hat {\lambda} ^ {k, l} = k ^ {2} \pi^ {2} + l ^ {2} \pi^ {2} \\ + \mathcal {O} (\Delta x ^ {2}, \Delta y ^ {2}) \\ \end{array}
$$

High Modes 

$$
k \approx n, l \approx m
$$

$$
\begin{array}{l} \hat {\lambda} ^ {k, l} = 4 (n + 1) ^ {2} + 4 (m + 1) ^ {2} \\ \mathrm {a s} \Delta x, \Delta y \rightarrow 0 \\ \end{array}
$$

# Condition Number of $A$

# Eigenvalue Problem in 2D

$$
\kappa_ {A} \rightarrow \frac {4 n ^ {2} + 4 m ^ {2}}{2 \pi^ {2}} \quad \text {a s} \Delta x, \Delta y \rightarrow 0
$$

$$
\mathrm {I f} m \approx n
$$

$$
\kappa_ {A} \rightarrow \frac {4 n ^ {2}}{\pi^ {2}}
$$

grows (in $\mathbb{R}^2$ ) as number of grid points. 

(better than in 1D, relatively speaking !!) 