# 第七章 数学物理方程的变分原理

# 变分问题介绍

# 1.1 古典变分问题

例1.1（最速降线问题）这是Bernoulli1696年提出的问题。如图7.1，设点A（0，0）和B（x1,y1）不同在一条与y轴平行的直线上。有一质点受重力作用从A到B沿曲线路径自由下滑，求质点下降最快的路径。问题中不考虑摩擦阻力。

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/42856561caafd102ef68a18de98143b951ae3dbadf55d5ad8dac9d1bc6944906.jpg)


下降最快即所需时间最短，

从A到B沿任一光滑曲线

$$
l: y = y (x), 0 \leq x \leq x _ {1}
$$

下滑所需时间：

# 1.1 古典变分问题

例1.1（最速降线问题）这是Bernoulli1696年提出的问题。如图7.1，设点A（0，0）和B（x1,y1）

不同在一条与y轴平行的直线上。有一质点受重力作用从A到B沿曲线路径自由下滑，求质点下降最快的路径。问题中不考虑摩擦阻力。

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/67f8b7d95fc2b0f5c525eabf958fb90ed451947cf3cb6cdca9159dda42d54337.jpg)


$$
l: y = y (x), 0 \leq x \leq x _ {1}
$$

从A点至曲线上任意一点P(x,y)速率达到了 dsv =dt

质点质量为m

# 1.1 古典变分问题

例1.1（最速降线问题）这是Bernoulli1696年提出的问题。如图7.1，设点A（0，0）和B（x1,y1）不同在一条与y轴平行的直线上。有一质点受重力作用从A到B沿曲线路径自由下滑，求质点下降最快的路径。问题中不考虑摩擦阻力。

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/834a9a652c50c56421097e0a972b27fbdaca8e7e69d519440d5c9689fa9dd01e.jpg)


从A点至曲线上任意一点P(x,y)速率达到了

,ds v P mgy= ,到达 点时失去的势能为 dt

1得到的动能为 2mv 。2

# 1.1 古典变分问题

# 例1.1（最速降线问题）这是

# Bernoulli1696年提出的问题。

如图7.1，设点A（0，0）和B（x1,y1）

ds从A点至曲线上任意一点P(x,y)

, v P mgydt= ,到达 点时失去的势能为 速率达到了

1 2得到的动能为 mv 。2

根据能量守恒原理有 ： $\frac { 1 } { 2 } m ( \frac { d s } { d t } ) ^ { 2 } = m g y$ 

$$
\sqrt {1 + y ^ {\prime 2}} \frac {d x}{d t} = \frac {d s}{d t} = \sqrt {2 g y}
$$

$$
d t = \sqrt {\frac {1 + y ^ {\prime 2}}{2 g y}} d x
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/8b811ef1c7428d8500c86951dd3f80b7889b72ce065fb4268a5050770cc612e2.jpg)


# 1.1 古典变分问题

从A沿l下滑至B所需时间为

$$
d t = \sqrt {\frac {1 + y ^ {\prime 2}}{2 g y}} d x
$$

$$
T = \int_ {0} ^ {x _ {1}} \sqrt {\frac {1 + y ^ {\prime 2}}{2 g y}} d x \tag {1.1}
$$

T y y x 是一个包含函数 = ( )的积分式

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/207d87f85b913fc374175b73d60aa0b41c56ba6ee750e50d6eefa4edfec9404f.jpg)


当y在某一个函数集合K中取定一个函数

从（1.1）式就得到一个确定的实数值T。

称T是K上的一个泛函，记成

$$
T (y) = \int_ {0} ^ {x _ {1}} \sqrt {\frac {1 + y ^ {\prime 2}}{2 g y}} d x
$$


即 $T=T(y)$ 是定义在函数集合 $K$ 上的泛函，可取（例如）

$$
K = \{ y \mid y\in C^{1}[0,x_1],\ y(0)=0,\ y(x_1)=y_1\}.
$$

最速降线问题可写成：在集合 $K$ 上最小化泛函 $T(y)$ ，即找一个 $y_0\in K$ 使 $T(y_0)=\min\limits_{y\in K}T(y)$ 。

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/a576be9eefdf8d6181a7c470cfde47945b3aef996415a09349f6d0ede7e22c48.jpg)


例1.2（最小曲面问题）如图7.2：在 $xy$ 平面上给定有界开区域 $\Omega$ 与其边界 $\partial\Omega$，在 $\partial\Omega$ 上给定 $u \big |_{\partial\Omega}=\varphi(x,y)$ ，其中 $\varphi$ 为已知函数。

这样就给出了三维空间中的一条封闭曲线C，最小曲面问题就是求张紧在曲线C的曲面中，其面积最小的曲面。

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/a5f9adb1c76b7816d8f67b1239f919369d32331917a13158aedc3b1d5fea3490.jpg)


记曲面 $u=u(x,y)$ ， $(x,y)\in \overline\Omega$

$$
u = u (x, y), (x, y) \in \overline {{\Omega}}
$$

对应u的曲面面积为

$$
S (\boldsymbol {u}) = \iint_ {\Omega} \sqrt {1 + \left(\frac {\partial u}{\partial x}\right) ^ {2} + \left(\frac {\partial u}{\partial y}\right) ^ {2}} d x d y \tag {1.4}
$$

例1.3（等周问题）在长度为一定的所有平面光滑封闭曲线中求所围面积为最大的曲线（图7.3）。

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/8d2d6e2a6556db3ca2369dceccb893542bba3ed99b2c6034e5a4c00326c0e251.jpg)


设曲线方程为 $\left\{ { \begin{array} { l } { x = x ( s ) , } \\ { y = y ( s ) , } \end{array} } \right.$ 

$$
s _ {1} \leq s \leq s _ {2}
$$

1 2 1 2且满足x s x s y s y s( ) ( ), ( ) ( ).= =

2 1 s s dx dy ds lds ds         + =    等周问题是在条件 下，

在函数集合

$$
\begin{array}{l} K = \left\{(x, y) \mid x, y \in C ^ {1} \left[ s _ {1}, s _ {2} \right], x \left(s _ {1}\right) = x \left(s _ {2}\right), y \left(s _ {1}\right) = y \left(s _ {2} \right. \right. \\ \max  s (x, y) = \frac {1}{2} \int_ {s _ {1}} ^ {s _ {2}} \left(x \frac {d y}{d s} - y \frac {d x}{d s}\right) d s, \\ \end{array}
$$

变分问题：

在一个函数集合中求泛函的极小或极大的问题。

函数集合K是一元函数的集合（例1.1 ）

(例1.1 )从A沿l下滑至B所需时间为

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/1a671bed9964350f0d2f0d7fa5418b38b65b527f907f8a4947f060669dc7599c.jpg)


$$
T (y) = \int_ {0} ^ {x _ {1}} \sqrt {\frac {1 + y ^ {\prime 2}}{2 g y}} d x
$$

0 , y K 求 使得0 = min ( )y Ky T y T ( ) . (1.3)

变分问题：

在一个函数集合中求泛函的极小或极大的问题。

K是一个多元函数（例1.2）的集合

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/5e6abf1aa800724015748efb727b6a91f211a6cdd841da3693be7a34564d434b.jpg)


$$
S (u) = \iint_ {\Omega} \sqrt {1 + \left(\frac {\partial u}{\partial x}\right) ^ {2} + \left(\frac {\partial u}{\partial y}\right) ^ {2}} d x d y \tag {1.4}
$$

$$
\left\{ \begin{array}{l} {\text {求} u _ {0} \in K, \text {使 得}} \\ {S \left(u _ {0}\right) \leq S \left(u\right), \forall u \in K.} \end{array} \right. \tag {1.5}
$$

变分问题：

在一个函数集合中求泛函的极小或极大的问题。

K是多个一元函数（例1.3）的集合

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/68b3924b-76f3-4f9b-8520-04dd5957c8d6/818b4867d95a226a8d34447a23f2bc906121c226b9d4f705c2f78525f9a6b6d0.jpg)


dxds    dyds 2等周问题是在条件 s ds l = 下，s在函数集合

$$
\max  s (x, y) = \frac {1}{2} \int_ {s _ {1}} ^ {s _ {2}} \left(x \frac {d y}{d s} - y \frac {d x}{d s}\right) d s,
$$

# 1.2 变分问题解的必要条件

$$
C _ {0} ^ {1} [ a, b ] = \left\{v \mid v \in C ^ {1} [ a, b ], v (a) = 0, v (b) = 0 \right\}.
$$

引理1.1（变分法的基本引理）

设u C a b  , ，  1 0 ( ) ( ) 0, ,且 u x v x dx v C a b=  

则有u x( ) 0 a x b( ) 

证 ：

设有  ( , ), a b 使u( ) 0,  不妨设u( ) 0. 

由u的连续性，存在   的一 个邻 域( 1 2, ( , ),)  a b

1 2    , ( ) u x在 上有  0.

# 1.2 变分问题解的必要条件

$$
C _ {0} ^ {1} [ a, b ] = \left\{\nu | \nu \in C ^ {1} [ a, b ], \nu (a) = 0, \nu (b) = 0 \right\}.
$$

引理1.1（变分法的基本引理）

设u C a b   , ，  1 0 ( ) ( ) 0, ,a 且 u x v x dx v C a b=  

则有 $u ( x ) \equiv \mathbf { 0 }$ $( a \leq x \leq b )$ 

证 ：设有 $\xi \in ( a , b )$ ,  使u( ) 0,   不妨设u( ) 0. 

由u的连续性，存在   的一 个邻 域( 1 2, ( , ),)  a b1 2   , ( ) 0.u x 在 上有

取函数 0, 1,a x 2 21 2( ) ( ) ( ) ,v x x x = − − 1 2,x  (1.6)0, 2 .x b 

# 1.2 变分问题解的必要条件

$$
C _ {0} ^ {1} [ a, b ] = \left\{\nu | \nu \in C ^ {1} [ a, b ], \nu (a) = 0, \nu (b) = 0 \right\}.
$$

引理1.1（变分法的基本引理）

设u C a b   , ，  1 0 ( ) ( ) 0, ,a b且 u x v x dx v C a b=  

则有u x a x b( ) 0 ( )  

$\begin{array} { r l } & { \overset { \triangledown } { \boldsymbol { \mathcal { X } } } \quad \mathbf { \mathcal { V } } ( \mathbf { \mathcal { X } } ) = \left\{ \begin{array} { l l } { \mathbf { 0 } , } \\ { ( x - \xi _ { 1 } ) ^ { 2 } ( x - \xi _ { 2 } ) ^ { 2 } , } \\ { \mathbf { 0 } , } \end{array} \right. } \\ & { \left[ \mathbf { \mathcal { a } } , b \right] , } \end{array}$ 1,a x 取函 1 2,x  (1.6)10v C 2 .x b 

( ) ( )b u x v x dx =而且 ob 2 2 1u x x ( )( ) −  2 2( ) x dx−   01

与假设条件矛盾，所以在a b u x, ( )有  0

# 1.2 变分问题解的必要条件

$$
C _ {0} ^ {1} [ a, b ] = \left\{\nu | \nu \in C ^ {1} [ a, b ], \nu (a) = 0, \nu (b) = 0 \right\}.
$$

引理1.1（变分法的基本引理）

设u C a b   , ，  1 0 ( ) ( ) 0, ,a 且 u x v x dx v C a b=  

则有u x a x b( ) 0 ( )  

   1 0 0 , ,n 引理中的条件v C a b v C a b 若改为 ，

结论仍然成立，其中

$$
C _ {0} ^ {n} [ a, b ] = \left\{\nu | \nu \in C ^ {n} [ a, b ], \nu^ {(k)} (a) = \nu^ {(k)} (b) = 0, k = 0, 1, \mathrm {L}, n - 1 \right\}.
$$

它的证明只要将（1.6）式中 v x( ) 的表达式的

幂指数2改为n+1或n+2

$$
C _ {0} ^ {1} [ a, b ] = \left\{\nu | \nu \in C ^ {1} [ a, b ], \nu (a) = 0, \nu (b) = 0 \right\}.
$$

讨论最简单的变分问题,推导其解满足的必要条件。函数集合

$$
K = \left\{y \mid y \in C ^ {1} [ a, b ], y (a) = y _ {a}, y (b) = y _ {b} \right\} \tag {1.7}
$$

考虑y和y有关的泛函 ( ) ( , , )J y F x y y dx=  (1.8)F为对各自变量的偏导数均连续的函数。

变分问题（Ⅰ）： 求y K  ，使得J y J w w K( ) ( ),  设y是变分问题（Ⅰ）的解，则 y K   10对一切C a b,因  ( ) ( ) 0, a b= = 有y K +   ， 其中任意实数 R。

函数集合

$$
K = \left\{y \mid y \in C ^ {1} [ a, b ], y (a) = y _ {a}, y (b) = y _ {b} \right\} \tag {1.7}
$$

考虑y和y有关的泛函 ( ) ( , , )J y F x y y dx=  (1.8)

F为对各自变量的偏导数均连续的函数。

变分问题（Ⅰ）： 求y K  ，使得J y J w w K( ) ( ), .  

设y是变分问题（Ⅰ）的解，  10则y K C a b 对一切 ,

因  ( ) ( ) 0, a b= = 有y K +   ，其中任意实数 R。

因y J y 是使 ( )达到极小的函数，所以有

$$
J (y) \leq J (y + \alpha \eta), \quad \forall \alpha \in ;. \tag {1.9}
$$

$$
J (y) = \int_ {a} ^ {b} F \left(x, y, y ^ {\prime}\right) d x \tag {1.8}
$$

因y J y 是使 ( )达到极小的函数，所以有

$$
J (y) \leq J (y + \alpha \eta), \quad \forall \alpha \in ;. \tag {1.9}
$$

把J y( ) +  看成 的一元函数，记

$$
\varphi (\alpha) = J (y + \alpha \eta) = \int_ {a} ^ {b} F (x, y (x) + \alpha \eta (x), y ^ {\prime} (x) + \alpha \eta^ {\prime} (x)) d x
$$

把函数 的一阶导数值称为泛函J的一阶变分， 记为 J。 2 设y C a b ( , )，有

$$
\delta J = \frac {d \varphi}{d \alpha} \bigg | _ {\alpha = 0} = \int_ {a} ^ {b} \left(\frac {\partial F (x , y , y ^ {\prime})}{\partial y} \eta + \frac {\partial F (x , y , y ^ {\prime})}{\partial y ^ {\prime}} \eta^ {\prime}\right) d x
$$

$$
J (y) = \int_ {a} ^ {b} F \left(x, y, y ^ {\prime}\right) d x \tag {1.8}
$$

因y J y 是使 ( )达到极小的函数，所以有

$$
\begin{array}{l} J (y) \leq J (y + \alpha \eta), \quad \forall \alpha \in ;. \tag {1.9} \\ \varphi (\alpha) = J (y + \alpha \eta) = \int_ {a} ^ {b} F (x, y (x) + \alpha \eta (x), y ^ {\prime} (x) + \alpha \eta^ {\prime} (x)) d x \\ \delta J = \frac {d \varphi}{d \alpha} \bigg | _ {\alpha = 0} = \int_ {a} ^ {b} \left(\frac {\partial F (x , y , y ^ {\prime})}{\partial y} \eta + \frac {\partial F (x , y , y ^ {\prime})}{\partial y ^ {\prime}} \eta^ {\prime}\right) d x \\ \end{array}
$$

有（1.9）式可推出

   ( ) =0 min 在 有 由一元函数取得极值的必要条件得

$$
\int_ {a} ^ {b} \left(\frac {\partial F}{\partial y} \eta + \frac {\partial F}{\partial y ^ {\prime}} \eta^ {\prime}\right) d x = 0, \quad \forall \eta \in C _ {0} ^ {1} [ a, b ]
$$

$$
K = \left\{y \mid y \in C ^ {1} [ a, b ], y (a) = y _ {a}, y (b) = y _ {b} \right\} \tag {1.7}
$$

变分问题（Ⅰ）： 求y K  ，使得J y J w w K( ) ( ),  

$$
J (y) \leq J (y + \alpha \eta), \quad \forall \alpha \in ;. \tag {1.9}
$$

有（1.9）式可推出

$$
\int_ {a} ^ {b} \left(\frac {\partial F}{\partial y} \eta + \frac {\partial F}{\partial y ^ {\prime}} \eta^ {\prime}\right) d x = 0, \quad \forall \eta \in C _ {0} ^ {1} [ a, b ]
$$

第二项经过分部积分， 可得到

$$
\int_ {a} ^ {b} \left[ \frac {\partial F}{\partial y} - \frac {d}{d x} \left(\frac {\partial F}{\partial y ^ {\prime}}\right) \right] \eta d x = 0, \quad \forall \eta \in C _ {0} ^ {1} [ a, b ] \quad 1. 1 0
$$

# 引理1.1（变分法的基本引理）

设u C a b   , ，  1 0 ( ) ( ) 0, ,a 且 u x v x dx v C a b=  

则有u x a x b( ) 0 ( )  

第二项经过分部积分， 可得到

$$
\int_ {a} ^ {b} \left[ \frac {\partial F}{\partial y} - \frac {d}{d x} \left(\frac {\partial F}{\partial y ^ {\prime}}\right) \right] \eta d x = 0, \quad \forall \eta \in C _ {0} ^ {1} [ a, b ]
$$

根据变分法的基本引理，得到

$$
\frac {\partial F}{\partial y} - \frac {d}{d x} \left(\frac {\partial F}{\partial y ^ {\prime}}\right) = 0 \quad (1. 1 1) \quad \text {E u l e r 方 程}
$$

$$
\frac {\partial F}{\partial y} - \frac {d}{d x} \left(\frac {\partial F}{\partial y ^ {\prime}}\right) = 0 \quad (1. 1 1) \quad \text {E u l e r 方 程}
$$

例1.4设K仍如（1.7）式所示，且

$$
J (y) = \frac {1}{2} \int_ {a} ^ {b} \left[ p (x) \left(\frac {d y}{d x}\right) ^ {2} + q (x) y ^ {2} - 2 f (x) y \right] d x,
$$

求J y K( )在 的极值。

  2 如果取到极值的函数y C a b  , ，则y满足的Euler

方程为 ( ) ( ) ( ) 0d dy p x q x y f xdx dx   + − =   

( ) , ( )a by a y y b y= = 构成一个定解问题。

例1.4设K仍如（1.7）式所示，且

$$
J (y) = \frac {1}{2} \int_ {a} ^ {b} \left[ p (x) \left(\frac {d y}{d x}\right) ^ {2} + q (x) y ^ {2} - 2 f (x) y \right] d x,
$$

求J y K( )在 的极值。

  2 如果取到极值的函数y C a b  , ，则y满足的Euler方程为 ( ) ( ) ( ) 0d dy p x q x y f xdx dx   − + − =    dx

$$
y (a) = y _ {a}, y (b) = y _ {b} \quad {\text {构 成 一 个 定 解 问 题}}.
$$

在例1.4中,  F x y y ( , , ) 式中头两项是

y y J y 及 的平方项 ( )为二次泛函。

$$
J (y) = \int_ {a} ^ {b} F \left(x, y, y ^ {\prime}\right) d x \tag {1.8}
$$

$$
K = \left\{y \mid y \in C ^ {1} [ a, b ], y (a) = y _ {a}, y (b) = y _ {b} \right\} \tag {1.7}
$$

$$
C _ {0} ^ {1} [ a, b ] = \left\{v \mid v \in C ^ {1} [ a, b ], v (a) = 0, v (b) = 0 \right\}.
$$

 10则y K C a b 对一切 ,

变分问题（Ⅰ）： 求y K  ，使得J y J w w K( ) ( ), .  

$$
K _ {1} = \left\{y \mid y \in C ^ {1} [ a, b ], \quad y (a) = y _ {a} \right\} \tag {1.12}
$$

同样求泛函（1.8）式的极小值。

变分问题（Ⅱ）： 1求y K  ， 1使得J y J w w K( ) ( ),   

$$
J (y) = \int_ {a} ^ {b} F \left(x, y, y ^ {\prime}\right) d x \tag {1.8}
$$

变分问题（Ⅱ）： 1 求y K  ， 1使得J y J w w K( ) ( ),   在解的必要条件推导过程中，规定

$$
\begin{array}{l} \eta \in \left\{\nu \mid \nu \in C ^ {1} [ a, b ], \quad \nu (a) = 0 \right\} (1.13) \\ \delta J = \frac {d \varphi}{d \alpha} \bigg | _ {\alpha = 0} = \int_ {a} ^ {b} \left(\frac {\partial F (x , y , y ^ {\prime})}{\partial y} \eta + \frac {\partial F (x , y , y ^ {\prime})}{\partial y ^ {\prime}} \eta^ {\prime}\right) d x \\ = \int_ {a} ^ {b} \left[ \frac {\partial F}{\partial y} - \frac {d}{d x} \left(\frac {\partial F}{\partial y ^ {\prime}}\right) \right] \eta d x + \left[ \eta \frac {\partial F}{\partial y ^ {\prime}} \right] _ {x = b} = 0 (1.14) \\ \end{array}
$$

这个式子对一切满足（1.13）式的 成立。

$$
J (y) = \int_ {a} ^ {b} F (x, y, y ^ {\prime}) d x \quad (1. 8) \quad w = y + \alpha \eta
$$

$$
C _ {0} ^ {1} [ a, b ] = \left\{\nu | \nu \in C ^ {1} [ a, b ], \nu (a) = 0, \nu (b) = 0 \right\}.
$$

变分问题（Ⅱ）： 1求y K  ， 1使得J y J w w K( ) ( ),   在解的必要条件推导过程中，规定

$$
\eta \in \left\{\nu \mid \nu \in C ^ {1} [ a, b ], \quad \nu (a) = 0 \right\} \tag {1.13}
$$

$$
\delta J = \left. \frac {d \varphi}{d \alpha} \right| _ {\alpha = 0} = \int_ {a} ^ {b} \left[ \frac {\partial F}{\partial y} - \frac {d}{d x} \left(\frac {\partial F}{\partial y ^ {\prime}}\right) \right] \eta d x + \left[ \eta \frac {\partial F}{\partial y ^ {\prime}} \right] _ {x = b} = 0 \tag {1.14}
$$

这个式子对一切满足（1.13）式的 成立。

 1 0 对   C a b, 也成立

$$
\begin{array}{l} \int_ {a} ^ {b} \left[ \frac {\partial F}{\partial y} - \frac {d}{d x} \left(\frac {\partial F}{\partial y ^ {\prime}}\right) \right] \eta d x = 0, \quad \forall \eta \in C _ {0} ^ {1} [ a, b ] \quad 1. 1 0 \\ \frac {\partial F}{\partial y} - \frac {d}{d x} \left(\frac {\partial F}{\partial y ^ {\prime}}\right) = 0 (1.11) \\ \delta J = \frac {d \varphi}{d \alpha} \bigg | _ {\alpha = 0} = \int_ {a} ^ {b} \left[ \frac {\partial F}{\partial y} - \frac {d}{d x} \left(\frac {\partial F}{\partial y ^ {\prime}}\right) \right] \eta d x + \left[ \eta \frac {\partial F}{\partial y ^ {\prime}} \right] _ {x = b} = 0 (1.14) \\ \end{array}
$$

得到（1.10）式，同样导出Euler方程（1.11），将它代回（1.14）式，有

$$
\left[ \eta \frac {\partial F}{\partial y ^ {\prime}} \right] _ {x = b} = 0, \quad \forall \eta \in \left\{\nu | \nu \in C ^ {1} [ a, b ], \quad \nu (a) = 0 \right\}
$$

$$
\delta J = \frac {d \varphi}{d \alpha} \bigg | _ {\alpha = 0} = \int_ {a} ^ {b} \left[ \frac {\partial F}{\partial y} - \frac {d}{d x} \left(\frac {\partial F}{\partial y ^ {\prime}}\right) \right] \eta d x + \left[ \eta \frac {\partial F}{\partial y ^ {\prime}} \right] _ {x = b} = 0 \tag {1.14}
$$

得到（1.10）式，同样导出Euler方程（1.11），将它代回（1.14）式，有

$$
\left[ \eta \frac {\partial F}{\partial y ^ {\prime}} \right] _ {x = b} = 0, \quad \forall \eta \in \left\{v \mid v \in C ^ {1} [ a, b ], \quad v (a) = 0 \right\}
$$

由于 的任意性，可得 $\left. { \frac { \partial F } { \partial y ^ { \prime } } } \right| _ { x = b } = 0$ （1.15）

（1.15）式是在右端点 x b = 上的边界条件，这种条件称为自然边界条件。

$$
J (y) = \int_ {a} ^ {b} F \left(x, y, y ^ {\prime}\right) d x \tag {1.8}
$$

变分问题（Ⅱ）： 1使得J y J w w K( ) ( ),   

变分问题（Ⅱ）的解y满足定解问题

$$
\left\{ \begin{array}{l} \frac {\partial F}{\partial y} - \frac {d}{d x} \left(\frac {\partial F}{\partial y ^ {\prime}}\right) = 0, \\ y (a) = y _ {a}, \quad \left. \frac {\partial F}{\partial y ^ {\prime}} \right| _ {x = b} = 0. \end{array} \right.
$$

例1.4设K仍如（1.7）式所示P190，且

$$
J (y) = \frac {1}{2} \int_ {a} ^ {b} \left[ p (x) \left(\frac {d y}{d x}\right) ^ {2} + q (x) y ^ {2} - 2 f (x) y \right] d x,
$$

求J y K( )在 的极值。

$$
\begin{array}{l} \left. \frac {\partial F}{\partial y ^ {\prime}} \right| _ {x = b} = 0 \tag {1.15} \\ J (y) = \frac {1}{2} \int_ {a} ^ {b} \left[ p (x) \left(\frac {d y}{d x}\right) ^ {2} + q (x) y ^ {2} - 2 f (x) y \right] d x, \\ \end{array}
$$

将变分问题（Ⅱ）应用到例1.4的二次泛函 J y( )，

我们得到（1.15） 式为

$$
\left[ p (x) \frac {d y}{d x} \right] _ {x = b} = 0
$$

规定 0p x p a x b( ) 0( ),   

在 x b= 的自然边界条件就是 $\left. { \frac { d y } { d x } } \right| _ { x = b } = 0$ 

第二类边界条件。

考虑依赖于多个函数的泛函p192

$$
J (y _ {1}, \mathbf {L}, y _ {n}) = \int_ {a} ^ {b} F (x, y _ {1}, \mathbf {L}, y _ {n}; y _ {1} ^ {\prime}, \mathbf {L}, y _ {n} ^ {\prime}) d x
$$

或者是依赖于高阶导数的泛函，例如

$$
J (y) = \int_ {a} ^ {b} F \left(x, y, y ^ {\prime}, y ^ {\prime \prime}\right)
$$

也可以考虑依赖于多元函数的泛函，例如

$$
J (u) = \iint_ {\Omega} F (x, y, u, u _ {x}, u _ {y}) d x d y
$$

可以分别讨论其极值必要条件，即Euler方程

# 1.3 ¡ n 中的变分问题

设向量 1 ( , , )T n n x x x =  L ¡ ，矩阵 A aij =    ¡ ,n n设det 0A →，方程组Ax b= 有惟一解，其中已知的向量b ¡ n

命题1.1 n x  ¡ ，满足 Ax b− = 0。

命题1.2 $\left( A x - b , y \right) = \mathbf { 0 } , \forall y \in \mathfrak { i }$ n

如果命题1.1成立，则命题1.2成立

若命题1.2成立成立，令 ,iy e= 1, , ,i n= L

y取各坐标轴上的单位向量，则可推出 Ax b− =0

命题1.1和命题1.2是相互等价的命题。

$$
J (\boldsymbol {x}) = \frac {1}{2} \sum_ {i, j = 1} ^ {n} a _ {i j} x _ {i} y _ {i} - \sum_ {i = 1} ^ {n} b _ {i} x _ {i}
$$

矩阵形式为 $J ( x ) = \frac { 1 } { 2 } ( A x , x ) - ( b , x )$ 设A对称，有 ( , ) ( , ), Ax y Ay x a=  ¡

$$
\begin{array}{l} J (x + \alpha y) = \frac {1}{2} (A x + \alpha A y, x + \alpha y) - (b, x + \alpha y) \\ = J (x) + \alpha (A x - b, y) + \frac {\alpha^ {2}}{2} (A y, y) \\ \end{array}
$$

命题1.3  n x  ¡ 使得 ( ) ( ),J x J y ny  ¡若命题1.3成立，即x使 J x( ) 达到极小，

命题1.2 $\left( A x - b , y \right) = \mathbf { 0 } , \forall y \in \mathfrak { i }$ n

$$
J (x + \alpha y) = J (x) + \alpha (A x - b, y) + \frac {\alpha^ {2}}{2} (A y, y)
$$

命题1.3  n x  ¡ 使得 ( ) ( ),J x J y ny  ¡

若命题1.3成立，即x使 J x( ) 达到极小，

则 的函数 在 取得极小值。

此时有 0( ) | 0d J x y + ==d

( , ) 0,Ax b y− = ny  ¡ 命题1.2成立。

命题1.3  n x  ¡ 使得 ( ) ( ),J x J y ny  ¡

命题1.2  ( , 0,)Ax b y− = ny  ¡

$$
J (x + \alpha y) = J (x) + \alpha (A x - b, y) + \frac {\alpha^ {2}}{2} (A y, y)
$$

设A对称正定，即 ( , ) 0,Ax x  0.x 

如果命题1.2成立， 就有

$$
J (x + \alpha y) = J (x) + \frac {\alpha^ {2}}{2} (A y, y), \quad \forall y \in ; ^ {n}, \alpha \in ;
$$

如果令 ，就有 ( ) ( ),J x J w nw  ¡

命题1.3 成立

命题1.1 n x  ¡ ，满足 Ax b− = 0。

命题1.2  ( , 0,)Ax b y− = ny  ¡

命题1.3  n x  ¡ 使得 ( ) ( ),J x J y ny  ¡

小结：若A对称则命题1.3可推出命题1.2。

如果A对称正定，则命题1.2可推出命题1.3，不论A是否对称，命题1.1和命题1.2总是等价的。如果A对称正定，则3个命题都是等价的。

关于矩阵的特征值问题，即求数 和非零向量

$$
x \in ; ^ {n}, \quad A x = \lambda x
$$

也有对应的变分问题。如果A对称正定，则A的

特征值是实的，满足 1 2 0     L n

kx 是对应 的特征向量

定义关于A的Rayleigh商为

$$
R (x) = \frac {(A x , x)}{(x , x)}, \quad x \in ; ^ {n}, \quad x \neq 0
$$

R x( ) 的最大、最小值对应A的最大、最小特征值

$$
\lambda_ {1} = R \left(x _ {1}\right) = \max  _ {x \neq 0} R (x) \quad R (x) = \frac {(A x , x)}{(x , x)}
$$

其他的特征值（k>1） ( )k k= = max ( )0x

$$
(x _ {1}, x) = \mathbf {L} = (x _ {k - 1}, x) = \mathbf {0}
$$

$$
\lambda_ {n} = R \left(x _ {n}\right) = \min  _ {x \neq 0} R (x)
$$

最大特征值 可以看成条件( , ) 1x x = 下求

( , ) Ax x 极值的条件极值问题的解，

可以看成条件 ( , ) 1 x x = 和 1( , ) 0x x =

下求 ( , ) Ax x 极值的条件极值问题的解。