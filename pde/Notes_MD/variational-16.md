# 2 一维数学物理问题的变分问题

弦振动方程

$$
\rho \frac {\partial ^ {2} u}{\partial t ^ {2}} = T \frac {\partial ^ {2} u}{\partial x ^ {2}} + F (x , t).
$$

其中 $T$ 为弦的张力， $F(x,t)$ 为外力密度。下面考虑平衡问题：若外力 $F$ 与边界条件都与时间无关，

分析弦平衡时的位移u(x)有两种途径:

1 利用力的平衡关系,知u(x)满足微分方程

$$
- T \frac {d ^ {2} u}{d x ^ {2}} = f (x) \quad \text {和 边 界 条 件} u (0) = u (l) = 0
$$

# 2 一维数学物理问题的变分问题

# 分析弦平衡时的位移u(x)有两种途径:

1 利用力的平衡关系,知u(x)满足微分方程

$$
- T \frac {d ^ {2} u}{d x ^ {2}} = f (x) \quad \text {和 边 界 条 件} u (0) = u (l) = 0
$$

$$
\left( \begin{array}{l} P _ {1} \end{array} \right) \quad \left\{ \begin{array}{l l} - T \frac {d ^ {2} u}{d x ^ {2}} = f (x) & \\ u (0) = u (l) = 0 & \end{array} \right. \tag {1}
$$

2 利用能量守恒定律，即利用力学中的极小位能原理或虚功原理也可得到（1）和（2）式。

# 2 一维数学物理问题的变分问题

# 利用能量方法讨论弦平衡问题

# 1 极小位能原理：

一个物体在静态平衡时具有极小位能

弦变形前（零位移近似）可取 $W _ {1} = \displaystyle \int_ {0} ^ {l} T \, d x$ ；变形后弧长近似上式给出弦变形后的能量

$$
W _ {2} = \int_ {0} ^ {s} T d s = \int_ {0} ^ {l} T \left(1 + u _ {x} ^ {2}\right) ^ {1 / 2} d x \approx \int_ {0} ^ {l} T \left(1 + \frac {1}{2} u _ {x} ^ {2}\right) d x
$$

弦的应变能（内能）和外力作功分别为

$$
W _ {\text {内}} = W _ {2} - W _ {1} = \frac {1}{2} \int_ {0} ^ {l} T u _ {x} ^ {2} d x \quad W _ {\text {外}} = \int_ {0} ^ {l} f u d x
$$

# 2 一维数学物理问题的变分问题

# 利用能量方法讨论弦平衡问题

# 一、 极小位能原理：

一个物体在静态平衡时具有极小位能弦的应变能（内能）和外力作功分别为

$$
W _ {\text {内}} = W _ {2} - W _ {1} = \frac {1}{2} \int_ {0} ^ {l} T u _ {x} ^ {2} d x \quad W _ {\text {外}} = \int_ {0} ^ {l} f u d x
$$

弦的总位能

$$
J (u) = W _ {\text {内}} - W _ {\text {外}} = \frac {1}{2} \int_ {0} ^ {l} T u _ {x} ^ {2} d x - \int_ {0} ^ {l} f u d x
$$

弦平衡时J(u)最小。

# 2 一维数学物理问题的变分问题弦的总位能

$$
J (\boldsymbol {u}) = W _ {\text {内}} - W _ {\text {外}} = \frac {1}{2} \int_ {0} ^ {l} T u _ {x} ^ {2} d x - \int_ {0} ^ {l} f u d x
$$

弦平衡时J(u)最小。

弦的平衡问题的变分提法：存在 $u \in K$ ，使泛函达到最小：
$$
K = \left\{ v \in C ^{1} [ 0, l ]\mid v (0) = v (l) = 0 \right\},\qquad 
J (u)=\min_{v\in K}\left\{ \frac{1}{2} \displaystyle \int_ {0}^{l} T v _ {x}^{2}\, dx - \int_ {0}^{l} fv\, dx\right\}.
$$

# （承接）变分提法与欧拉方程的联系

$$
F = \frac {1}{2} T u _ {x} ^ {2} - f u,
$$

设拉格朗日密度 \(F=\frac{1}{2}Tu_{x}^2-fu\) ，则由变分极小化 \(J\) 得到相应的 Euler–Lagrange 方程即 \((P_1)\) 的微分方程形式。

$$
\left( \begin{array}{l} P _ {1} \end{array} \right) \quad \left\{ \begin{array}{l l} - T \frac {d ^ {2} u}{d x ^ {2}} = f (x) & \\ u (0) = u (l) = 0 & \end{array} \right. \quad \left( \begin{array}{l} 1 \\ 2 \end{array} \right)
$$

$$
\begin{array}{r} {(P _ {3}) \left\{ \begin{array}{l l} {\text {求} u \in K, \text {使 得} J (u) =  \frac {1}{2} \int_ {0} ^ {l} T u _ {x} ^ {2} d x - \int_ {0} ^ {l} f u d x \text {最 小}} \\ {K = \left\{u \mid u \in C ^ {1} [ 0, l ], u (0) = u (l) = 0 \right\}} \end{array} \right.} \end{array}
$$

设极小解为u

$$
\varphi (\alpha) = J (u + \alpha \eta), \forall \alpha \in R, \forall \eta \in K _ {0}
$$

$$
K _ {0} = \left\{v \mid v \in C ^ {1} [ 0, l ], v (0) = v (l) = 0 \right\}
$$

$$
\varphi (\alpha) = J (u + \alpha \eta)
$$

$$
= \frac {1}{2} \int_ {0} ^ {l} T \left(u _ {x} + \alpha \eta_ {x}\right) ^ {2} d x - \int_ {0} ^ {l} f (u + \alpha \eta) d x
$$

$$
\varphi^ {\prime} (\alpha) \mid_ {\alpha = 0} = \int_ {0} ^ {l} T u _ {x} \eta_ {x} d x - \int_ {0} ^ {l} f \eta d x = 0 \tag {3}
$$

$$
\varphi^ {\prime \prime} (\alpha) \mid_ {\alpha = 0} = \int_ {0} ^ {l} T \left(\eta_ {x}\right) ^ {2} d x > 0 \tag {4}
$$

$$
\begin{array}{r} {(P _ {3}) \left\{ \begin{array}{l l} {\text {求} u \in K, \text {使 得} J (u) =  \frac {1}{2} \int_ {0} ^ {l} T u _ {x} ^ {2} d x - \int_ {0} ^ {l} f u d x \text {最 小}} \\ {K = \left\{u \mid u \in C ^ {1} [ 0, l ], u (0) = u (l) = 0 \right\}} \end{array} \right.} \end{array}
$$

设极小解为u

$$
\varphi (\alpha) = J (u + \alpha \eta), \forall \alpha \in R, \forall \eta \in K _ {0}
$$

$$
K _ {0} = \left\{v \mid v \in C ^ {1} [ 0, l ], v (0) = v (l) = 0 \right\}
$$

$$
\varphi^ {\prime} (\alpha) \mid_ {\alpha = 0} = \int_ {0} ^ {l} T u _ {x} \eta_ {x} d x - \int_ {0} ^ {l} f \eta d x = 0 \tag {3}
$$

$$
\rightarrow \int_ {0} ^ {l} \left[ - T \frac {d ^ {2} u}{d x ^ {2}} - f (x) \right] \eta d x = 0, \forall \eta \in K _ {0} \tag {5}
$$

---变分方程

# 2 一维数学物理问题的变分问题

# 利用能量方法讨论弦平衡问题

# 二、 虚功原理：

若物体所受合力对一切的虚位移所做的功为零，则物体处于静态平衡态。

因弦所受合力为-Td²u -f(x)，利用虚功原理，dx²弦平衡时必有 2 2 00 [ ( )] 0 ,l d u T f x vdx v Kdx − − =  

$$
\rightarrow \int_ {0} ^ {l} T u _ {x} v _ {x} d x = \int_ {0} ^ {l} f v d x, \forall v \in K _ {0}
$$

---虚功方程即变分方程

# 二、 虚功原理：

$$
\int_ {0} ^ {l} T u _ {x} v _ {x} d x = \int_ {0} ^ {l} f v d x, \forall v \in K _ {0}
$$

---虚功方程即变分方程引入二次泛函D(u,v)和线性泛函F(v)的记号

$$
D (u, v) = \int_ {0} ^ {l} T u _ {x} v _ {x} d x, \quad F (v) = \int_ {0} ^ {l} f v d x,
$$

于是虚功方程有下述形式

求u∈K,使得(P）[D(u,v)= F(v),∀v ∈ K,

# 极小位能原理有下述形式

(求,（-最[K={u $\vert u \in C ^ { 1 } [ 0 , l ] , u ( 0 ) = u ( l ) = 0$ }

$$
D (u, v) = \int_ {0} ^ {l} T u _ {x} v _ {x} d x, \qquad F (v) = \int_ {0} ^ {l} f v d x,
$$

[求u∈K,使得J(a)=D(uu)-F(a)最小(P 2[K={u| $u \in C ^ { 1 } [ 0 , l ] , u ( 0 ) = u ( l ) = 0 \big \}$ 

# 2.1 两点边值问题的变分形式

上考虑一般的二阶自伴微分算子L，其定义为

$$
L u = - \frac {d}{d x} \left(p (x) \frac {d u}{d x}\right) + q (x) u \tag {2.1}
$$

  1其中p C a b , ${ | { \mathrm { , } } p ( x ) \geq p _ { 0 } > 0 , q \in C [ a , b ] , q }$ x  ( ) 0.

$$
L u = - \frac {d}{d x} \left(p (x) \frac {d u}{d x}\right) + q (x) u \tag {2.1}
$$

设f ∈c[ab] 对于方程 Lu = f

再加上两端边界条件，就得微分方程的一个两点边值问题

$$
\begin{array}{l l} (P _ {1}) & \left\{ \begin{array}{l l} - \frac {d}{d x} \left(p (x) \frac {d u}{d x}\right) + q (x) u = f (x), & a <   x <   b, \\ u (a) = 0, & \frac {d u (b)}{d x} = 0. \end{array} \right. \end{array}
$$

$$
L u = - \frac {d}{d x} \left(p (x) \frac {d u}{d x}\right) + q (x) u
$$

$$
\begin{array}{l l} (P _ {1}) & \left\{ \begin{array}{l l} - \frac {d}{d x} \left(p (x) \frac {d u}{d x}\right) + q (x) u = f (x), & a <   x <   b, \\ u (a) = 0, & \frac {d u (b)}{d x} = 0. \end{array} \right. \end{array}
$$

1( ) P 中的微分方程写成 Lu- f = 0.

方程两边乘函数 再积分 记函数的内积为

$$
(u, v) = \int_ {a} ^ {b} u (x) v (x) d x \quad \text {则 有} \quad (L u - f, v) = 0
$$

$$
- \int_ {a} ^ {b} v \frac {d}{d x} \left(p \frac {d u}{d x}\right) d x + \int_ {a} ^ {b} q u v d x - \int_ {a} ^ {b} f v d x = 0
$$

$$
- \int_ {a} ^ {b} v \frac {d}{d x} \left(p \frac {d u}{d x}\right) d x + \int_ {a} ^ {b} q u v d x - \int_ {a} ^ {b} f v d x = 0
$$

上式第一项经过分部积分利用边界条件 ( ) 0du b = ，

$$
\left. \int_ {a} ^ {b} p \frac {d u}{d x} \frac {d v}{d x} d x + (p \frac {d u}{d x} v) \right| _ {x = a} + \int_ {a} ^ {b} q u v d x - \int_ {a} ^ {b} f v d x = 0
$$

对函数v x( ) 如果规定 v a( ) 0 = 上式第二项就 不出现。

$$
S _ {0} ^ {1} = \left\{v \mid v \in C ^ {1} [ a, b ], \quad v (a) = 0 \right\}
$$

$$
D (u, v) = \int_ {a} ^ {b} \left(p \frac {d u}{d x} \frac {d v}{d x} + q u v\right) d x \quad F (v) = \int_ {a} ^ {b} f v d x
$$

$$
\begin{array}{l} S _ {0} ^ {1} = \left\{\nu \mid \nu \in C ^ {1} [ a, b ], \quad \nu (a) = 0 \right\} \\ \left. \int_ {a} ^ {b} p \frac {d u}{d x} \frac {d v}{d x} d x + (p \frac {d u}{d x} v) \right| _ {x = a} + \int_ {a} ^ {b} q u v d x - \int_ {a} ^ {b} f v d x = 0 \\ \end{array}
$$

对函数v x( ) 如果规定 上式第二项就 不出现。

$$
D (u, v) = \int_ {a} ^ {b} \left(p \frac {d u}{d x} \frac {d v}{d x} + q u v\right) d x \quad F (v) = \int_ {a} ^ {b} f v d x
$$

是以下问题则

$$
\begin{array}{r l} {\left(P _ {2}\right)} & {\quad \left\{ \begin{array}{l l} {\text {求} u \in S _ {0} ^ {1}, \text {使 得}} \\ {D (u, v) - F (v) = 0, \quad \forall v \in S _ {0} ^ {1}.} \end{array} \right.} \end{array}
$$

$$
S _ {0} ^ {1} = \left\{\nu \mid \nu \in C ^ {1} [ a, b ], \quad \nu (a) = 0 \right\}
$$

是以下问题 $F ( \nu ) = \int _ { a } ^ { b } f \nu d x$ 2.4

10 ,u S求 使 得2 ( )P( , ) ( ) 0,D u v F v − = 10 . v S  

从(2.4)式看到，f是给定的函数，当任取一个函数1v S  时，对应一个实数值

1S0到 的映射 F是一个泛函，它满足

$$
\begin{array}{l} F (v _ {1} + v _ {2}) = F (v _ {1}) + F (v _ {2}) \\ F (c v) = c F (v), \quad \forall c \in \square \\ \end{array}
$$

$$
S _ {0} ^ {1} = \left\{v \mid v \in C ^ {1} [ a, b ], \quad v (a) = 0 \right\}
$$

$$
D (u, v) = \int_ {a} ^ {b} \left(p \frac {d u}{d x} \frac {d v}{d x} + q u v\right) d x \quad F (v) = \int_ {a} ^ {b} f v d x
$$

$$
F (v _ {1} + v _ {2}) = F (v _ {1}) + F (v _ {2}) \quad F (c v) = c F (v), \quad \forall c \in \square
$$

它是一个线性泛函。

同理，D是 1 1 S S 0 0 到 的映射。

u v 或 有一个固定时， D u v ( , ) 是另一者的线性泛函。

我们称D是 上的双线性泛函

D u v D v u ( , ) ( , ) = 称D是 $S _ { \mathbf { 0 } } ^ { 1 } \times S _ { \mathbf { 0 } } ^ { 1 }$ 上的对称双线性泛函

同时，D还满足 D(u,u)≥0, AueS

$$
S _ {0} ^ {1} = \left\{v \mid v \in C ^ {1} [ a, b ], \quad v (a) = 0 \right\}
$$

$\left\{ \begin{array} { l l } { \displaystyle { - \frac { d } { d x } \left( p ( x ) \frac { d u } { d x } \right) + q ( x ) u = f } } \\ { \displaystyle { u ( a ) = 0 } , \quad \displaystyle { \frac { d u ( b ) } { d x } = 0 } . } \end{array} \right.$ ( )p x  ( ), ,x a x b  1( )P 

2 ( )P $\begin{array} { r } { \left\{ \begin{array} { l } { \vec { \mathfrak { X } } u \in S _ { 0 } ^ { 1 } , \not \in \mathcal { Z } ^ { \sharp } \overline { { \mathbb { H } } } ^ { \sharp } } \\ { D ( u , \nu ) - F ( \nu ) = 0 , } \end{array} \right. } \end{array}$ 10 .v S  

如果 u 是问题 也是问题则 u ( ) P 的 解 。

$$
D (u, v) = \int_ {a} ^ {b} \left(p \frac {d u}{d x} \frac {d v}{d x} + q u v\right) d x \partial y \quad F (v) = \int_ {a} ^ {b} f v d x
$$

10 ,u S  求 使 得2( )P( , ) ( ) 0,D u v F v − = 10 .v S 

（ 2 ） 则 有

$$
\begin{array}{l} \int_ {a} ^ {b} p \frac {d u}{d x} \frac {d v}{d x} d x + \int_ {a} ^ {b} q u v d x - \int_ {a} ^ {b} f v d x = 0, \quad \forall v \in S _ {0} ^ {1}. \\ \left. - \int_ {a} ^ {b} v \frac {d}{d x} \left(p \frac {d u}{d x}\right) d x + \left(p \frac {d u}{d x} v\right) \right| _ {x = b} - \left. \left(p \frac {d u}{d x} v\right) \right| _ {x = a} \tag {2.5} \\ + \int_ {a} ^ {b} q u v d x - \int_ {a} ^ {b} f v d x = 0, \quad \forall v \in S _ {0} ^ {1}. \\ \end{array}
$$

$$
\begin{array}{l} L u = - \frac {d}{d x} \left(p (x) \frac {d u}{d x}\right) + q (x) u \quad (u, v) = \int_ {a} ^ {b} u (x) v (x) d x \\ - \int_ {a} ^ {b} v \frac {d}{d x} \left(p \frac {d u}{d x}\right) d x + \left(p \frac {d u}{d x} v\right) \Bigg | _ {x = b} - \left(p \frac {d u}{d x} v\right) \Bigg | _ {x = a} \tag {2.5} \\ + \int_ {a} ^ {b} q u v d x - \int_ {a} ^ {b} f v d x = 0, \quad \forall v \in S _ {0} ^ {1}. \\ \end{array}
$$

 1 10 0 S V v v S v b 的一子集 =  = | , ( ) 0 ,  

2.5对一切 v V 也成立

$$
\int_ {a} ^ {b} (L u - f, v) = 0, \quad \forall v \in C _ {0} ^ {1} [ a, b ].
$$

$$
\begin{array}{r l} {(P _ {1})} & {\left\{ \begin{array}{l l} {- \frac {d}{d x} \bigg (p (x) \frac {d u}{d x} \bigg) + q (x) u = f (x),} & {a <   x <   b,} \\ {u (a) = 0, \quad \frac {d u (b)}{d x} = 0.} & {(u, v) = \int_ {a} ^ {b} u (x) v (x) d x} \end{array} \right.} \end{array}
$$

引理1.1（变分法的基本引理）

设u ∈C[a,b],且["u(x)(x)dx=0,vC[a,b]

则有u x a x b ( ) 0 ( )   

$$
\int_ {a} ^ {b} (L u - f, v) = 0, \quad \forall v \in C _ {0} ^ {1} [ a, b ].
$$

根据变分法的基本引理，就得到

所以 u满足问题

因u∈S,有u(a)=0,u满足(P)中的左端边界条件。

$$
\begin{array}{l} \left. - \int_ {a} ^ {b} v \frac {d}{d x} \left(p \frac {d u}{d x}\right) d x + \left(p \frac {d u}{d x} v\right) \right| _ {x = b} - \left. \left(p \frac {d u}{d x} v\right) \right| _ {x = a} \tag {2.5} \\ + \int_ {a} ^ {b} q u v d x - \int_ {a} ^ {b} f v d x = 0, \quad \forall v \in S _ {0} ^ {1}. \\ \int_ {a} ^ {b} (L u - f, v) = 0, \quad \forall v \in C _ {0} ^ {1} [ a, b ]. \\ \end{array}
$$

所以 满足问题 中的微分方程 Lu- f = 0,

把方程 代回（2.5）式，

得到 $( p { \frac { d u } { d x } } \nu ) \Bigg | _ { x = b } = 0 ,$ 10 . v S  

$$
\begin{array}{l} \begin{array}{l l} (P _ {1}) & \left\{ \begin{array}{l l} - \frac {d}{d x} \left(p (x) \frac {d u}{d x}\right) + q (x) u = f (x), & a <   x <   b, \\ u (a) = 0, & \frac {d u (b)}{d x} = 0. \end{array} \right. \end{array} \\ \left. \left(p \frac {d u}{d x} v\right) \right| _ {x = b} = 0, \quad \forall v \in S _ {0} ^ {1}. \\ \end{array}
$$

0p x p ( ) 0,   而且 是任意的，所以有 =0x=b

1u P 满足( ) 的右端边界条件

以上导出了u是问题（P)的解。

$$
\begin{array}{l l} \left(P _ {1}\right) & \left\{ \begin{array}{l l} - \frac {d}{d x} \left(p (x) \frac {d u}{d x}\right) + q (x) u = f (x), & a <   x <   b, \\ u (a) = 0, & \frac {d u (b)}{d x} = 0. \end{array} \right. \end{array}
$$

10 ,u S求 使 得2 ( )P( , ) ( ) 0,D u v F v − = 10 .v S 

2在( ) P 的提法中，边界上只要求

2在u P 满足( ) 的解，且

1则u P 满足( ) 其中的右端边界条件 =0x=kb

是自然满足的自然边界条件。是关于力的边界条件而第一类边界条件 我们称为约束边界条件。 

$$
D (u, v) = \int_ {a} ^ {b} \left(p \frac {d u}{d x} \frac {d v}{d x} + q u v\right) d x \quad F (v) = \int_ {a} ^ {b} f v d x
$$

（3）考虑

$$
\begin{array}{l} J (u) = \frac {1}{2} D (u, u) - F (u) \tag {2.6} \\ = \frac {1}{2} \int_ {a} ^ {b} \left[ p \left(\frac {d u}{d x}\right) ^ {2} + q u ^ {2} - 2 f u \right] d x \\ \end{array}
$$

代表一维问题在 的边界条件下的总势能。

从最小势能原理出发考虑下面问题

10 ,u S 求 使 得3 ( )P ( ) ( ),J u J w  10 . w S  

$$
J (u) = \frac {1}{2} D (u, u) - F (u) = \frac {1}{2} \int_ {a} ^ {b} \left[ p \left(\frac {d u}{d x}\right) ^ {2} + q u ^ {2} - 2 f u \right] d x
$$

10 ,u S 求 使 得3 ( )P ( ) ( ),J u J w  10 . w S   （2.6）

设u是问题 $( P _ { 3 } )$ 的解，则对任意的令w u v = +  , $w \in S _ { \mathbf { 0 } } ^ { 1 } ,$ 所 以 J u J u v ( ) ( )  +  。

$$
\begin{array}{l} J (u + \alpha v) = \frac {1}{2} D (u + \alpha v, u + \alpha v) - F (u + \alpha v) \\ = J (u) + \alpha \left[ D (u, v) - F (v) \right] + \frac {\alpha^ {2}}{2} D (v, v) \\ \end{array}
$$

10 ,u S 求 使 得3 ( )P ( ) ( ),J u J w  10 . w S  

$$
\begin{array}{l} J (u + \alpha v) = \frac {1}{2} D (u + \alpha v, u + \alpha v) - F (u + \alpha v) \\ = J (u) + \alpha \left[ D (u, v) - F (v) \right] + \frac {\alpha^ {2}}{2} D (v, v) ^ {(2. 7)} \\ \end{array}
$$

J u v ( ) +  可以看成 的一元函数

3u P 是( ) 的解 所以函数 达到极小，

有 ( )d J u v+ =0 D(u,v)-F(v)=0, ∀v ∈ Sld 0

$$
\begin{array}{l} \begin{array}{r l} {\left(P _ {3}\right)} & {\quad \left\{ \begin{array}{l l} {\text {求} u \in S _ {0} ^ {1}, \text {使 得}} \\ {J (u) \leq J (w), \quad \forall w \in S _ {0} ^ {1}.} \end{array} \right.} \end{array} \\ \begin{array}{r l} {\left(P _ {2}\right)} & {\quad \left\{ \begin{array}{l l} {\text {求} u \in S _ {0} ^ {1}, \text {使 得}} \\ {D (u, v) - F (v) = 0, \quad \forall v \in S _ {0} ^ {1}.} \end{array} \right.} \end{array} \\ D (u, v) - F (v) = 0, \quad \forall v \in S _ {0} ^ {1} \\ \end{array}
$$

结论是： 的解，则 的解。

（4）反过来，设 的解，则

10 ,u S 求 使得2 ( )P( , ) ( ) 0,D u v F v − = 10 .v S 

$$
\begin{array}{l} J (u + \alpha v) = \frac {1}{2} D (u + \alpha v, u + \alpha v) - F (u + \alpha v) \\ = J (u) + \alpha \left[ D (u, v) - F (v) \right] + \frac {\alpha^ {2}}{2} D (v, v) \tag {2.7} \\ \end{array}
$$

（4）反过来，设 的解，则

因（2.7）式对任意的 成立，所以有

$$
J (u + \alpha v) = J (u) + \frac {\alpha^ {2}}{2} D (v, v)
$$

一切 1 有D(v,v)≥0 所以J(u)≤J(u+av)

$$
\begin{array}{r l} {\left(P _ {2}\right)} & {\quad \left\{ \begin{array}{l l} {\text {求} u \in S _ {0} ^ {1}, \text {使 得}} \\ {D (u, v) - F (v) = 0, \quad \forall v \in S _ {0} ^ {1}.} \end{array} \right.} \end{array}
$$

10 ,u S 求 使得3 ( )P ( ) ( ),J u J w  10 . w S  

所以J(u)≤J(u+av)

对一切 成立，可以写成

J(u)≤J(w), Aw eS1. 导出了u是(P)的解。

$$
\begin{array}{l l} (P _ {1}) & \left\{ \begin{array}{l l} - \frac {d}{d x} \left(p (x) \frac {d u}{d x}\right) + q (x) u = f (x), & a <   x <   b, \\ u (a) = 0, & \frac {d u (b)}{d x} = 0. \end{array} \right. \end{array}
$$

u S  求 使 得( )P( , ) ( ) 0,D u v F v − = .v S 

10 ,u S  求 使 得3( )P ( ) ( ),J u J w  10  . w S  

（ 1 ） 是等价的。

（ 2 ） 若u是(P)的解，则u必是(P)的解。

（ 3 ） 若u是 $( P _ { 2 } )$ 的解，且则u必是(P)的解。 

$$
\begin{array}{l l} (P _ {1}) & \left\{ \begin{array}{l l} - \frac {d}{d x} \left(p (x) \frac {d u}{d x}\right) + q (x) u = f (x), & a <   x <   b, \\ u (a) = 0, & \frac {d u (b)}{d x} = 0. \end{array} \right. \end{array}
$$

10 ,u S 求 使 得2 ( )P( , ) ( ) 0,D u v F v − = 10 .v S 

10 ,u S 求 使 得3 ( )P ( ) ( ),J u J w  10 . w S  

$( P _ { 1 } )$ 是微分方程的两点边值问题，对应力学的平衡方程定解问题，

2( ) P 称为对应 的Galerkin 变分问题，

它对应力学的虚功原理。

3( ) P 称为Ritz变分问题，它对应最小势能原理。

# 2.2 非齐次约束边界条件的处理

第一类非齐次边界条件的定解问题，讨论对应的变分问题。设定解问题为

$$
\left\{ \begin{array}{l l} - \frac {d}{d x} \left(p (x) \frac {d u}{d x}\right) + q (x) u = f (x), & a <   x <   b, \\ u (a) = u _ {0}, \quad \frac {d u (b)}{d x} = 0. \end{array} \right.
$$

为了导出它对应的变分问题，方法之一是先使边界条件齐次化，令 ，则

$$
\left\{ \begin{array}{l l} - \frac {d}{d x} \left(p (x) \frac {d u}{d x}\right) + q (x) u = f (x), & a <   x <   b, \\ u (a) = u _ {0}, \quad \frac {d u (b)}{d x} = 0. \end{array} \right.
$$

令 则

$$
\left\{ \begin{array}{l} - \frac {d}{d x} \left(p (x) \frac {d \bar {u}}{d x}\right) + q (x) \bar {u} = f (x) - q (x) u _ {0}, \quad a <   x <   b, \\ u (a) = 0, \quad \frac {d \bar {u} (b)}{d x} = 0. \end{array} \right.
$$

$$
S = \left\{v \mid v \in C ^ {1} [ a, b ], \quad v (a) = u _ {0} \right\}
$$

[求u∈S,使得Galerkin变分问题是[D(u,v)-F(v)=0,

$$
\left\{ \begin{array}{l l} - \frac {d}{d x} \left(p (x) \frac {d \overline {{u}}}{d x}\right) + q (x) \overline {{u}} = f (x) - q (x) u _ {0}, & a <   x <   b, \\ u (a) = 0, \quad \frac {d \overline {{u}} (b)}{d x} = 0. \end{array} \right.
$$

$$
S = \left\{\nu \mid \nu \in C ^ {1} [ a, b ], \quad \nu (a) = u _ {0} \right\}
$$

[求u∈S,使得Galerkin变分问题是[D(u,v)-F(v)=0, Av eS.

,u S 求 使 得而Ritz变分问题是( ) ( ),J u J w  . w S  

10以上的S D u v F v J u ， ( , ) ( ) ( ) 和 的定义与上一 小节相同

# 2.3 第二、三类边界条件

$$
\left\{ \begin{array}{l} - \frac {d}{d x} \left(p (x) \frac {d u}{d x}\right) + q (x) u = f (x), \quad a <   x <   b, \\ u (a) = 0, \quad p (b) \frac {d u (b)}{d x} + \alpha u (b) = g. \end{array} \right.
$$

0其中p x p q x ( ) 0, ( ) 0    ，常数  0,

在x b = 的边界条件包括了第二类和第三类边界条件完全可以推导出对应的变分问题，

$$
D (u, v) = \int_ {a} ^ {b} \left[ p \frac {d u}{d x} \frac {d v}{d x} + q u v \right] d x + \alpha u (b) v (b)
$$

$$
F (v) = \int_ {a} ^ {b} f v d x + g v (b)
$$