# 第五章 椭圆型方程的差分方法

# 主要内容

■ Poisson方程

差分方程的性质

边界条件的处理

变系数方程

双调和方程

特征值问题

# 第一节 Poisson 方程

考虑方程

$$
\Delta u = \frac {\partial^ {2} u}{\partial x ^ {2}} + \frac {\partial^ {2} u}{\partial y ^ {2}} = f (x, y), (x, y) \in D \tag {1}
$$

其中 $D$ 是平面内有界区域, 边界用 $\partial D$ 表示。

不妨设 $\pmb{D}$ 为矩形区域, $D = \{(x, y) | 0 < x < 1, 0 < y < 1\}$ , 方程满足边值条件 (Dirichlet条件)

$$
u (x, y) = a (x, y) \quad (x, y) \in \partial \Omega
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/c763622a-db61-41ee-8d7f-a09a4a0bdbf5/12f560103c8120f483193eb256d0bcac24b45e0ec3e9c138a1cdff86f4bb5c4d.jpg)


区域D的内点

$$
D _ {h} = \{(x _ {i}, y _ {j}), x _ {i} = i h _ {1}, 1 \leq i \leq \mathrm {I}; y _ {j} = j h _ {2}, 1 \leq j \leq \mathrm {J} \}
$$

区域D的边界点

$$
\partial D = \{(x _ {i}, y _ {j}), x _ {i} = i h _ {1}, y _ {j} = j h _ {2}, \mathrm {i = 0 , I + 1 , 0 \leq j \leq J + 1 ;}
$$

$$
j = 0, \quad J + 1, 0 \leq i \leq I + 1 \}
$$

# 1.1:五点差分格式

现在假定 $\left(x_{i}, y_{j}\right)$ 为内点，沿 $x, y$ 方向分别用二阶中心差商代替 $u_{xx}, u_{yy}$ ，则得差分方程

$$
\Delta_ {h} u _ {i j} = - \left[ \frac {u _ {i + 1 , j} - 2 u _ {i j} + u _ {i - 1 , j}}{h _ {1} ^ {2}} + \frac {u _ {i , j + 1} - 2 u _ {i j} + u _ {i , j - 1}}{h _ {2} ^ {2}} \right] = f _ {i j} \tag {1.1}
$$

差分格式为：

$$
\left\{ \begin{array}{l l} \Delta_ {h} u _ {i j} = f _ {i j}, & (x _ {i}, y _ {j}) \in D _ {h} \\ u _ {i j} = a _ {i j}, & (x _ {i}, y _ {j}) \in \partial D _ {h} \end{array} \right.
$$

利用Taylor展式

$$
\begin{array}{l} \frac {u (x _ {i + 1} , y _ {j}) - 2 u (x _ {i} , y _ {j}) + u (x _ {i - 1} , y _ {j})}{h _ {1} ^ {2}} \\ = \frac {\partial^ {2} u (x _ {i} , y _ {j})}{\partial x ^ {2}} + \frac {h _ {1} ^ {2}}{1 2} \frac {\partial^ {4} u (x _ {i} , y _ {j})}{\partial x ^ {4}} + \frac {h _ {1} ^ {4}}{3 6 0} \frac {\partial^ {6} u (x _ {i} , y _ {j})}{\partial x ^ {6}} + O \left(h _ {1} ^ {6}\right), \tag {1.2} \\ \end{array}
$$

$$
\begin{array}{l} \frac {u (x _ {i} , y _ {j + 1}) - 2 u (x _ {i} , y _ {j}) + u (x _ {i} , y _ {j - 1})}{h _ {1} ^ {2}} \\ = \frac {\partial^ {2} u \left(x _ {i} , y _ {j}\right)}{\partial y ^ {2}} + \frac {h _ {2} ^ {2}}{1 2} \frac {\partial^ {4} u \left(x _ {i} , y _ {j}\right)}{\partial y ^ {4}} + \frac {h _ {2} ^ {4}}{3 6 0} \frac {\partial^ {6} u \left(x _ {i} , y _ {j}\right)}{\partial y ^ {6}} + O \left(h _ {2} ^ {6}\right), (1. 3) \\ \end{array}
$$

可得差分算子 $\Delta_{h}$ 的截断误差

$$
R _ {i j} (u) = \Delta u (x _ {i}, y _ {j}) - \Delta_ {h} u (x _ {i}, y _ {j})
$$

$$
\begin{array}{l} = - \frac {1}{1 2} [ h _ {1} ^ {2} \frac {\partial^ {4} u (x _ {i} , y _ {j})}{\partial x ^ {4}} + h _ {2} ^ {2} \frac {\partial^ {4} u (x _ {i} , y _ {j})}{\partial y ^ {4}} ] + O (h ^ {4}), \\ = O (h _ {1} ^ {2} + h _ {2} ^ {2}) = O (h ^ {2}), \\ \end{array}
$$

由于差分方程(1.1)中只出现 $u$ 在 $(i, j)$ 及其四个邻点上的值，故称为五点差分格式

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/c763622a-db61-41ee-8d7f-a09a4a0bdbf5/d4ccc17fc4461b6bc596059c0e18176e8a028f3c020394927ed3c1ca347f77cf.jpg)


特别取正方形网格: $h_{1} = h_{2} = h$ , 且设 $a = b$ 则差分方程(1.1)简化为

$$
4 u _ {i j} - (u _ {i + 1, j} + u _ {i - 1, j} + u _ {i, j + 1} + u _ {i, j - 1}) = h ^ {2} f _ {i j}.
$$

若 $f \equiv 0(\text{Laplace方程})$ , 则有

$$
4 u _ {i j} = (u _ {i + 1, j} + u _ {i - 1, j} + u _ {i, j + 1} + u _ {i, j - 1})
$$

单位正方形中的内部结点上的 $J^{2}$ 个线性方程 (1.2) 写成矩阵形式为

$$
A U = K
$$

其中， $A$ 是 $J^2$ 阶方阵

$$
A = \left[ \begin{array}{c c c c c} B & - I & & & \\ - I & B & - I & & \\ & \ddots & \ddots & \ddots & \\ & & - I & B & - I \\ & & & - I & B \end{array} \right]
$$

$I$ 是 $J$ 阶单位方阵; $B$ 是 $J$ 阶方阵。

$$
B = \left[ \begin{array}{c c c c c} 4 & - 1 & & & \\ - 1 & 4 & - 1 & & \\ & \ddots & \ddots & \ddots & \\ & & - 1 & 4 & - 1 \\ & & & - 1 & 4 \end{array} \right]
$$

# 1.2: 九点差分格式

注 若将(1.2)，(1.3)两式相加，则得

$$
\begin{array}{l} \Delta_ {h} u (x _ {i}, y _ {j}) = \Delta u (x _ {i}, y _ {j}) \\ + \frac {1}{1 2} [ h _ {1} ^ {2} \frac {\partial^ {4} u (x _ {i} , y _ {j})}{\partial x ^ {4}} + h _ {2} ^ {2} \frac {\partial^ {4} u (x _ {i} , y _ {j})}{\partial y ^ {4}} ] + O (h ^ {4}) \\ \end{array}
$$

$$
= \Delta u (x _ {i}, y _ {j})
$$

$$
+ \frac {1}{1 2} \left(h _ {1} ^ {2} \frac {\partial^ {2}}{\partial x ^ {2}} + h _ {2} ^ {2} \frac {\partial^ {2}}{\partial y ^ {2}}\right) \left(\frac {\partial^ {2} u (x _ {i} , y _ {j})}{\partial x ^ {2}} + \frac {\partial^ {2} u (x _ {i} , y _ {j})}{\partial y ^ {2}}\right)
$$

$$
\Delta u (x _ {i}, y _ {j}) = f (x _ {i}, y _ {j})
$$

$$
- \frac {h _ {1} ^ {2} + h _ {2} ^ {2}}{1 2} \frac {\partial^ {4} u (x _ {i} , y _ {j})}{\partial x ^ {2} \partial y ^ {2}} + O (h ^ {4})
$$

又

$$
\frac {\partial^ {4} u (x _ {i} , y _ {j})}{\partial x ^ {2} \partial y ^ {2}} = \frac {u _ {x x} ^ {\prime \prime} (x _ {i} , y _ {j + 1}) - 2 u _ {x x} ^ {\prime \prime} (x _ {i} , y _ {j}) + u _ {x x} ^ {\prime \prime} (x _ {i} , y _ {j - 1})}{h _ {2} ^ {2}} + 0 \left(h _ {2} ^ {2}\right)
$$

$$
\begin{array}{l} = \frac {1}{h _ {1} ^ {2} h _ {2} ^ {2}} \left[ u \left(x _ {i + 1}, y _ {j + 1}\right) - 2 u \left(x _ {i}, y _ {j + 1}\right) + u \left(x _ {i - 1}, y _ {j + 1}\right) \right. \\ - 2 (u (x _ {i + 1}, y _ {j}) - 2 u (x _ {i}, y _ {j}) + u (x _ {i - 1}, y _ {j})) + u (x _ {i + 1}, y _ {j - 1}) \\ - 2 u (x _ {i}, y _ {j - 1}) + u (x _ {i - 1}, y _ {j - 1}) ] + 0 (h ^ {2}). \\ \end{array}
$$

因此

$$
\begin{array}{l} \Delta_ {h} u (x _ {i}, y _ {j}) + \frac {1}{1 2} \frac {h _ {1} ^ {2} + h _ {2} ^ {2}}{h _ {1} ^ {2} h _ {2} ^ {2}} [ 4 u (x _ {i}, y _ {j}) - 2 u (x _ {i - 1}, y _ {j}) + u (x _ {i}, y _ {j - 1}) \\ + u (x _ {i + 1}, y _ {j}) + u (x _ {i}, y _ {j + 1}) + u (x _ {i - 1}, y _ {j - 1}) + u (x _ {i + 1}, y _ {j - 1}) \\ + u (x _ {i + 1}, y _ {j + 1}) + u (x _ {i - 1}, y _ {j + 1}) ] \\ = - f (x _ {i}, y _ {j}) - \frac {1}{1 2} \left(h _ {1} ^ {2} \frac {\partial^ {2} f (x _ {i} , y _ {j})}{\partial x ^ {2}} + h _ {2} ^ {2} \frac {\partial^ {2} f (x _ {i} , y _ {j})}{\partial y ^ {2}}\right) + 0 \left(h ^ {4}\right) \\ \end{array}
$$

舍去截断误差项，便得到逼近方程的九点差分格式：

$$
\begin{array}{l} - \Delta_ {h} u _ {i j} - \frac {1}{1 2} \frac {h _ {1} ^ {2} + h _ {2} ^ {2}}{h _ {1} ^ {2} h _ {2} ^ {2}} [ 4 u _ {i j} - 2 (u _ {i - 1, j} + u _ {i, j - 1} + u _ {i + 1, j} + u _ {i, j + 1}) \\ + u _ {i - 1, j - 1} + u _ {i + 1, j - 1} + u _ {i + 1, j + 1} + u _ {i - 1, j + 1}) ] \\ = f _ {i j} + \frac {1}{1 2} \left(h _ {1} ^ {2} \frac {\partial^ {2} f (x _ {i} , y _ {j})}{\partial x ^ {2}} + h _ {2} ^ {2} \frac {\partial^ {2} f (x _ {i} , y _ {j})}{\partial y ^ {2}}\right) \\ \end{array}
$$

其截断误差的阶为 $0(h^{4})$ .

附：九点差分格式推导 方法二

取 $\xi = h \frac{\partial}{\partial x}, \eta = h \frac{\partial}{\partial y}, D = \frac{\partial^2}{\partial x \partial y}$ 

$$
\xi^ {2} + \eta^ {2} = h ^ {2} \Delta , \xi \eta = h ^ {2} D, \xi^ {4} + \eta^ {4} = h ^ {4} (\Delta^ {2} - 2 D ^ {4})
$$

由Taylor公式有：

$$
u (x + h) = (1 + h \frac {d}{d x} + \frac {h ^ {2}}{2} \frac {d ^ {2}}{d x ^ {2}} + \dots + \frac {h ^ {n}}{n !} \frac {d ^ {n}}{d x ^ {n}}) u (x) = e ^ {h \frac {d}{d x}} u (x)
$$

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/c763622a-db61-41ee-8d7f-a09a4a0bdbf5/42a2892245a6651eaf7bbb133dfe3eb492d76c62e5fd70ae1dbf26dca3c0e15b.jpg)


令 $S_{1} = u_{1} + u_{2} + u_{3} + u_{4}, S_{2} = u_{5} + u_{6} + u_{7} + u_{8}$ 

$$
\begin{array}{l} S _ {1} = u _ {1} + u _ {2} + u _ {3} + u _ {4} \\ = (1 + \xi + \frac {\xi^ {2}}{2} + \frac {\xi^ {3}}{3 !} + \frac {\xi^ {4}}{4 !} + 1 + \eta + \frac {\eta^ {2}}{2} + \frac {\eta^ {3}}{3 !} + \frac {\eta^ {4}}{4 !} \\ + 1 - \xi + \frac {\xi^ {2}}{2} - \frac {\xi^ {3}}{3 !} + \frac {\xi^ {4}}{4 !} + 1 - \eta + \frac {\eta^ {2}}{2} - \frac {\eta^ {3}}{3 !} + \frac {\eta^ {4}}{4 !}) \\ = 4 u _ {0} + h ^ {2} \Delta u _ {0} + \frac {1}{1 2} h ^ {4} (\Delta^ {2} - 2 D ^ {4}) u _ {0} + O (h ^ {6}) \\ \end{array}
$$

$$
S _ {2} = 4 u _ {0} + 2 h ^ {2} \Delta u _ {0} + \frac {1}{6} h ^ {4} (\Delta^ {2} + 4 D ^ {4}) u _ {0} + O (h ^ {6})
$$

从上面2个式子消除 $D^{4} u_{0}$ , 并忽略高阶项, 得

$$
4 S _ {1} + S _ {2} - 2 0 u _ {0} = 6 h ^ {2} f _ {i j} + \frac {1}{2} h ^ {4} \Delta f _ {i j}
$$

积分法

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/c763622a-db61-41ee-8d7f-a09a4a0bdbf5/ca1d0bb6344df58c5c56cb30472c75b62945f8af5ebf29b84dd922863b80c450.jpg)


$$
\int_ {D} \Delta u d x d y = \int_ {D} f d x d y
$$

由 Gauss公式有

$$
\int_ {D} \Delta u d x d y = \int_ {\partial D} \frac {\partial u}{\partial n} d s = (\int_ {l _ {1}} + \int_ {l _ {2}} + \int_ {l _ {3}} + \int_ {l _ {4}}) \frac {\partial u}{\partial n} d s
$$

从上面的公式把边界积分转化为在四条边上的积分,于是积分分成4段进行, 也即转化为定积分。

对于在边 $l_{1}$ 上, 因为在矩形边 $l_{1}$ 上的外法向就是 $y$ 的负方向, 弧长的微分 $d s = d x$ , 于是对此定积分用矩形公式近似计算, 并且用差商代替微商, 得到

$$
\begin{array}{l} \int_ {l _ {1}} \frac {\partial u}{\partial n} d s = - \int_ {l _ {1}} \frac {\partial u}{\partial y} d x = - \frac {\partial u}{\partial y} (x _ {i}, y _ {j - \frac {1}{2}}) h + O (h ^ {3}) \\ = - \frac {u _ {i , j} - u _ {i , j - 1}}{h} h + O \left(h ^ {3}\right) \\ \end{array}
$$

对于在边 $l_{2}$ 上, 因为在矩形边 $l_{2}$ 上的外法向就是 $x$ 的正方向,弧长的微分 $d s = d y$ , 于是对此定积分用矩形公式近似计算, 并且用差商代替微商, 得到

$$
\begin{array}{l} \int_ {l _ {2}} \frac {\partial u}{\partial n} d s = \int_ {l _ {2}} \frac {\partial u}{\partial x} d y = \frac {\partial u}{\partial x} (x _ {i + \frac {1}{2}}, y _ {j}) h + O (h ^ {3}) \\ = \frac {u _ {i + 1 , j} - u _ {i , j}}{h} h + O (h ^ {3}) \\ \end{array}
$$

类似的计算出另外两条边上的积分如下:

$$
\int_ {l _ {3}} \frac {\partial u}{\partial n} d s = - \int_ {l _ {3}} \frac {\partial u}{\partial y} d x = - \frac {\partial u}{\partial y} \left(x _ {i}, y _ {j + \frac {1}{2}}\right) h + O \left(h ^ {3}\right)
$$

$$
= - \frac {u _ {i , j + 1} - u _ {i , j}}{h} h + O (h ^ {3})
$$

$$
\begin{array}{l} \int_ {l _ {4}} \frac {\partial u}{\partial n} d s = \int_ {l _ {4}} \frac {\partial u}{\partial x} d y = \frac {\partial u}{\partial x} \left(x _ {i - \frac {1}{2}}, y _ {j}\right) h + O \left(h ^ {3}\right) \\ = \frac {u _ {i , j} - u _ {i - 1 , j}}{h} h + O \left(h ^ {3}\right) \\ \end{array}
$$

# 于是可以得到差分格式，即5点格式

$$
u _ {i + 1, j} + u _ {i - 1, j} + u _ {i, j + 1} + u _ {i, j - 1} - 4 u _ {i j} = h ^ {2} f _ {i j}
$$

注：(1) 关于积分区域的选取，选择积分区域的方法不同，得到的格式也有所不同，在此方法中积分区域是围绕节点的矩形。也可以直接用差分的网格作为积分区域，或者在非矩形网格上（例如三角形网格）分割积分区域。如下图:

![image](https://cdn-mineru.openxlab.org.cn/result/2026-04-29/c763622a-db61-41ee-8d7f-a09a4a0bdbf5/8321ebab10b339435ceb0117b8c9078d1f8381e7bc3074d266f03af4de6f3313.jpg)


(2)采用不同的数值积分，比如梯形公式等等也可以得到不同的差分格式。

# 1.3:极坐标下的差分格式

如果求解域是圆环、环形域或扇行域,

则采用极坐标是方便的，此时Poisson方程形如

$$
- \Delta_ {r, \theta} u = - [ \frac {1}{r} \frac {\partial}{\partial r} (r \frac {\partial u}{\partial r}) + \frac {1}{r ^ {2}} \frac {\partial^ {2} u}{\partial \theta^ {2}} ] = f (r, \theta), (2)
$$

方程(2)的系数于 $r = 0$ 奇异, 因此只当 $r > 0$ 时有意义。为了定出有兴趣的解, 需补充 $u$ 于 $r = 0$ 有界的条件。从而 $u$ 满足

$$
\lim  _ {r \rightarrow 0 ^ {+}} r \frac {\partial u}{\partial r} = 0. \tag {1.8}
$$

关于变量 $r, \theta$ 分别取等步长 $h_r$ 和 $h_\theta$ . 令

$$
\begin{array}{l} r _ {i} = (i + 0. 5) h _ {r}, i = 0, 1, 2, \dots , \\ \theta_ {j} = (j + 1) h _ {\theta}, j = 0, 1, \dots , J - 1, h _ {\theta} = 2 \pi / J \\ \end{array}
$$

对任一点 $(r_{i}, \theta_{j})(i > 0)$ , 用中心差商公式

$$
\begin{array}{l} \left[ \frac {1}{r} \frac {\partial}{\partial r} \left(r \frac {\partial u}{\partial r}\right) \right] _ {\left(r _ {i}, \theta_ {j}\right)} \approx \frac {1}{r _ {i}} \frac {r _ {i + \frac {1}{2}} u _ {i + 1 , j} - \left(r _ {i + \frac {1}{2}} + r _ {i - \frac {1}{2}}\right) u _ {i j} + r _ {i - \frac {1}{2}} u _ {i - 1 , j}}{h _ {r} ^ {2}} \\ \left[ \frac {1}{r ^ {2}} \frac {\partial^ {2} u}{\partial \theta^ {2}} \right] _ {(r _ {i}, \theta_ {j})} \approx \frac {1}{r _ {i} ^ {2}} \frac {u _ {i + 1 , j} - 2 u _ {i , j} + u _ {i , j - 1}}{h _ {\theta} ^ {2}} \\ \left(u _ {i 0} = u _ {i J}, u _ {i, - 1} = u _ {i J - 1}\right) \\ \end{array}
$$

代到(2),则就 $i \geq 1$ 得到逼近它的差分方程

$$
\begin{array}{l} - \left[ \frac {1}{r _ {i}} \frac {r _ {i + \frac {1}{2}} u _ {i + 1 , j} - (r _ {i + \frac {1}{2}} + r _ {i - \frac {1}{2}}) u _ {i j} + r _ {i - \frac {1}{2}} u _ {i - 1 , j}}{h _ {r} ^ {2}} \right. \\ + \frac {1}{r _ {i} ^ {2}} \frac {u _ {i + 1 , j} - 2 u _ {i , j} + u _ {i , j - 1}}{h _ {\theta} ^ {2}} ] = f (r _ {i}, \theta_ {j}). \\ \end{array}
$$

今用积分插值法导出点 $\left(r_{0}, \theta_{j}\right)$ 的差分方程, 以 $r$ 乘方程(1), 并对 $r$ 由 $\varepsilon$ 到 $h_{r}$ 积分, 对 $\theta$ 由 $\theta_{j - \frac{1}{2}}$ 到 $\theta_{j + \frac{1}{2}}$ 积分, 然后令 $\varepsilon \rightarrow 0$ , 且注意条件(1.8), 则得

$$
\begin{array}{l} - [ h _ {r} \int_ {\theta_ {j - \frac {1}{2}}} ^ {\theta_ {j + \frac {1}{2}}} \frac {\partial}{\partial r} u (h _ {r}, \theta) d \theta + \\ \int_ {0} ^ {h _ {r}} \frac {\partial}{\partial \theta} u (r, \theta_ {j + \frac {1}{2}}) - \frac {\partial}{\partial \theta} u (r, \theta_ {j - \frac {1}{2}}) ] d r ] \\ = \int_ {0} ^ {h _ {r}} r d r \int_ {\theta_ {j - \frac {1}{2}}} ^ {\theta_ {j + \frac {1}{2}}} f (r, \theta) d \theta . \\ \end{array}
$$

用中矩形公式代替上述积分，则

$$
\begin{array}{l} - \left[ h _ {r} h _ {\theta} \frac {\partial}{\partial r} u \left(h _ {r}, \theta_ {j}\right) + 2 \left(\frac {\partial}{\partial \theta} u \left(\frac {h _ {r}}{2}, \theta_ {j + \frac {1}{2}}\right) - \frac {\partial}{\partial \theta} u \left(\frac {h _ {r}}{2}, \theta_ {j - \frac {1}{2}}\right) \right] \right. \\ \approx \frac {1}{2} h ^ {2} r h _ {\theta} f (\frac {h _ {r}}{2}, \theta_ {j}). \\ \end{array}
$$

再用中心差商代替微商, 就得出点 $\left(r_{0}, \theta_{j}\right)$ 的差分方程:

$$
- \left[ h _ {r} h _ {\theta} \frac {u _ {1 j} - u _ {0 j}}{h _ {r}} + 2 \frac {u _ {0 , j + 1} - 2 u _ {0 j} + u _ {0 , j - 1}}{h _ {\theta}} \right] = \frac {1}{2} h ^ {2} r h _ {\theta} f _ {0 j}.
$$

或(两端除以 $\frac{1}{2} h_{r}^{2} h_{\theta}$ )

$$
- [ \frac {2}{h _ {r}} \frac {u _ {1 j} - u _ {0 j}}{h _ {r}} + \frac {4}{h _ {r} ^ {2}} \frac {u _ {0 , j + 1} - 2 u _ {0 j} + u _ {0 , j - 1}}{h _ {\theta} ^ {2}} ] = f _ {0 j}.
$$