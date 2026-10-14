# 第一阶段：模型推导与最小求解器

版本 V0.1 | 2026-10-02 | JIajun FU

本阶段完成一维非线性浅水模型的推导、线性化、无量纲化和一阶有限体积求解器。验证对象是已知数学解与守恒性质；尚未使用 NOAA 实验验证，也未模拟真实海啸事件。

## 1. 几何、符号与模型范围

考虑固定宽度、单位横向宽度的一维水域。x 是水平坐标，z 是向上为正的垂直坐标；底床为 z=b(x)，自由表面为 z=eta(x,t)。底床静止，平均静水位定为 0。

| 符号 | 含义 | 单位 |
|---|---|---|
| b | 底床高程，水下通常为负 | m |
| eta | 自由表面高程，不是总水深 | m |
| h=eta-b | 总水深 | m |
| D=-b | 静水深，要求 D>0 | m |
| u | 深度平均水平速度 | m/s |
| q=hu | 单位宽度体积流量 | m²/s |
| g | 重力加速度，默认 9.81 | m/s² |
| rho | 常数水密度；除以 rho 后不出现于方程 | kg/m³ |

静水时 eta=0，h=D=-b，u=q=0。波幅 A 与深度 h 是不同变量。例如 D=50 m、eta=1 m 时，总水深 h=51 m。

假设：不可压缩、密度恒定；波的水平尺度远大于水深；垂向加速度较小，采用静水压力；水平速度的垂向变化可忽略；底床固定；忽略摩擦、旋转、风应力和色散。一维化意味着不描述横向折射、绕射或海湾二维聚焦。

“浅水”指水深相对于波长或波尺度很小，不要求水深的绝对数值很小。海啸长波可以在深海中满足这一比例条件。

本代码要求每个物理单元始终有水。一旦接近干涸，就明确报错；上岸爬高与干湿边界属于后续阶段。

## 2. 质量守恒：从不可压缩到一维方程

二维不可压缩连续方程为

\[
\partial_x v+\partial_z w=0,
\]

其中 v(x,z,t)、w(x,z,t) 是水平和垂直速度。定义深度平均速度

\[
u(x,t)=\frac1h\int_{b(x)}^{\eta(x,t)}v(x,z,t)\,dz.
\]

用 Leibniz 法则积分连续方程：

\[
\partial_x\int_b^{\eta}v\,dz
-v_s\eta_x+v_b b_x+w_s-w_b=0.
\]

下标 s、b 表示表面和底床值。自由表面与静止不可穿透底床的运动学条件分别为

\[
w_s=\eta_t+v_s\eta_x,\qquad w_b=v_b b_x.
\]

代入后得到 eta_t+partial_x(hu)=0；又因为 b 不随时间变化，eta_t=h_t，因此

\[
\boxed{h_t+(hu)_x=0},\qquad \boxed{h_t+q_x=0}.
\]

等价的控制体积解释：区间 [x_L,x_R] 内单位宽度水体积变化率等于左侧流入减右侧流出。

\[
\frac{d}{dt}\int_{x_L}^{x_R}h\,dx=q(x_L,t)-q(x_R,t).
\]

所以开放边界下总质量可以变化。应该检查“质量变化+边界流出积分=0”，而不是要求总质量不变。

## 3. 动量方程：静水压力与底床源项

垂直压力平衡给出

\[
p_z=-\rho g,\qquad p(\eta)=0,
\]

其中采用相对于大气压力的表压。故 p=rho g(eta-z)。压力沿水深积分为

\[
\int_b^{\eta}p\,dz=\frac12\rho gh^2.
\]

在深度均匀速度假设下，水平动量通量除以密度为 hu²；积分压力通量为 gh²/2。底床斜面产生水平压力力，单位水平长度、除以密度后为 -gh b_x。因此保守动量方程为

\[
\boxed{(hu)_t+\left(hu^2+\frac12gh^2\right)_x=-ghb_x}.
\]

令 q=hu，得到求解器采用的形式：

\[
\boxed{q_t+\left(\frac{q^2}{h}+\frac12gh^2\right)_x=-ghb_x}.
\]

也可在光滑、h>0 时结合质量方程得到原始变量形式

\[
u_t+u u_x+g(h+b)_x=0
\quad\Longleftrightarrow\quad
u_t+u u_x+g\eta_x=0.
\]

注意：原始变量的代数推导仅对光滑解成立。出现激波时，需要保守弱形式和适当熵条件；有限体积通量正是为此设计。

底床源项不能省略：静水时 h=-b，压力梯度 g h h_x 与 -gh b_x 恰好平衡。若计算时破坏该平衡，会在本来静止的水面上产生虚假波。

## 4. 双曲结构与传播速度

写成

\[
U_t+F(U)_x=S(U,b),\quad
U=\begin{pmatrix}h\\q\end{pmatrix},\quad
F=\begin{pmatrix}q\\q^2/h+gh^2/2\end{pmatrix},\quad
S=\begin{pmatrix}0\\-ghb_x\end{pmatrix}.
\]

对 h>0，通量 Jacobian 为

\[
F'(U)=\begin{pmatrix}0&1\\gh-u^2&2u\end{pmatrix},
\]

两个特征速度为

\[
\lambda_-=u-\sqrt{gh},\qquad\lambda_+=u+\sqrt{gh}.
\]

因此时间步由流速与重力波速共同决定。在静水深 D 上，小振幅波的基准速度为 c0=sqrt(gD)。它不是任意有限振幅波的精确速度。

## 5. 线性化：后续比较模型的基础

令 h=D(x)+eta，其中 D=-b。围绕 eta=u=0 展开，并忽略 eta*u、u² 等二阶小量：

\[
\boxed{\eta_t+(Du)_x=0},\qquad
\boxed{u_t+g\eta_x=0}.
\]

若使用线性流量 q_lin=Du，则

\[
\eta_t+(q_{lin})_x=0,\qquad(q_{lin})_t+gD\eta_x=0.
\]

q_lin 与非线性 q=(D+eta)u 只有在一阶近似下相等。比较两种模型时，要明确初始化匹配方式。

恒定水深 D 下，消去 u：

\[
\eta_{tt}-gD\eta_{xx}=0.
\]

右行解为 eta(x,t)=f(x-c0 t)，速度 u=(c0/D)eta，线性流量 q_lin=c0 eta。这用于小振幅波速测试。变水深时，(Du)_x 不能漏掉 D_x u。

**本次仅推导线性模型，没有实现独立线性求解器。** 当前波速测试是非线性求解器在小振幅极限下与线性理论比较。

## 6. 无量纲化与适用条件

选定参考深度 D0、水平尺度 L 和 c0=sqrt(gD0)，令

\[
X=x/L,\quad T=t c0/L,\quad H=h/D0,\quad
V=u/c0,\quad Q=q/(D0c0),\quad B=b/D0.
\]

得到

\[
H_T+Q_X=0,\qquad
Q_T+\left(Q^2/H+H^2/2\right)_X=-H B_X.
\]

两个重要参数为 epsilon=A/D0 和 delta=D0/L。epsilon 小支持线性化；delta 小支持浅水、近似无色散的长波假设。模型由有量纲方程无量纲化后不显式包含 delta，是因为建模过程中已经忽略了相关高阶项，不能据此认为 delta 不重要。

局部 eta/D(x) 也应监测：离岸小波幅进入很浅区域后，线性化仍可能失效。两个 SWE 模型互相接近也不证明色散效应可以忽略。

## 7. 有限体积离散

均匀网格 I_i=[x_(i-1/2),x_(i+1/2)]，单元宽 dx。未知量是**单元平均值**，而不是点值：

\[
U_i(t)=\frac1{\Delta x}\int_{I_i}U(x,t)\,dx.
\]

初始波使用 8 点 Gauss–Legendre 积分求单元平均。底床由每单元一个离散高程表示；示例中的光滑底床按中心采样，作为分片常数近似。

平底时，显式 Euler 更新为

\[
U_i^{n+1}=U_i^n-\frac{\Delta t}{\Delta x}
\left(\widehat F_{i+1/2}-\widehat F_{i-1/2}\right).
\]

空间与时间都是一阶。较粗网格会平滑波峰，出现数值扩散；这不是人为加入的真实摩擦。

### 7.1 HLL 数值通量

给定界面左右状态 U_L、U_R，选波速界 s_L<=0<=s_R：

\[
\widehat F=
\frac{s_RF(U_L)-s_LF(U_R)+s_Ls_R(U_R-U_L)}{s_R-s_L}.
\]

把 0 包括在波速界中，使纯右行或纯左行状态自然退化到相应的迎风物理通量。两侧全为零的退化界面返回零通量。

代码波速界包含左右 u±sqrt(gh) 以及 Roe 平均速度：

\[
\widetilde u=\frac{\sqrt{h_L}u_L+\sqrt{h_R}u_R}{\sqrt{h_L}+\sqrt{h_R}},
\qquad\widetilde c=\sqrt{g(h_L+h_R)/2}.
\]

s_L 为 0、u_L-c_L、u_R-c_R、u_tilde-c_tilde 的最小值；s_R 为相应正向量的最大值。

### 7.2 静水重构与源项修正

采用 Audusse 等（2004）的静水重构思想。界面底床取

\[
b^*=\max(b_L,b_R),\quad
h_L^*=\max(0,h_L+b_L-b^*),\quad
h_R^*=\max(0,h_R+b_R-b^*).
\]

保持原单元速度，构造 U_L^*=(h_L^*,h_L^*u_L) 和 U_R^*=(h_R^*,h_R^*u_R)，计算基础 HLL 通量 F*。

同一界面针对左右单元分别修正动量通量：

\[
F^L=F^*+\begin{pmatrix}0\\g(h_L^2-(h_L^*)^2)/2\end{pmatrix},
\qquad
F^R=F^*+\begin{pmatrix}0\\g(h_R^2-(h_R^*)^2)/2\end{pmatrix}.
\]

更新为

\[
\boxed{U_i^{n+1}=U_i^n-\frac{\Delta t}{\Delta x}
\left(F^L_{i+1/2}-F^R_{i-1/2}\right)}.
\]

源项已通过这些修正纳入；**不要再额外叠加中心差分 -gh b_x**，否则重复计算。

静水平衡解释：eta 常数、u=0 时，左右重构水深相同；基础质量通量为 0，基础动量通量为 g(h*)²/2。对单元 i，两侧修正后的动量通量都变为 g h_i²/2，故差为 0。浮点运算下可能有舍入级漂移。

左右界面质量通量相同，因而总质量可通过界面求和抵消；动量则受到真实底床力影响，不能在变底床时要求总动量恒定。

界面 h* 的非负截断属于重构算法，不能与“将演化后的负水深直接改为正数”混淆。本版本不做后者，也不声称实现通用干湿边界。

### 7.3 时间步与输出

\[
\Delta t=\mathrm{CFL}\frac{\Delta x}{a_{max}},\qquad\mathrm{CFL}=0.45,
\]

其中 a_max 同时覆盖单元的 |u|+sqrt(gh) 和实际界面波速界。代码允许 0<CFL<=0.5；这是一种保守实施选择，不是未经证明的任意输入稳定性保证。

每步检查有限值和最低水深。不满足完整有水假设则报错，不裁剪。为精确落在指定输出时刻，最后一小步会缩短。

## 8. 初始与边界条件

### 8.1 避免初始扰动一分为二

如果只给水面扰动但设 u=0，在线性平底模型中通常会分裂成左右两支波。右行小振幅初始化需配对 u=(c0/D)eta。

平底非线性右行简单波可保持一个 Riemann 不变量：

\[
u-2\sqrt{gh}=-2\sqrt{gD},\quad
u=2(\sqrt{gh}-\sqrt{gD}).
\]

这保证初始速度与波面符合右行简单波关系。在未形成激波前，波形会因非线性传播速度不同而变化；Gaussian 脉冲不是浅水方程的精确、永久保形孤立波。

### 8.2 三种边界

- periodic：两端循环连接，用于周期行波及收敛实验。
- reflective：幽灵单元水深和底床复制、q 反号，表示固壁；边界质量通量为零。
- outflow：复制端点状态与底床，即零梯度幽灵单元。它不是精确吸收边界，后续还需量化反射。本次平底传播演示让脉冲远离边界。

## 9. 独立验证参考

平底非线性简单波中，沿右行特征线 h 保持初值 h0(xi)，位置满足

\[
x=\xi+[3\sqrt{g h_0(\xi)}-2\sqrt{gD}]t.
\]

给定 x、t 用 Newton 解 xi，得到 h(x,t)=h0(xi)。这是独立于有限体积代码的参考，不是“拿同一个算法跑细网格来证明自身正确”。

本例 h0=D+A sin(2pi xi/L)，D=10 m、A=0.2 m、L=1000 m、t=20 s，未发生特征交叉。参考解和数值初值都按单元平均比较。程序检查一个保守的非交叉充分条件，并在接近破波时拒绝使用该参考。

还用湿底 dam-break 的左稀疏波、右激波解检查间断传播。它是求解器基准，不是本项目最终海啸初始条件。

## 10. 怎样复现与修改

解压后进入项目根目录：

```bash
python -m pip install -r requirements.txt
python scripts/run_verification.py
python scripts/run_experiments.py
```

Windows 可用 py 替代 python。第二条运行 11 项测试，第三条重新生成 figures/ 和 results/。完整指标见 docs/VERIFICATION_RESULTS.md 与 results/verification_metrics.json。

可从 examples/minimal_demo.py 修改 amplitude、depth 或网格数，但应一次只改一项，并检查质量残差、最低水深和 CFL。改变边界或大幅波高可能需要新的验证。

## 11. 本阶段能支持哪些结论

当前测试支持：给定基准中，离散静水平衡保持到舍入误差；质量闭合正确计入边界；小振幅传播速度接近理论；平滑非线性行波误差约按一阶下降；湿底 dam-break 误差在两级网格中下降。

它们不能支持：真实海啸预测准确、真实近岸增幅可靠、任意输入全局稳定、任意地形非负性、干湿边界已处理，或研究论文的原创性已确认。

下一阶段：独立线性求解器、边界反射量化，再加入完整有水的陆架地形做受控模型比较。

## 参考资料

1. Audusse, E., Bouchut, F., Bristeau, M.-O., Klein, R., Perthame, B. (2004). A Fast and Stable Well-Balanced Scheme with Hydrostatic Reconstruction for Shallow Water Flows. SIAM J. Scientific Computing 25(6), 2050–2065. DOI: https://doi.org/10.1137/S1064827503431090 . 作者预印本：https://www.math.univ-paris13.fr/~audusse/articles/hydro.pdf
2. Clawpack Riemann book, The shallow water equations：https://www.clawpack.org/riemann_book/html/Shallow_water.html

方法来自上述公开数学资料；推导文字与代码为本项目重新整理和实现。
