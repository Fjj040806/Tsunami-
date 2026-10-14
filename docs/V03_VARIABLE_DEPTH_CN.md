# V0.3：变水深方程、良平衡离散与静水保持验证

JIajun FU；版本 0.3.0；2026-10-02。

本阶段完成变水深线性及非线性模型的推导、全湿二阶变底算法实现，以及静水、质量收支和光滑精确解测试。这里的“验证”是数值算法 verification；没有进行实验数据 validation，也还没有得到大陆架浅化实验的研究结论。

## 1. 几何约定与非线性方程

取静止参考水面为 z=0，床面高程为 b(x)<0，静水深 D(x)=-b(x)>0。扰动水面为 eta(x,t)，总水深 h=eta-b=D+eta，单位宽度流量 q=hu。床面固定，故 h_t=eta_t。

质量守恒仍为

$$h_t+q_x=0.$$

深度积分的压力项与床面作用给出动量平衡律

$$q_t+\partial_x\left(\frac{q^2}{h}+\frac12gh^2\right)=-gh b_x=gh D_x.$$

在光滑全湿区域，把 q=hu 代入，并使用质量方程，可化为

$$\eta_t+\partial_x[(D+\eta)u]=0,\qquad
u_t+u u_x+g\eta_x=0.$$

床面导数没有消失：它包含在第一式中的 D(x)，以及保守动量方程的源项中。对于不连续床面，不能随意把 b_x 当普通函数；本项目的非线性离散采用静水重构给出界面床面作用的离散解释。

静水解为 u=q=0、eta=C（任意使全域保持湿润的常数）。此时 h=C-b，压力梯度满足

$$\partial_x\left(\frac12gh^2\right)=gh h_x=-gh b_x.$$

通量梯度和源项必须相抵消。仅分别对两者做“看起来二阶”的近似，不保证离散静水平衡。

## 2. 变水深线性化

在 eta=0、u=0、h=D(x) 附近，舍去 (eta u)_x 和 u u_x，得到

$$\eta_t+\partial_x(Du)=0,\qquad u_t+g\eta_x=0.$$

特别地，第一式是 eta_t+D u_x+D_x u=0，不能漏掉 D_x u。将 Q=Du 作为线性流量时，对固定 D 有 Q_t+gD eta_x=0；它一般不等于 Q_t+(gD eta)_x=0。代码直接演化单元平均的 (eta,u)，通量是 (Du,g eta)。

光滑 D 下，连续线性系统有能量恒等式

$$\partial_t\left(\frac12 g\eta^2+\frac12 D u^2\right)
+\partial_x(gD\eta u)=0.$$

这是连续方程的性质；本阶段没有把它宣称为离散算法的精确能量守恒性质。

在一个水深跳跃界面，采用公共水面 eta_star 与公共流量 Q_star，设两侧 c_l=sqrt(gD_l)、c_r=sqrt(gD_r)，则入射特征给出

$$\eta_* = \frac{c_l\eta_l+c_r\eta_r+Q_l-Q_r}{c_l+c_r},$$

$$Q_* = \frac{c_r Q_l+c_l Q_r+c_lc_r(\eta_l-\eta_r)}{c_l+c_r}.$$

通量为 (Q_star,g eta_star)。当 D_l=D_r 时，它还原为 V0.2 的常系数精确迎风通量；代码已有退化一致性测试。

若左侧入射波在界面处的振幅为 A，右侧无入射波，水面连续与流量连续给出水面振幅系数

$$R=\frac{c_l-c_r}{c_l+c_r},\qquad
T=\frac{2c_l}{c_l+c_r},\qquad
R^2+\frac{c_r}{c_l}T^2=1.$$

这些是线性台阶参考关系。当前已检查界面通量满足这一代数关系；尚未用完整波包传播实验测量反射与透射系数。

## 3. 非线性二阶静水重构

每个单元从平均 h、q、b 构造 eta=h+b 和速度 q/h。eta 与速度使用 MC 限制器重构；光滑床面使用中心斜率重构 b。得到单元内部左、右迹 (eta_i,l,u_i,l,b_i,l)、(eta_i,r,u_i,r,b_i,r)，并取 h_i,l/r=eta_i,l/r-b_i,l/r。

光滑床面采用中心斜率有助于使相邻两侧床面迹之差具有 O(dx³) 的光滑展开，这是避免静水修正降低精度的重要条件。MC 限制器仍会影响水面和速度极值附近的局部阶数。对与网格面严格对齐的台阶，选择 bed_reconstruction='constant'，不跨台阶进行中心床面重构；不对台阶附近宣称全局二阶精度。

在单元界面，令 b_star=max(b_l,b_r)，h_l_star=eta_l-b_star、h_r_star=eta_r-b_star，以 (h_star,h_star u) 调用 HLL 通量 F_star。两侧压力修正分别为

$$F_l=F_*+\begin{pmatrix}0\\\frac g2(h_l^2-h_{l,*}^2)\end{pmatrix},\qquad
F_r=F_*+\begin{pmatrix}0\\\frac g2(h_r^2-h_{r,*}^2)\end{pmatrix}.$$

单元内部还必须保留中心地形源项

$$S_i^c=\begin{pmatrix}0\\
\frac g2(h_{i,l}+h_{i,r})(b_{i,l}-b_{i,r})\end{pmatrix}.$$

半离散式为

$$\frac{dU_i}{dt}=-\frac{F_{i+1/2,l}-F_{i-1/2,r}}{\Delta x}+\frac{S_i^c}{\Delta x}.$$

该中心源项与压力修正遵循 Audusse et al. (2004) 的式 (3.8)–(3.9)。本项目重构限制和全湿检查是具体实现选择，未继承或证明所有湿干、正性及熵性质。

**静水相消证明。** eta=C、u=0 时，每个修正界面的动量通量恢复到本单元对应迹的 gh²/2。因此通量差为 g(h_i,r²-h_i,l²)/2。而 h_i,r-h_i,l=b_i,l-b_i,r，故 S_i^c 正好等于该压力通量差，半离散导数为零。SSP-RK2 的两个阶段均使用同一良平衡算子，因此保持这一平衡，实际运算只剩浮点误差。

若只加界面压力修正而漏掉 S_i^c，本次光滑静水场的最大瞬时流量导数绝对值约为 1.29259 m²/s²。这个量来自删项后的半离散残差分析，不是实际传播实验，也不是速度加速度。

## 4. 时间推进、全湿范围和质量收支

采用 SSP-RK2 与 CFL=0.4。非线性速度估计包括单元、重构物理迹和 HLL 界面；线性采用相应 sqrt(gD)。第二阶段重新检查 CFL，必要时减小试探步长。

物理单元、阶段状态和重构迹均须满足全湿条件。为保持本阶段范围清晰，静水重构后的界面也要求正水深；若近干即报错，不进行水深裁剪。中心床面重构可能在强变化处超调，故只用于光滑且充分湿润的地形；台阶使用常数床面重构。当前不存在对任意地形与近干波动的正性保证。

质量通量在相邻单元之间共享。地形源项的质量分量为零。开边界质量积分采用两 RK 阶段通量的平均；制造解的外加质量源同样以两阶段加权积分。残差定义为

$$\mathcal R(t)=M(t)-M(0)+\int_0^t(F_h^{right}-F_h^{left})\,dt-
\int_0^t\int f_h\,dx\,dt,$$

并除以初始水体积（单位宽度）。很小的收支残差只说明离散账目一致，不等同于物理模型准确或波形误差很小。零梯度 outflow 边界不是严格无反射边界。

## 5. 验证结果

### 5.1 静水保持

L=1000 m，N=400，eta=0.3 m，q=u=0，周期边界，运行 100 s。光滑水深为 D=10+2 cos(2 pi x/L)，使用单元平均；台阶为 x<500 时 D=10 m，否则 D=4 m，台阶与网格面重合。

| 地形 | 模型 | 最大水面漂移 / m | 最大速度 / (m/s) |
|---|---|---:|---:|
|光滑|线性|3.33067e-16|9.45465e-16|
|光滑|非线性|1.34892e-14|2.37825e-14|
|台阶|线性|0|4.48369e-17|
|台阶|非线性|7.21645e-16|0|

上述四个实验的质量残差在本次双精度运行中均为 0；不应解读为任意平台或计算中必定精确为 0。另有自动测试覆盖周期、开边界和反射边界下的两种静水地形。

### 5.2 光滑变底制造解

使用 method of manufactured solutions：先指定光滑解析场，再从方程残差求外加源项。它能检查运动状态中的地形、通量、时间推进实现，避免仅用静水测试遗漏错误。

取 D=10+2 cos(kx)，eta=A sin(k(x-ct))，u=(c/10)eta，A=0.2 m，k=2pi/1000，c=sqrt(9.81×10)。两模型拥有相同的解析 eta、u，但各自所需的外加源不同，因此它们的重合不是模型等价性的证据。

记 h=D+eta，d=D_x，e=eta_x，v=u_x，a=eta_t，w=u_t，则制造源为

$$f_L=\begin{pmatrix}a+du+Dv\\w+ge\end{pmatrix},$$

$$f_{NL}=\begin{pmatrix}a+(d+e)u+hv\\
au+hw+(d+e)u^2+2huv+ghe\end{pmatrix}.$$

初值、参考解、床面和制造源均以 8 点 Gauss 积分取单元平均。另用复步长微分独立检查解析源公式确实等于保守方程残差；不是重复同一公式做自我验证。参考解没有调用数值求解器。

t=20 s 时的平均 L1 水面误差为：

| N | 线性误差 / m | 非线性误差 / m |
|---:|---:|---:|
|50|5.44199e-4|5.30496e-4|
|100|1.33775e-4|1.30323e-4|
|200|3.05227e-5|3.14585e-5|
|400|6.93444e-6|7.40595e-6|
|800|1.53347e-6|1.73156e-6|

最后一对网格的观测阶为 2.177 与 2.097，支持该光滑基准平均 L1 的近二阶收敛。这不是所有范数、所有地形、冲击波或长期传播的统一精度定理。制造源测试不验证真实海啸源、地形或观测。

### 5.3 自动测试

累计 27 项测试成功（V0.2 的 19 项加 8 项变底测试）。新增覆盖：两地形静水与三边界；变底制造解收敛；解析源的独立残差检查；线性退化到平底的一致性；台阶界面传输与能量系数恒等式；运动状态质量和 CFL；非线性平底简单波精确解收敛；输入检查和数组不修改。精确解收敛测试使用观测阶 >1.7 作为软件回归门槛，不作为数学证明。

## 6. 使用与下一步

```bash
python -m pip install -r requirements.txt
python scripts/run_verification.py
python scripts/run_v03.py
```

新代码在 src/tsunami/variable.py；解析制造解在 manufactured.py；配置在 configs/variable_verification.json。结果表、场数据及运行记录在 results；两张科学图在 figures。旧版本实验保留在 v0.1_results/v0.1_figures 与 v0.2_results/v0.2_figures。English 方法草稿在 paper/methods_v0.3.tex。

本阶段已达到进入地形传播实验的基础门槛；没有完成全部 V0.3 地形研究。下一步先发送低幅波包穿过单个对齐台阶，测量反射与透射，检查与本报告 R/T 参考关系一致，并使边界回波在测量时窗之外。然后构建平滑大陆架，进行波幅、水深比、地形过渡长度与传播时间的参数研究。Green 定律只在缓变、小幅、反射弱等假设下作为近似参考。实际实验数据验证、湿干爬高和最终论文整合仍在后续阶段。

## 7. 方法引用与定位

Audusse, E., Bouchut, F., Bristeau, M.-O., Klein, R., and Perthame, B. (2004). A Fast and Stable Well-Balanced Scheme with Hydrostatic Reconstruction for Shallow Water Flows. SIAM J. Sci. Comput. 25(6), 2050–2065. DOI: https://doi.org/10.1137/S1064827503431090。作者机构公开版本：https://publications.imp.fu-berlin.de/478/1/file_2004_siam.pdf。

Delestre et al. (2012), A limitation of the hydrostatic reconstruction technique for Shallow Water equations, https://arxiv.org/abs/1206.4986，指出部分坡度、网格、水深组合下的局限；这也是我们保留全湿范围与独立网格检查的理由。离散公式不是新方法，项目贡献仍是自主实现、验证与后续受控研究。本文件是阶段报告，不是可投稿原创论文。
