# `schakel2011_sommerfeld.py` 公式与理论逻辑梳理

本文档围绕 `schakel2011_sommerfeld.py` 中实际进入计算的公式，按“由孔隙反应输运参数到界面电磁响应波形”的逻辑链梳理理论来源、代码实现和模型取舍。公式主要来自三篇附件论文：

- Schakel & Smeulders (2010, JASA): *Seismoelectric reflection and transmission at a fluid/porous-medium interface*。本文档简称 **Schakel & Smeulders 2010**。
- Schakel et al. (2011, JAP): *Laboratory measurements and theoretical modeling of seismoelectric interface response and coseismic wave fields*。本文档简称 **Schakel et al. 2011 JAP**。
- Schakel et al. (2011, Geophysics): *Seismoelectric interface response: Experimental results and forward model*。本文档简称 **Schakel et al. 2011 Geophysics**。

## 1. 这段代码解决的物理问题

`schakel2011_sommerfeld.py` 的中心问题不是重新推导一套新震电理论，而是把反应输运输出的孔隙结构和流体化学演化，放入 Schakel 系列论文的流体-多孔介质界面震电正演框架中，合成零偏移距附近的界面 EM 波形。

它的完整计算链可以写成：

$$
\{\phi,k_0,\alpha_\infty,c_H\}(t_d)
\rightarrow
\{k(\omega),L(\omega),\sigma(\omega),\bar{\varepsilon}(\omega)\}
\rightarrow
\{R^E,T^{TM},T^{Pf}\}
\rightarrow
\hat{\phi}^E(\omega,z)
\rightarrow
u(z,t).
$$

其中 $t_d$ 是反应输运或溶蚀时间，通常为秒到小时；$t$ 是声波激发后的波形时间，通常为微秒。代码中二者严格分开：`Time_s` 只用于逐个溶蚀状态更新材料参数，`waveform time` 只用于频域响应反变换。

从故事线看，主角是 **界面电流不平衡产生的 interface EM response**。孔隙率、渗透率、曲折度、H+ 浓度、动态渗透率、电动耦合系数、界面矩阵和 Sommerfeld 积分都是支撑这个主角的证据链。

## 2. 符号、相位和坐标约定

三篇 Schakel 论文均采用频域时间因子：

$$
\exp(i\omega t).
$$

空间传播项采用：

$$
\exp(-i\mathbf{k}\cdot\mathbf{x}).
$$

代码中的 `complex_sqrt_branch()` 和 `_frequency_response_for_receivers()` 都围绕这一相位约定写成。为了保证波沿传播方向衰减，复垂向波数的平方根分支选择为使衰减项具有正确符号。对于流体侧反射 EM 波，论文写为 $\operatorname{Im}(k_z^E)<0$，并在波场中使用 $\exp(i k_z^E z_r)$；代码相应使用

$$
k_z^E=\sqrt{\omega^2s_E^2-k_x^2},
$$

并用 `exp(1j * z_fluid * k_z^E)` 传播到流体侧接收点。

坐标上，代码把界面置于 $z=0$：

- $z<0$ 为上覆流体侧，接收的是反射 interface EM 响应；
- $z>0$ 为多孔介质侧，接收的是透射 TM 电势项；
- 声源距界面的正距离为 `z_s`，但在 Sommerfeld 公式中使用有符号坐标 $z_s^{paper}<0$，代码写为 `source_z_m = -abs(cfg.z_s)`。

## 3. 反应输运变量如何进入 Schakel 模型

### 3.1 输入变量

对每个溶蚀时刻 $t_d$，代码从反应输运表中读取：

$$
\phi(t_d),\quad k_0(t_d),\quad \alpha_\infty(t_d),\quad c_H(t_d).
$$

对应代码变量为：

- `Porosity` $\rightarrow \phi$；
- `Permeability_mD` $\rightarrow k_0$；
- `Tortuosity` $\rightarrow \alpha_\infty$；
- `OutletHConc` $\rightarrow c_H$。

渗透率单位换算为：

$$
k_0[\mathrm{m^2}]
=k_0[\mathrm{mD}]\times 9.869233\times10^{-16}.
$$

H+ 浓度若以 $\mathrm{mol/cm^3}$ 给出，换算为：

$$
c_H[\mathrm{mol/L}]=1000\,c_H[\mathrm{mol/cm^3}].
$$

随后计算：

$$
c_H^\ast=\max(c_H,10^{-7}\ \mathrm{mol/L}),
\qquad
pH=-\log_{10}c_H^\ast.
$$

这些单位换算和下限处理是代码的反应输运接口，不是 Schakel 论文中的新公式。

### 3.2 电解质浓度和 zeta 电位

若输入表提供 `ElectrolyteConcentration_molL`，代码直接使用该浓度。否则采用一阶近似：

$$
C=C_{\mathrm{bg}}+c_H^\ast.
$$

这是当前 RT-SE 映射假设，不是 Schakel 论文的实测约束。

zeta 电位使用 Schakel & Smeulders 2010 Appendix A Eq. (A5)：

$$
\zeta=
\left[0.010+0.025\log_{10}(C)\right]
\frac{pH-2}{5}.
$$

在代码中对应 `electrochemistry_from_h()`。

## 4. 动态渗透率、耦合系数和动态电导率

这一组公式来自 Schakel & Smeulders 2010 Appendix A。代码函数为 `dynamic_coefficients()`。

### 4.1 动态渗透率

Schakel & Smeulders 2010 Eq. (A1) 给出 Johnson et al. 型动态渗透率：

$$
k(\omega)
=
k_0
\left[
\sqrt{1+i\frac{\omega}{\omega_t}\frac{M}{2}}
+i\frac{\omega}{\omega_t}
\right]^{-1}.
$$

过渡角频率来自 Schakel & Smeulders 2010 Eq. (A2)：

$$
\omega_t=
\frac{\phi\eta}{\alpha_\infty k_0\rho_f}.
$$

相似参数和孔隙几何长度满足 Schakel & Smeulders 2010 Eq. (A3)：

$$
M=\frac{8\alpha_\infty k_0}{\phi\Lambda^2},
\qquad
\Lambda=\sqrt{\frac{8\alpha_\infty k_0}{\phi M}}.
$$

代码默认 `M_similarity=1.0`，因此用反应输运给定的 $\phi,k_0,\alpha_\infty$ 反推 $\Lambda$。

### 4.2 Debye 长度

Schakel & Smeulders 2010 Eq. (A11) 给出：

$$
d=
\left[
\sum_l
\frac{(e z_l)^2N_l}
{\varepsilon_0\varepsilon_f k_BT}
\right]^{-1/2}.
$$

代码采用二元对称电解质：

$$
z_1=1,\quad z_2=-1,\quad N_1=N_2=1000\,C\,N_A.
$$

其中 $C$ 单位为 $\mathrm{mol/L}$，乘以 $1000N_A$ 后得到 $\mathrm{m^{-3}}$。

### 4.3 电动耦合系数

Schakel & Smeulders 2010 Eq. (A4) 给出：

$$
L(\omega)=
-
\frac{\phi}{\alpha_\infty}
\frac{\varepsilon_0\varepsilon_f\zeta}{\eta}
\left(1-\frac{2d}{\Lambda}\right)
\left[
1+
2i\frac{\omega}{M\omega_t}
\left(1-\frac{2d}{\Lambda}\right)^2
\left(1+d\sqrt{\frac{i\omega\rho_f}{\eta}}\right)^2
\right]^{-1/2}.
$$

这一步是本文模型的机制核心之一：溶蚀改变 $k_0$、$\phi$、$\alpha_\infty$ 和 $C/pH$，这些量共同改变 $L(\omega)$，进而改变界面处机械波到电磁波的转换强度。

### 4.4 流体电导率、过量电导和动态体电导率

Schakel & Smeulders 2010 Eq. (A7)：

$$
\sigma_f=\sum_l(ez_l)^2b_lN_l.
$$

若输入表提供 `FluidConductivity_S_m`，代码以该实测/校准值覆盖上述估计值。

Schakel & Smeulders 2010 Eq. (A8)：

$$
C_{\mathrm{em}}
=
2d\sum_l(ez_l)^2b_lN_l
\left[
\exp\left(-\frac{ez_l\zeta}{2k_BT}\right)-1
\right].
$$

Schakel & Smeulders 2010 Eq. (A10)：

$$
P_{\mathrm{os}}
=
\frac{8k_BTd^2}{\varepsilon_0\varepsilon_f\zeta^2}
\sum_lN_l
\left[
\exp\left(-\frac{ez_l\zeta}{2k_BT}\right)-1
\right].
$$

Schakel & Smeulders 2010 Eq. (A9)：

$$
C_{\mathrm{os}}(\omega)
=
\frac{(\varepsilon_0\varepsilon_f)^2\zeta^2}{2d\eta}
P_{\mathrm{os}}
\left[
1+\frac{2}{P_{\mathrm{os}}}d\sqrt{\frac{i\omega\rho_f}{\eta}}
\right]^{-1}.
$$

Schakel & Smeulders 2010 Eq. (A6)：

$$
\sigma(\omega)=
\frac{\phi\sigma_f}{\alpha_\infty}
\left[
1+\frac{2(C_{\mathrm{em}}+C_{\mathrm{os}}(\omega))}
{\sigma_f\Lambda}
\right].
$$

### 4.5 有效介电常数

Schakel & Smeulders 2010 正文 Eq. (20) 定义：

$$
\bar{\varepsilon}(\omega)
=
\varepsilon
-i\frac{\sigma(\omega)}{\omega}
+i\frac{\eta L^2(\omega)}
{\omega k(\omega)}.
$$

代码中体介电参数采用 Schakel & Smeulders 2010 Eq. (7) 后的混合表达：

$$
\varepsilon
=
\varepsilon_0
\left[
\frac{\phi(\varepsilon_f-\varepsilon_s)}{\alpha_\infty}
+\varepsilon_s
\right].
$$

### 4.6 上覆流体电导率和孔隙流体电导率的区别

代码中有两个容易混淆的电导率：

- $\sigma(\omega)$：多孔介质内部的动态体电导率，由 Schakel & Smeulders 2010 Appendix A 计算，进入 $\bar{\varepsilon}(\omega)$ 和多孔介质波模；
- $\sigma_{fl}$：上覆自由流体的电导率，进入流体侧 EM 慢度 $s_E^2$、流体侧 $k_z^E$ 和 Sommerfeld 传播。

在 `se_coefficients()` 中，上覆流体 EM 慢度写为：

$$
s_E^2(\omega)
=
\mu_0\varepsilon_0\varepsilon_{fl}
-i\frac{\mu_0\sigma_{fl}}{\omega}.
$$

这与 Schakel & Smeulders 2010 Appendix B 后对 $s_E^2$ 的定义、以及 Schakel et al. 2011 Geophysics Eq. (5) 上方对 $c_E(\omega)$ 的定义一致。代码的模型开关为：

$$
\sigma_{fl}=
\begin{cases}
\sigma_{fl}^{default}, & \texttt{upper\_fluid\_conductivity\_mode=constant},\\
\max\{\operatorname{Re}(\sigma_f),\sigma_{fl}^{default}\}, & \texttt{upper\_fluid\_conductivity\_mode=dynamic\_pore\_fluid}.
\end{cases}
$$

默认模式 `constant` 更接近 Schakel 论文中表格给定上覆流体参数的做法；`dynamic_pore_fluid` 是项目级敏感性假设，表示上覆流体电导率也随反应输运估计的孔隙流体电导率变化。

## 5. 孔弹性系数和波模慢度

这一组公式来自 Schakel & Smeulders 2010 正文。代码函数为 `biot_elastic_coefficients()` 和 `wave_slownesses()`。

### 5.1 Biot 弹性系数

Schakel & Smeulders 2010 Eq. (8)-(10) 给出 $A,Q,R$，代码进一步定义 $P=A+2G$：

$$
A=
\frac{
(1-\phi)^2K_sK_f-(1-\phi)K_bK_f+\phi K_sK_b
}
{K_f(1-\phi-K_b/K_s)+\phi K_s}
-\frac{2G}{3},
$$

$$
Q=
\frac{\phi\left[K_s(1-\phi)-K_b\right]K_f}
{K_f(1-\phi-K_b/K_s)+\phi K_s},
$$

$$
R=
\frac{\phi^2K_sK_f}
{K_f(1-\phi-K_b/K_s)+\phi K_s},
\qquad
P=A+2G.
$$

这里 $K_s$ 是矿物体积模量，$K_b$ 是排水骨架体积模量，$K_f$ 是孔隙流体体积模量，$G$ 是骨架剪切模量。

### 5.2 频率相关密度

Schakel & Smeulders 2010 Eq. (11)-(13)：

$$
\rho_{12}(\omega)
=
\phi\rho_f
\left[
1+i\frac{\phi\eta}{\omega\rho_fk(\omega)}
\right],
$$

$$
\rho_{11}(\omega)=(1-\phi)\rho_s-\rho_{12}(\omega),
\qquad
\rho_{22}(\omega)=\phi\rho_f-\rho_{12}(\omega).
$$

Schakel & Smeulders 2010 Eq. (26)-(29) 把电动耦合写入有效密度：

$$
\bar{\rho}_{11}=\rho_{11}-E_K,\qquad
\bar{\rho}_{12}=\rho_{12}+E_K,\qquad
\bar{\rho}_{22}=\rho_{22}-E_K,
$$

$$
E_K(\omega)
=
\frac{\eta^2\phi^2L^2(\omega)}
{k^2(\omega)\bar{\varepsilon}(\omega)\omega^2}.
$$

代码中特别保留了 $\phi^2$，与论文 Eq. (29) 一致。

### 5.3 纵波慢度

Schakel & Smeulders 2010 Eq. (24)-(25) 给出快 P 波和慢 P 波的平方慢度 $s^2_{Pf},s^2_{Ps}$：

$$
s_l^2(\omega)
=
\frac{-d_1(\omega)}{2d_2}
\mp
\frac{1}{2}
\sqrt{
\left[\frac{d_1(\omega)}{d_2}\right]^2
-4\frac{d_0(\omega)}{d_2}
},
$$

其中 $l=Pf,Ps$，

$$
d_0=\bar{\rho}_{11}\bar{\rho}_{22}-\bar{\rho}_{12}^2,
$$

$$
d_1=-
\left[
P\bar{\rho}_{22}
+R\bar{\rho}_{11}
-2Q\bar{\rho}_{12}
\right],
$$

$$
d_2=PR-Q^2.
$$

代码按绝对值大小排序两个根，较快的作为 `s2_Pf`，较慢的作为 `s2_Ps`。

### 5.4 横向波慢度

Schakel & Smeulders 2010 Eq. (30) 给出横向 EM/TM 波和 SV 波的平方慢度，仍用 Eq. (24) 的二次根结构，但系数变为：

$$
d_0=
\mu\bar{\varepsilon}
\frac{\bar{\rho}_{11}\bar{\rho}_{22}-\bar{\rho}_{12}^2}{G},
$$

$$
d_1=
-\mu\bar{\varepsilon}\bar{\rho}_{22}
-\frac{\rho_{11}\rho_{22}-\rho_{12}^2}{G},
$$

$$
d_2=\rho_{22}.
$$

代码按绝对值大小排序，较小根作为 `s2_TM`，另一根作为 `s2_SV`。

### 5.5 流固幅值比和电场幅值比

Schakel & Smeulders 2010 Eq. (36)-(37)：

$$
\beta^m(\omega)
=
\frac{\hat{\phi}_f^m}{\hat{\phi}_s^m}
=
\frac{\bar{\rho}_{11}-Ps_m^2}
{Qs_m^2-\bar{\rho}_{12}},
\qquad m=Pf,Ps,
$$

$$
\beta^n(\omega)
=
\frac{\hat{\psi}_f^n}{\hat{\psi}_s^n}
=
\frac{Gs_n^2-(1-\phi)\rho_s}{\phi\rho_f},
\qquad n=TM,SV.
$$

Schakel & Smeulders 2010 Eq. (38)-(39)：

$$
\alpha^m(\omega)
=
\frac{\hat{\phi}_E^m}{\hat{\phi}_s^m}
=
\frac{\eta\phi L(\omega)}
{k(\omega)\bar{\varepsilon}(\omega)}
\left[1-\beta^m(\omega)\right],
$$

$$
\alpha^n(\omega)
=
\frac{\hat{\psi}_E^n}{\hat{\psi}_s^n}
=
\frac{\mu\eta\phi L(\omega)}
{k(\omega)\left[\mu\bar{\varepsilon}(\omega)-s_n^2(\omega)\right]}
\left[1-\beta^n(\omega)\right].
$$

这四组幅值比将机械位移势和电磁势联系起来，是后面把 $T^{TM}$ 转为电势项的关键。

## 6. 界面边界条件和六阶线性系统

这一部分来自 Schakel & Smeulders 2010 Section III 和 Appendix B。代码函数为 `se_coefficients()`。

### 6.1 波场展开

流体入射 P 波、流体反射 P 波、流体反射 EM 波、多孔介质中的 $Pf,Ps,TM,SV$ 波分别用 Schakel & Smeulders 2010 Eq. (40)-(44) 的平面波形式表示。核心传播相位为：

$$
\hat{\phi}_I^{fl}
=
\tilde{\phi}_I^{fl}
\exp[-i(k_1x_1+k_3^{fl}x_3)],
$$

$$
\hat{\phi}_R^{fl}
=
\tilde{\phi}_R^{fl}
\exp[-i(k_1x_1-k_3^{fl}x_3)],
$$

$$
\hat{\Psi}_R^{fl}
=
(0,\tilde{\psi}_R^{fl}\exp[-i(k_1x_1-k_3^Ex_3)],0)^T,
$$

$$
\hat{\phi}_s^m
=
\tilde{\phi}_s^m
\exp[-i(k_1x_1+k_3^mx_3)],
$$

$$
\hat{\Psi}_s^n
=
(0,\tilde{\psi}_s^n\exp[-i(k_1x_1+k_3^nx_3)],0)^T.
$$

### 6.2 反射和透射系数定义

Schakel & Smeulders 2010 Eq. (45) 定义六个未知系数：

$$
R^E=\frac{\tilde{\psi}_R^{fl}}{\tilde{\phi}_I^{fl}},
\quad
R^M=\frac{\tilde{\phi}_R^{fl}}{\tilde{\phi}_I^{fl}},
\quad
T^{Pf}=\frac{\tilde{\phi}_s^{Pf}}{\tilde{\phi}_I^{fl}},
$$

$$
T^{Ps}=\frac{\tilde{\phi}_s^{Ps}}{\tilde{\phi}_I^{fl}},
\quad
T^{TM}=\frac{\tilde{\psi}_s^{TM}}{\tilde{\phi}_I^{fl}},
\quad
T^{SV}=\frac{\tilde{\psi}_s^{SV}}{\tilde{\phi}_I^{fl}}.
$$

其中代码最关心的是：

- `R_E`: 流体侧反射 EM 矢量势系数；
- `T_TM`: 多孔介质侧透射 TM 矢量势系数；
- `T_Pf`: 多孔介质侧快 P 波透射系数，用于可选 Pf 共震电势项。

### 6.3 边界条件

Schakel & Smeulders 2010 Eq. (31)-(35) 给出 open-pore 边界条件：

$$
(1-\phi)\hat{u}_3+\phi\hat{U}_3=\hat{U}_3^{fl},
$$

$$
\hat{p}=\hat{p}^{fl},
$$

$$
\hat{\sigma}_{13}=\hat{\sigma}_{33}=0,
$$

$$
\hat{H}_2=\hat{H}_2^{fl},
$$

$$
\hat{E}_1=\hat{E}_1^{fl}.
$$

它们分别表示：法向体积位移连续、压力连续、骨架界面剪应力和法向应力消失、切向磁场连续、切向电场连续。

### 6.4 六阶矩阵

Schakel & Smeulders 2010 Eq. (46) 和 Appendix B Eq. (B1) 给出：

$$
\mathbf{A}
\begin{bmatrix}
R^E & R^M & T^{Pf} & T^{Ps} & T^{TM} & T^{SV}
\end{bmatrix}^T
=
\begin{bmatrix}
k_3^{fl} & \phi\rho_{fl} & 0 & 0 & 0 & 0
\end{bmatrix}^T.
$$

代码逐项填入 Appendix B Eq. (B2)-(B7)。为便于核对，矩阵的物理分组如下：

- 第 1 行，体积位移连续，对应 Appendix B Eq. (B2)；
- 第 2 行，压力连续，对应 Appendix B Eq. (B3)；
- 第 3-4 行，界面应力条件，对应 Appendix B Eq. (B4)-(B5)；
- 第 5 行，切向磁场连续，对应 Appendix B Eq. (B6)；
- 第 6 行，切向电场连续，对应 Appendix B Eq. (B7)。

其中 Appendix B 中的两个组合项为：

$$
N_1(\omega)=
P-\frac{Q(1-\phi)}{\phi}
+
\left[
Q-\frac{R(1-\phi)}{\phi}
\right]\beta^{Pf}(\omega),
$$

$$
N_2(\omega)=
P-\frac{Q(1-\phi)}{\phi}
+
\left[
Q-\frac{R(1-\phi)}{\phi}
\right]\beta^{Ps}(\omega).
$$

代码解出的向量 `x` 依次为：

$$
x=[R^E,R^M,T^{Pf},T^{Ps},T^{TM},T^{SV}]^T.
$$

## 7. 从 2010 位移归一化系数到 2011 压力归一化系数

`schakel2011_sommerfeld.py` 不是直接把 2010 的 $R^E$ 放进 2011 Sommerfeld 积分，而是进行压力归一化。这个转换来自 Schakel et al. 2011 JAP Table I，也在 Schakel et al. 2011 Geophysics Eq. (3) 下方文字中说明：Schakel & Smeulders 2010 的系数是 displacement normalized，而 2011 forward model 中的 $R^E$ 是 pressure normalized，单位为 V/Pa。

代码中：

$$
R^E_{\mathrm{press}}
=
\frac{R^E_{2010}}{\rho_{fl}\omega^2}.
$$

对应 `pressure_normalized_re_from_coeff()`。

Schakel et al. 2011 JAP Table I 还给出 TM 透射电势项：

$$
T^{TM}_{f,\mathrm{press}}
=
\frac{\alpha^{TM}T^{TM}_{2010}}
{\rho_{fl}\omega^2}.
$$

JAP Table I 给出的是压力归一化幅值关系，不包含当前代码中的前导负号。代码采用的标量电势方向转换写为：

$$
\Phi^{TM}_{\mathrm{term}}
=
-
\frac{\alpha^{TM}T^{TM}}
{\rho_{fl}\omega^2}.
$$

Pf 共震项为：

$$
\Phi^{Pf}_{\mathrm{term}}
=
-
\frac{\alpha^{Pf}T^{Pf}}
{\rho_{fl}\omega^2}.
$$

对应 `pressure_normalized_porous_terms_from_coeff()`。因此，论文来源应写为“$\alpha T/(\rho_f\omega^2)$ 的压力归一化来自 Schakel et al. 2011 JAP Table I”，而前导负号应写为“当前代码/复现实用的标量电势符号约定”。

## 8. 2011 Sommerfeld 声源和流体侧反射电势

这一部分是 `schakel2011_sommerfeld.py` 相比 Liu 2018 风格脚本最主要的变化。公式来源为 Schakel et al. 2011 Geophysics Eq. (1)-(5)，同一结构也出现在 Schakel et al. 2011 JAP Eq. (1)-(5)。

### 8.1 圆形活塞声源

Schakel et al. 2011 Geophysics Eq. (1)：

$$
\hat{p}(\omega,R,\theta)
=
\frac{A(\omega)}{R}
\exp(-ikR)D(\theta).
$$

Schakel et al. 2011 Geophysics Eq. (2)：

$$
D(\theta)
=
\frac{J_1(ka\sin\theta)}
{ka\sin\theta}.
$$

其中 $a$ 是圆形换能器半径，$J_1$ 是一阶第一类 Bessel 函数。代码把这个方向性函数嵌入 Sommerfeld 积分后，实角分支的核中出现：

$$
J_0(kr_r\sin\theta)J_1(ka\sin\theta),
$$

并在积分前系数中保留 $1/a$。这等价于把 Eq. (2) 与 Sommerfeld 变量替换后的 $k\sin\theta\,d\theta$ 合并。

### 8.2 代码中的声源谱

Schakel et al. 2011 Geophysics/JAP 使用实验测得压力波形作为 $A(\omega)$。当前代码默认没有使用实验压力记录，而是把 Ricker 子波作为物理压力历史：

$$
p(t)=p_0\,\mathrm{Ricker}(t-t_p;f_0),
$$

其中代码的 Ricker 为：

$$
\mathrm{Ricker}(\tau;f_0)
=
\left[1-2(\pi f_0\tau)^2\right]
\exp[-(\pi f_0\tau)^2].
$$

然后用 Schakel 相位约定下的因果傅里叶核计算：

$$
P(\omega)
=
\int_0^T p(\tau)\exp(-i\omega\tau)\,d\tau.
$$

最后放入 2011 压力源形式：

$$
A(\omega)=P(\omega)R_{\mathrm{ref}}W(f).
$$

其中 $R_{\mathrm{ref}}$ 默认为 `abs(z_s)`，$W(f)$ 是代码加入的平滑带通窗。这个 Ricker 压力源是代码的 RT-SE 建模选择，不是三篇论文中的新公式；论文实验模式 `fig4_digitized` 才对应 Schakel et al. 2011 Geophysics Fig. 4 的压力波形近似。

代码中的平滑带通窗为：

$$
W(f)=0,\qquad f\le f_{\min}\ \mathrm{or}\ f\ge f_{\max},
$$

低频端：

$$
W(f)=
\frac{1}{2}
-\frac{1}{2}
\cos\left[
\pi\frac{f-f_{\min}}{\Delta f}
\right],
\qquad f_{\min}<f<f_{\min}+\Delta f,
$$

高频端：

$$
W(f)=
\frac{1}{2}
-\frac{1}{2}
\cos\left[
\pi\frac{f_{\max}-f}{\Delta f}
\right],
\qquad f_{\max}-\Delta f<f<f_{\max},
$$

中间通带 $W(f)=1$。这里 $\Delta f$ 对应 `schakel_bandpass_taper_hz`。

Ricker 源起始端还乘以一个半余弦 ramp：

$$
R_{\mathrm{ramp}}(n)
=
\frac{1}{2}
\left[
1-\cos\left(\pi\frac{n}{N_{\mathrm{ramp}}-1}\right)
\right],
\qquad n=0,\ldots,N_{\mathrm{ramp}}-1.
$$

可选的 `fig4_digitized` 模式使用代码中对 Schakel et al. 2011 Geophysics Fig. 4 压力记录的视觉近似：

$$
p_{\mathrm{kPa}}(\tau)
=
58
\exp\left[
-\frac{1}{2}
\left(\frac{\tau_{\mu s}-1.78}{0.78}\right)^2
\right]
\sin\left[2\pi\,0.56(\tau_{\mu s}-0.75)\right]
+
10
\exp\left[
-\frac{1}{2}
\left(\frac{\tau_{\mu s}-3.05}{0.28}\right)^2
\right].
$$

这只是为了在没有原始实验压力记录时复现 Fig. 4 的视觉波形，论文写作中应标为近似或诊断模式。

此外，`fig4_digitized` 分支在构造 $A(\omega)$ 时使用：

$$
A_{\mathrm{fig4}}(\omega)
=
P_{\mathrm{fig4}}(\omega)\,W(f)\,
\frac{|z_s|}{0.5}.
$$

其中 `/0.5` 是旧复现脚本保留下来的幅值缩放/诊断因子，不是 Schakel et al. 2011 Geophysics 或 JAP 的理论公式。若使用该模式，应在图注或方法说明中标为 legacy calibration choice。

若 SciPy 的 Bessel 函数不可用，代码还提供数值备用定义：

$$
J_0(x)\approx
\frac{1}{N}
\sum_{j=1}^N
\exp\left(ix\cos\theta_j\right),
$$

$$
J_1(x)\approx
\frac{1}{\pi}
\int_0^\pi
\cos(\theta-x\sin\theta)\,d\theta.
$$

正常运行时优先使用 `scipy.special.j0` 和 `scipy.special.j1`。

### 8.3 流体侧 Sommerfeld 积分

Schakel et al. 2011 Geophysics Eq. (3) 先以径向波数 $k_r$ 写为：

$$
\hat{\phi}^E(\omega,r_r,z_r)
=
-iA(\omega)
\int_0^\infty
\frac{k_r}{k_z}
D(k_r)J_0(k_rr_r)
\exp(ik_zz_s)
R^E(k_r)
\exp(ik_z^Ez_r)\,dk_r.
$$

令

$$
k_r=k\sin\theta,\qquad k_z=k\cos\theta,
$$

得到 Schakel et al. 2011 Geophysics Eq. (4) 的复角积分：

$$
\hat{\phi}^E
=
-iA(\omega)
\int_0^{\pi/2+i\infty}
D(\theta)k\sin\theta J_0(kr_r\sin\theta)
\exp(ikz_s\cos\theta)
R^E(\theta)
\exp(ik_z^E(\theta)z_r)\,d\theta.
$$

其中

$$
k_z^E(\theta)
=
\omega
\sqrt{
\frac{1}{c_E^2(\omega)}
-\frac{\sin^2\theta}{c_P^2}
},
\qquad
\operatorname{Im}(k_z^E)<0.
$$

流体 EM 和声波速度为 Schakel et al. 2011 Geophysics Eq. (5) 上方定义：

$$
c_E(\omega)
=
\left[
\mu\varepsilon_0\varepsilon_{fl}
-i\frac{\mu\sigma_{fl}}{\omega}
\right]^{-1/2},
\qquad
c_P=\sqrt{\frac{K_f}{\rho_f}}.
$$

Schakel et al. 2011 Geophysics Eq. (5) 将积分路径拆成实角分支和倏逝分支：

$$
\hat{\phi}^E
=
-\frac{iA(\omega)}{a}
\int_0^{\pi/2}
J_0(kr_r\sin\theta)J_1(ka\sin\theta)
\exp(ikz_s\cos\theta)
R^E(\theta)
\exp(ik_z^E(\theta)z_r)
d\theta
$$

$$
\quad
+
\frac{A(\omega)}{a}
\int_0^\infty
J_0\!\left(kr_r\sqrt{\gamma^2+1}\right)
\frac{
J_1\!\left(ka\sqrt{\gamma^2+1}\right)
}
{\sqrt{\gamma^2+1}}
\exp(kz_s\gamma)
R^E(\gamma)
\exp(ik_z^E(\gamma)z_r)
d\gamma.
$$

代码对应关系：

- `theta = linspace(0, pi/2, n_theta)` 对应实角分支；
- `gamma = linspace(0, schakel_gamma_max, n_theta/2)` 对应倏逝分支；
- `J0(k * radial_m * ...)` 是 $J_0$；
- `J1(k * a * ...)` 是 $J_1$；
- `exp(1j * k * source_z_m * cos(theta))` 是 $\exp(ikz_s\cos\theta)$；
- `exp(k * source_z_m * gamma)` 是 $\exp(kz_s\gamma)$；
- `pressure_normalized_re_from_coeff()` 提供 $R^E_{\mathrm{press}}$；
- 最终 `response = -(1j*A/a)*first + (A/a)*second` 对应 Eq. (5) 的两个分支。

因为 `source_z_m<0`，倏逝分支的 $\exp(kz_s\gamma)$ 随 $\gamma$ 衰减。

## 9. 多孔介质侧透射电势项

多孔介质侧公式来自 Schakel et al. 2011 JAP Eq. (6)-(8)。JAP 论文的实验几何是有限厚度多孔样品，因此 Eq. (8) 包含：

1. 前界面生成并向下传播的 TM 电势项；
2. 与透射 Pf 波同行的共震电势项；
3. Pf 波在样品背面反射后产生的后界面多次项。

当前 `schakel2011_sommerfeld.py` 是单界面 RT-SE 正演，不是实验 slab 复现。因此代码只默认保留前界面 TM 项，并可选择性加入 Pf 共震项。对实角分支，代码实际计算的多孔侧频域响应为：

$$
\hat{\phi}_{por}^{(real)}(\omega,z_r)
=
\int_0^{\pi/2}
J_0(kr_r\sin\theta)J_1(ka\sin\theta)
\exp(ikz_s\cos\theta)
\left[
\Phi^{TM}_{\mathrm{term}}(\theta)
\exp(-ik_z^{TM}(\theta)z_r)
+
\Phi^{Pf}_{\mathrm{term}}(\theta)
\exp(-ik_z^{Pf}(\theta)z_r)
\right]d\theta,
\qquad z_r>0.
$$

对倏逝分支，代码使用与流体侧相同的 Sommerfeld 路径替换：

$$
\hat{\phi}_{por}^{(ev)}(\omega,z_r)
=
\int_0^{\gamma_{\max}}
J_0\left(kr_r\sqrt{\gamma^2+1}\right)
\frac{
J_1\left(ka\sqrt{\gamma^2+1}\right)
}
{\sqrt{\gamma^2+1}}
\exp(kz_s\gamma)
\left[
\Phi^{TM}_{\mathrm{term}}(\gamma)
\exp(-ik_z^{TM}(\gamma)z_r)
+
\Phi^{Pf}_{\mathrm{term}}(\gamma)
\exp(-ik_z^{Pf}(\gamma)z_r)
\right]d\gamma.
$$

最终组合权重与 Schakel et al. 2011 Eq. (5)/(8) 的结构一致：

$$
\hat{\phi}_{por}(\omega,z_r)
=
-\frac{iA(\omega)}{a}
\hat{\phi}_{por}^{(real)}
+
\frac{A(\omega)}{a}
\hat{\phi}_{por}^{(ev)}.
$$

其中前界面 TM 项用 2010 系数写入代码为：

$$
T_f^{TM}
=
-
\frac{\alpha^{TM}T^{TM}}
{\rho_{fl}\omega^2}.
$$

代码还保留一个默认关闭的诊断开关 `include_porous_pf_coseismic`。若打开，则额外加入：

$$
T_f^{Pf}
=
-
\frac{\alpha^{Pf}T^{Pf}}
{\rho_{fl}\omega^2},
\qquad
\hat{\phi}^{Pf}_t
\sim
T_f^{Pf}\exp(-ik_z^{Pf}z_r).
$$

若 `include_porous_pf_coseismic=False`，则 $\Phi^{Pf}_{\mathrm{term}}$ 在积分中被置为 0。

Schakel et al. 2011 JAP Eq. (8) 的完整实验 slab 形式还包含 $w_s$、$S^E$ 和背界面反射/多次项。当前代码没有这些项。二者关系可以概括为：

| 对象 | Schakel et al. 2011 JAP Eq. (8) | 当前代码 |
| --- | --- | --- |
| 前界面 TM interface EM | 保留 | 默认保留 |
| Pf 共震电势 | 保留 | 可选，默认关闭 |
| 背界面反射和样品宽度 $w_s$ | 保留 | 不实现 |
| 多次界面响应 $S^E$ | 保留 | 不实现 |
| 模型解释 | 有限厚度实验样品 | 单界面 RT-SE 正演 |

因此文档和论文使用时应明确写作：**本代码是单界面模型，不是 Schakel et al. 2011 JAP 有限厚度实验样品的完整多次反射模型**。

## 10. 频率响应到时间域波形

代码对每个正频率计算接收阵列频域响应 $\hat{u}(z,\omega)$，然后用 Schakel 相位约定反变换：

频率网格为线性采样：

$$
f_j=f_{\min}
+j\frac{f_{\max}-f_{\min}}{N_f-1},
\qquad j=0,\ldots,N_f-1.
$$

$$
u(z,t)
=
2\,\operatorname{Re}
\int_{f_{\min}}^{f_{\max}}
\hat{u}(z,f)
\exp(i2\pi ft)
\,df.
$$

这里的因子 2 表示只显式积分正频率，并用实部恢复实值时间波形。频率网格来自：

$$
f\in[f_{\min},f_{\max}],
$$

代码默认由 `schakel_bandpass_low_hz` 和 `schakel_bandpass_high_hz` 控制，并在源谱上施加平滑带通窗。

界面声波到时为：

$$
T_0=\frac{z_s}{v_f},
\qquad
v_f=\sqrt{\frac{K_{fl}}{\rho_{fl}}}.
$$

输出波形从 $t=0$ 开始，并在图上标记 $T_0$。这不是接收端门控；如果 $T_0$ 前出现有限带宽旁瓣，代码会保留并在 `waveform_t0_causality_diagnostics.csv` 中报告。

## 11. 接收点、界面行和峰值指标

接收点由：

$$
z_j=z_{\min}+j\Delta z,\qquad z_{\min}\le z_j\le z_{\max}.
$$

代码将它们分成：

- $z_j<0$: 流体侧 `R_E` side；
- $z_j>0$: 多孔介质侧 `T_E` side；
- $z_j=0$: 仅用于绘图的界面插值行。

界面行不是 Schakel 边界条件直接求得的接收点。代码用最近流体侧和多孔侧接收行的平均值给出：

$$
u(z=0,t)
=
\frac{1}{2}
\left[
u(z_{-},t)+u(z_{+},t)
\right].
$$

所有峰值、极性、收敛和 T0 因果性定量统计都排除 $z=0$ 行。

时间序列表中，孔弹模型有效性由代码标记为：

$$
\mathrm{valid}
=
(\phi_{\min}<\phi_{\mathrm{raw}}<\phi_{\max}).
$$

对于 $|R^E|$、$|T^{TM}|$、$|L|$、$|\sigma|$ 和波形峰值等趋势图，代码使用第一个有效状态作为归一化基准：

$$
X_{\mathrm{norm}}(t_d)
=
\frac{X(t_d)}{X(t_{d,0})},
$$

其中 $t_{d,0}$ 是第一个 `valid_poroelastic=True` 的溶蚀时刻。若基准值为 0 或非有限值，代码不强行归一化，而是写为 `NaN`。

峰值幅度定义为：

$$
A_{\max}
=
\max_{z\ne0,t}|u(z,t)|.
$$

流体侧和多孔侧分别为：

$$
A_{\max}^{R_E}
=
\max_{z<0,t}|u(z,t)|,
\qquad
A_{\max}^{T_E}
=
\max_{z>0,t}|u(z,t)|.
$$

绘图时，电势单位从 V 转为 mV：

$$
u_{\mathrm{mV}}=10^3u_{\mathrm{V}}.
$$

## 12. 数值收敛、极性和 T0 诊断

代码中的诊断公式不是 Schakel 论文的理论公式，而是为了检查正演结果是否可信。

### 12.1 收敛诊断

对不同积分级别 $N$，代码比较峰值相对变化：

$$
\Delta A_{\mathrm{rel}}^{(N)}
=
\frac{A_{\max}^{(N)}-A_{\max}^{(N_{\mathrm{prev}})}}
{A_{\max}^{(N_{\mathrm{prev}})}}.
$$

并比较峰值时间漂移：

$$
\Delta t_{\mathrm{peak}}^{(N)}
=
t_{\mathrm{peak}}^{(N)}
-t_{\mathrm{peak}}^{(N_{\mathrm{prev}})}.
$$

### 12.2 上下极性反转

Schakel et al. 2011 JAP 在 Eq. (8) 后的实验解释中指出，多孔介质内第一个界面脉冲相对流体侧计算具有反极性，这与垂直振荡电偶极子解释一致。代码用对称距离接收点检查：

$$
u(-d,t^\ast)\,u(+d,t^\ast)<0
$$

则标记为极性反转。

### 12.3 T0 前后能量比

代码报告：

$$
\mathrm{ratio}_{pre/post}
=
\frac{\max_{z\ne0,t<T_0}|u(z,t)|}
{\max_{z\ne0,t\ge T_0}|u(z,t)|}.
$$

该指标用于判断有限带宽旁瓣是否过强。它不是物理门控，也不改变波形。

## 13. 代码公式覆盖清单

| 代码位置 | 公式/理论对象 | 论文来源 | 本文档位置 |
| --- | --- | --- | --- |
| `h_conc_to_molL()` | H+ 单位换算、pH 下限 | 代码接口假设 | 第 3.1 节 |
| `electrochemistry_from_h()` | $\zeta(C,pH)$ | Schakel & Smeulders 2010 Eq. (A5) | 第 3.2 节 |
| `dynamic_coefficients()` | $k(\omega),\omega_t,M,\Lambda$ | Schakel & Smeulders 2010 Eq. (A1)-(A3) | 第 4.1 节 |
| `dynamic_coefficients()` | $d,L(\omega),\sigma_f,C_{\mathrm{em}},C_{\mathrm{os}},\sigma(\omega)$ | Schakel & Smeulders 2010 Eq. (A4), (A6)-(A11) | 第 4.2-4.4 节 |
| `dynamic_coefficients()` | $\bar{\varepsilon}(\omega)$ | Schakel & Smeulders 2010 Eq. (20) | 第 4.5 节 |
| `_upper_fluid_sigma()` | 上覆流体 $\sigma_{fl}$ 的 constant/dynamic 选择 | Schakel 上覆流体 EM 慢度公式 + 项目级开关 | 第 4.6 节 |
| `biot_elastic_coefficients()` | $A,Q,R,P$ | Schakel & Smeulders 2010 Eq. (8)-(10) | 第 5.1 节 |
| `wave_slownesses()` | $\rho_{ij},\bar{\rho}_{ij},E_K$ | Schakel & Smeulders 2010 Eq. (11)-(13), (26)-(29) | 第 5.2 节 |
| `wave_slownesses()` | $s^2_{Pf},s^2_{Ps},s^2_{TM},s^2_{SV}$ | Schakel & Smeulders 2010 Eq. (24)-(25), (30) | 第 5.3-5.4 节 |
| `wave_slownesses()` | $\beta,\alpha$ 幅值比 | Schakel & Smeulders 2010 Eq. (36)-(39) | 第 5.5 节 |
| `se_coefficients()` | 六未知量、边界条件、矩阵方程 | Schakel & Smeulders 2010 Eq. (31)-(46), Appendix B | 第 6 节 |
| `complex_sqrt_branch()` | 复波数衰减分支选择 | Schakel 相位约定和 $\operatorname{Im}(k_z^E)<0$ | 第 2 节 |
| `pressure_normalized_re_from_coeff()` | $R^E/(\rho_f\omega^2)$ | Schakel et al. 2011 JAP Table I; Schakel et al. 2011 Geophysics Eq. (3) 说明 | 第 7 节 |
| `pressure_normalized_porous_terms_from_coeff()` | $\alpha^{TM}T^{TM}/(\rho_f\omega^2)$ 和 $\alpha^{Pf}T^{Pf}/(\rho_f\omega^2)$ 的归一化；代码额外采用负号 | 归一化来自 Schakel et al. 2011 JAP Table I, Eq. (6)-(8)；负号来自代码符号约定 | 第 7、9 节 |
| `schakel_source_reference_distance_m()` | $R_{\mathrm{ref}}=|z_s|$ 或用户指定参考距离 | 代码压力源标定 | 第 8.2 节 |
| `schakel_source_A_spectrum()` | $A(\omega)$ 压力源谱 | Schakel et al. 2011 Geophysics/JAP Eq. (1)，代码 Ricker/Fig. 4 近似假设 | 第 8.1-8.2 节 |
| `_bandpass_taper()` / `_source_time_taper()` | 平滑带通窗、源起始半余弦 ramp | 代码数值实现 | 第 8.2 节 |
| `_source_time_table()` | Schakel Fig. 4 压力波形视觉近似 | Schakel et al. 2011 Geophysics Fig. 4；代码近似 | 第 8.2 节 |
| `_j0()` / `_j1()` | Bessel 函数及无 SciPy 时的积分近似 | Sommerfeld 积分所需数值函数 | 第 8.2 节 |
| `_frequency_response_for_receivers()` | Sommerfeld 实角和倏逝分支 | Schakel et al. 2011 Geophysics Eq. (3)-(5); Schakel et al. 2011 JAP Eq. (3)-(5) | 第 8.3 节 |
| `_frequency_response_for_receivers()` | 多孔侧前界面 TM 和可选 Pf 项 | Schakel et al. 2011 JAP Eq. (6)-(8) | 第 9 节 |
| `_frequency_grid()` | 线性频率网格 $f_j$ | 代码数值积分实现 | 第 10 节 |
| `synthesize_waveforms_schakel2011()` | 正频率反变换、$T_0$ | Schakel 相位约定；代码数值实现 | 第 10 节 |
| `compute_time_series()` | 首个有效状态归一化、孔弹有效性掩膜 | 代码结果整理 | 第 11-12 节 |
| peak/diagnostic functions | 峰值、极性、T0 前后比、收敛 | 代码诊断；JAP 极性解释 | 第 11-12 节 |

## 14. 与论文一致的部分和主动偏离的部分

### 与论文一致

1. 动态渗透率、电动耦合系数、动态电导率和有效介电常数遵循 Schakel & Smeulders 2010 Appendix A 与 Eq. (20)。
2. 孔弹性系数、有效密度、快慢纵波、TM/SV 横向波、$\alpha/\beta$ 幅值比遵循 Schakel & Smeulders 2010 正文公式。
3. 界面反射/透射系数由 Schakel & Smeulders 2010 的 open-pore 边界条件和 Appendix B 六阶矩阵求解。
4. 流体侧反射电势采用 Schakel et al. 2011 Geophysics/JAP 的 Sommerfeld Eq. (5) 结构，包括有限圆形换能器 Bessel 方向性和倏逝分支。
5. 2010 系数到 2011 压力归一化系数的转换遵循 Schakel et al. 2011 JAP Table I。

### 代码主动偏离或简化

1. 声源谱默认不是实验记录，而是物理压力 Ricker 波形的因果傅里叶积分；这是 RT-SE 正演选择，不是论文实验复现。
2. 代码强制 `offset_D=0`，关注零偏移距/同轴接收线。
3. 多孔介质侧只默认保留单界面前界面 TM 项；Schakel et al. 2011 JAP Eq. (8) 中的后界面 slab 多次项没有实现。
4. Pf 共震项在代码中作为诊断选项存在，但默认关闭，以突出 interface EM response。
5. 反应输运到电解质浓度和 zeta 电位的映射包含项目级假设，论文写作时应与 Schakel 理论公式分开表述。

## 15. 可用于论文方法部分的机制表述

可以把本模型写成以下机制链：

> Dissolution changes pore connectivity and electrolyte chemistry, represented by $\phi$, $k_0$, $\alpha_\infty$, and $c_H$. These variables modify the dynamic permeability $k(\omega)$, electrokinetic coupling coefficient $L(\omega)$, and dynamic conductivity $\sigma(\omega)$. Through the Schakel & Smeulders (2010) fluid/porous-medium boundary system, the modified electrokinetic state changes the pressure-normalized interface conversion coefficients. The Schakel et al. (2011) Sommerfeld integral then maps these angle- and frequency-dependent coefficients into a receiver-depth-dependent interface electromagnetic waveform.

用中文表达为：

> 溶蚀首先改变孔隙连通性和流体化学状态，这些变化通过 $\phi,k_0,\alpha_\infty,c_H$ 进入动态渗透率、电动耦合系数和动态电导率；随后，Schakel & Smeulders (2010) 的流体-多孔介质界面边界条件把这些材料参数变化转化为压力归一化的界面 EM 转换系数；最后，Schakel et al. (2011) 的 Sommerfeld 积分把角度-频率相关的转换系数合成为接收点随深度变化的 interface EM 波形。

这条线避免把公式写成列表，而是把每组公式放在“溶蚀如何改变界面 EM response”的机制链中。
