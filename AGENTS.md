# 震电效应模拟 (Seismoelectric Simulation) - 项目说明

本项目用于把孔隙尺度反应输运 (Reactive Transport) 输出映射到震电效应 (Seismoelectric Effect, SE) 正演模型中，重点考察溶蚀导致的孔隙率、渗透率、流体电导率、曲折度和界面/通道几何变化，是否会形成可解释的 interface EM response。

理论框架主要基于 Schakel & Smeulders (2010) 的流体-多孔介质界面边界条件，以及 Liu et al. (2018) 风格的频率-波数谱积分机制。当前项目的公式代码已基本确定，后续研究目标是发展“反应输运-水文/电导参数-震电响应”的多物理耦合论文框架，目标期刊为 JGR: Solid Earth。

## 研究主线

核心问题：

1. 将反应输运输出的孔隙率、渗透率、曲折度、H+ 浓度/流体电导率等变量映射为动态渗透率、电动耦合系数和动态电导率。
2. 在 Schakel & Smeulders (2010) 的界面边界条件下计算反射/透射震电转换系数。
3. 采用 Liu et al. (2018) 类型的频率-波数积分合成有限偏移距波形。
4. 分析溶蚀过程中水力连通性增强、界面电流不平衡变化与 interface EM response 幅值、极性、空间分布之间的关系。

需要始终区分两个时间尺度：

- `dissolution time` / `Time_s`：反应输运或溶蚀演化时间，通常为秒到小时量级。
- `waveform time`：声波激发后的传播和接收时间，通常为微秒量级。

## 主要代码

- `seismoelectric_offset_liu2018_spectral.py`：当前主模型脚本。读取反应输运表格，计算 Schakel/Pride 动态系数与界面转换系数，并用 Liu 2018 风格的频率-波数谱积分合成有限偏移距 VSEP/interface EM 波形。后续论文结果应优先以该脚本为核心。
- `liu2018_fig2b_dipole_comparison.py`：用于复现/解释 Liu et al. (2018) Fig. 2(b) 类型的空间幅度分布。比较完整 Pride 理论谱积分波形与 Liu Eq. (4) 的垂直电偶极子解释模型。注意：偶极子模型只用于解释空间指向性，不应乘回主正演波形。
- `seismoelectric_offset.py`：早期有偏移距版本。包含快速可视化和谱积分相关逻辑，可作为历史对照，但当前有限偏移距结果以 `seismoelectric_offset_liu2018_spectral.py` 为准。
- `seismoelectric_nooffset.py`：零偏移距/简化 VSEP 原型。适合快速检查反应输运参数到转换系数的趋势，不应替代完整频率-波数正演。

## 输入数据

- `global_evolution.xlsx`：反应输运导出的全局演化数据。主要输入列包括孔隙率、渗透率、曲折度、H+ 浓度等，脚本会据此构造震电模型参数。
- `uploaded_files/`：论文和附件材料目录。
  - `uploaded_files/manuscript.pdf`：此前反应输运-NMR 研究的主文稿。
  - `uploaded_files/supporting information.pdf`：此前反应输运-NMR 研究的补充材料。
  - `uploaded_files/震电理论公式/Seismoelectric reflection and transmission at a fluid_porous-medium interface.pdf`：Schakel & Smeulders (2010)，界面边界条件和转换系数的主要公式来源。
  - `uploaded_files/震电理论公式/Seismoelectric interface electromagnetic wave characteristics for the finite offset Vertical Seismoelectric Profiling configuration_ Theoretical modeling and experiment verification.pdf`：Liu et al. (2018)，有限偏移距 VSEP 频率-波数积分和电偶极子解释模型来源。
  - `uploaded_files/JGR Solid Earth - 2023 - Hu - Water Table and Permeability Estimation From Multi‐Channel Seismoelectric Spectral Ratios.pdf`：JGR Solid Earth 参考论文，用于论文定位、写作风格和震电应用背景参考。

## 文档与图示

- `doc/公式说明_中文版.md`：中文公式核对版。系统整理了从反应输运变量到动态系数、Schakel 界面边界值问题、Liu 频率-波数积分和偶极子解释模型的完整链条。
- `doc/公式说明.md`：英文/原始公式说明文档。
- `doc/代码验证标准.md`：当前代码审核与验证准则。重点包括因果性、T0 后出信号、界面上下极性反转、幅值随偏移距/空间位置变化，以及相位约定。
- `doc/震电结构布置图有偏移距.png`：有偏移距震电/VSEP 几何示意图。
- `doc/震电结构布置图无偏移距.png`：无偏移距震电/VSEP 几何示意图。
- `临时.md`：临时研究记录或待整理笔记。若其中内容被正式采用，应转移到 `doc/` 或论文草稿中。

## 结果目录

- `se_results_offset/`：当前主模型输出目录，包含有限偏移距谱积分波形、转换系数随溶蚀时间变化、峰值幅度趋势和参数记录。
  - `seismoelectric_timeseries_results.csv`：各溶蚀时刻的转换系数、峰值幅度等汇总结果。
  - `parameters_used.csv`：运行时使用的参数记录。
  - `run_summary.csv`：本次运行摘要。
  - `coefficients_vs_dissolution_time.png`：静态/界面转换系数随溶蚀时间变化。
  - `dynamic_coefficients_vs_dissolution_time.png`：动态电动参数随溶蚀时间变化。
  - `peak_amplitude_vs_dissolution_time.png`：峰值波形幅度随溶蚀时间变化。
  - `waveform_snapshot_spectral.csv/.npz/.png`：某一溶蚀时刻的频率-波数谱积分波形快照。
  - `waveform_spatial_peak_diagnostics.csv/.png`：空间峰值诊断结果。
- `se_results_offset copy/`：历史输出副本，文件结构与 `se_results_offset/` 类似。除非需要对比旧结果，否则不建议作为当前论文结果来源。
- `liu2018_fig2b_comparison/`：Liu Fig. 2(b) 对照实验输出，包含 Pride 理论波形、电偶极子空间分布、极性和归一化幅值图。
- `offset_zr_peak_comparison/`：不同偏移距下 `z/r` 或空间峰值趋势比较结果，用于检查有限偏移距几何下的峰值变化。

## 测试与缓存

- `tests/test_waveform_diagnostics.py`：主谱积分模型的单元测试，检查复波数分支、相位约定、因果 Ricker 频谱、Liu 电势系数转换等关键逻辑。
- `tests/test_liu2018_fig2b_dipole_comparison.py`：Liu Fig. 2(b) 对照脚本测试，检查默认积分采样数、电偶极子几何和极坐标图设置。
- `__pycache__/`、`tests/__pycache__/`：Python 缓存目录，无需人工维护。
- `.vscode/`：本地编辑器配置。
- `.git/`：Git 版本控制目录。

## 运行参考

主模型示例：

```bash
python seismoelectric_offset_liu2018_spectral.py --input global_evolution.xlsx --outdir se_results_offset
```

Liu Fig. 2(b) 对照示例：

```bash
python liu2018_fig2b_dipole_comparison.py --input global_evolution.xlsx --outdir liu2018_fig2b_comparison
```

测试示例：

```bash
python -m unittest discover -s tests
```

## 关键注意事项

- 相位约定需与 Schakel (2010) 和 Liu (2018) 保持一致：时间项采用 `exp(i omega t)`，空间传播采用 `exp(-i k dot x)`，因果声源谱使用 `exp(-i omega tau)` 傅里叶核。
- 复波数分支应保证 `exp(-i k z)` 沿传播方向衰减。
- 有偏移距波形应由频率-波数积分自然产生，不能为了贴合预期趋势而人为修改结果。
- interface EM signal 应在界面入射时刻 T0 之后自然出现；若 T0 前存在明显信号，需要优先检查相位、频谱窗、因果源和数值积分设置。
- 论文分析时建议同时报告：反应输运参数演化、动态电动参数演化、界面转换系数演化、波形峰值幅度演化、空间峰值/极性分布。
- 你在读论文的时候核对公式的时候，我建议截图读论文，这样准确率会高一点
