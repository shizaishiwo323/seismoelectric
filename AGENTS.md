# 震电效应模拟 (Seismoelectric Simulation) - 项目说明

本项目主要用于基于孔隙尺度的反应输运（Reactive Transport）数据，模拟流体-多孔介质界面产生的震电效应（Seismoelectric Effect）。模型主要基于 **Schakel & Smeulders (2010)** 的边界条件架构以及 **Liu et al. (2018)** 风格的频率-波数谱积分机制来合成波形。

## 领域背景与核心理论

- **二元时间尺度 (极重要)**：项目涉及两个完全独立的时间维度：
  - **反应输运时间 $t_d$**（较长周期）：该随着岩石溶解，影响岩石物性（孔隙率、渗透率、流体浓度等）。
  - **声波/电磁波形传播时间 $t$**（微秒级周期 $t$）：在特定的 $t_d$ 截面的岩石属性下，声波激发的波形变化时间。
- **理论公式与符号约定**：所有的数学模型基础、中间过渡变量（如反射系数 $R_E$、透射系数 $T_{TM}$、各种频率响应函数）以及单位换算（如渗透率从 $mD$ 转换至 $m^2$ 等），**必须严格参考** [doc/公式说明.md](doc/公式说明.md) 。

## 项目结构与开发约定

- **核心运行脚本**：项目主要通过命令行参数 (CLI) 执行。
  - 主力频-谱积分版脚本：[seismoelectric_offset_liu2018_spectral.py](seismoelectric_offset_liu2018_spectral.py)
  - 其他基准脚本：[seismoelectric_nooffset.py](seismoelectric_nooffset.py)、[seismoelectric_offset.py](seismoelectric_offset.py)
- **依赖栈**：主要作为数据处理及科学计算模块，强依赖 `numpy`，`pandas`，`scipy` 及 `matplotlib`。
- **运行命令示例（带参数）**：
  ```bash
  python seismoelectric_offset_liu2018_spectral.py --input uploaded_files/<data>.xlsx
  ```

## AI 助手开发指导

1. **查阅文档先于推导**：当被要求修改震电转换系数、动电耦合或是修改波形绘制逻辑时，第一时间参考 [doc/公式说明.md](doc/公式说明.md) 以确保物理公式的准确性。
2. **保持命令行接口兼容**：在改进或增加如绘图选项、自定义变量输入等新功能时，向现有的 `argparse` 中添加可选参数，不要硬编码数据文件路径。
3. **因果性（Causality）约束**：确保在构建声源和合成接收波形时保持系统因果性（如在理论到达时间 T0 之前的界面 EM 响应要严格置零）。