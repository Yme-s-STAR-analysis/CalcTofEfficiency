# TOF match efficiency

`version: 4.1`

`author: yghuang`

## Usage

noooop

## Change log

05.04.2024 by yghuang (v4.1):

> 使用最新的组件

2023 Oct. 15 by yghuang (4.0):

> 更新了Centrality util.
>
> vz范围从30提升到了50.
>
> 合并不再保存为TH2D了，而是保留TEfficiency，然后后续用其他包保存为参数化的形式（拟合）. 

2023 Aug. 15 by yghuang (3.0):

> 现在用了新的StCFMult，CentralityUtil，TpcShiftUtil等工具。
>
> 现在区分不同的vz范围计算效率。

2023 Jan. 14 by yghuang (2.0):

> 增加二维效率。
>
> 计算系统误差时，现在在readPicoDst.C里面修改cut，不需要额外改代码并cons。

2022 Dec. 26 by yghuang (2.0):

> 根据新的规范进行重写。

2021 September 22 by yghuang (1.2):

> Update: The effCalc part is finished.

2021 September 22 by yghuang (1.1):

> Update: Use TH1D (only pT) instead of TH3D (y, phi and pT).

2021 September 21 by yghuang (1.0):

> A release version.

