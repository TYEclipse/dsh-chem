# dsh-chem

[DeepSeek Harness](https://github.com/deepseek-ai/deepseek-harness)（`dsh`）的化学计量数学工具箱——**零运行时依赖**，基于 IUPAC 数据的纯算术计算。

专门替代 LLM 高频出错的化学心算：摩尔质量加总（下标数错、漏掉结晶水）、元素符号大小写（Co 是钴、CO 是一氧化碳）、稀释换算（C1V1=C2V2）、pH 对数/反对数运算。

## 七个工具

| 工具 | 功能 | 示例 |
|------|------|------|
| `formula_mass` | 化学式解析 + 摩尔质量 + 元素明细（支持嵌套括号、结晶水 `CuSO4·5H2O`、离子电荷 `Fe3+`/`SO4^2-`、前置系数 `2H2O`） | `"H2SO4"` → 98.072 g/mol |
| `mass_moles` | 质量(g) ↔ 物质的量(mol) ↔ 粒子数（阿伏伽德罗常数 6.02214076e23），六种模式 | 98 g H2SO4 → 0.99927 mol |
| `solution` | 摩尔浓度 / 配制称量 / 稀释 C1V1=C2V2 / 百分比浓度（w/v、v/v、质量%） | 2 M × 62.5 mL = 0.5 M × 250 mL |
| `ph` | pH/pOH/[H+]/[OH-] 换算 + 缓冲溶液（Henderson–Hasselbalch），25 °C、Kw=1e-14 | pKa 4.76、[A-]/[HA]=0.5 → pH 4.459 |
| `balance_equation` | 配平方程式：最小整数系数（全程 BigInt 有理数运算，不出 2.9999997）；含电荷守恒、结晶水、元素原子数双侧核对；配不平就明说 | `C3H8 + O2 -> CO2 + H2O` → `C3H8 + 5 O2 -> 3 CO2 + 4 H2O` |
| `stoichiometry` | 限量试剂计算：给出方程式与**每种**反应物的实际量（mol 或 g）→ 反应进度、哪种反应物先耗尽、各产物生成量与剩余量；给 `product` + `actualYieldG` 还可算理论产量与产率。每个反应物都必须给量（不给就是错，不假设“过量”）；恰好配比时**两种都算限量** | `C3H8` 0.5 mol + `O2` 2 mol → 限量 O2、生成 CO2 52.8108 g、剩 C3H8 4.4097 g |
| `gas_law` | 理想气体：`ideal`（给定三者解第四个）、`molar_volume`（标准状态 22.414 L/mol）、`density`（ρ=P·M/RT）、`molar_mass`（由密度反推 M）；固定单位 atm/L/mol/K，附带 kPa、mmHg、mL、°C 换算 | 1 mol、273.15 K、22.4 L → 1.00062 atm |

两条容易踩的确定性规则：**配平**里有电荷的物种必须写成带空格的形式（`Na+ + Cl-`），因为紧贴符号的 `+` 是电荷而不是分隔符；**气体**里只少给一个量才能解（`ideal` 必须恰好给三个），缺项与负值一律报错而不是猜。

## 电荷歧义的确定性规则

`Fe3+` 与 `NH4+` 语法上不可分（数字跟正负号），采用规则：

- **数字前只有同一种元素** → 数字是电荷量：`Fe3+` = Fe 带 3+、`S2-` = S 带 2-、`Ca2+` = Ca 带 2+
- **否则** → 数字是下标、裸符号是单位电荷：`NH4+` = N + H4 带 `+`、`NO3-` = N + O3 带 `-`
- **显式 `^` 形式恒为电荷**：`SO4^2-`、`O2^+`
- 单个小写字母符号自动纠正（`h2o` → H2O）；多字母小写（`nacl`）拒绝并提示大小写

## 安装

```bash
dsh plugin --profile <profile> add github:TYEclipse/dsh-chem
```

（PATH 需有 `pnpm`——`dsh plugin` 内部转发给它。）

## 验证

```bash
dsh plugin --profile <profile> list      # 出现 dsh-chem
dsh --dump-config --profile <profile>    # 出现 == dsh-chem 层
```

## 开发

```bash
pnpm install
pnpm build     # tsc -> dist/（已提交；git 安装不跑构建）
pnpm test      # vitest，136 测试
pnpm lint      # oxlint src test
```

数值锚点全部由独立 Python oracle（`test/anchors-chem.py`，自建元素表+自建分词器+自建 Fraction 配平器）生成，并与公开教科书值交叉验证（H2O 18.015、H2SO4 98.08、NaCl 58.44、Ca(OH)2 74.09、CuSO4·5H2O 249.68、C6H12O6 180.16、KMnO4 158.03、NH3 17.03；丙烷 1/5/3/4；高锰酸钾+盐酸 2/16/2/2/5/8；标准状态 22.414 L/mol）。运行：`python3 test/anchors-chem.py`。

## 许可证

MIT
