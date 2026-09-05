# dsh-chem

[DeepSeek Harness](https://github.com/deepseek-ai/deepseek-harness)（`dsh`）的化学计量数学工具箱——**零运行时依赖**，基于 IUPAC 数据的纯算术计算。

专门替代 LLM 高频出错的化学心算：摩尔质量加总（下标数错、漏掉结晶水）、元素符号大小写（Co 是钴、CO 是一氧化碳）、稀释换算（C1V1=C2V2）、pH 对数/反对数运算。

## 四个工具

| 工具 | 功能 | 示例 |
|------|------|------|
| `formula_mass` | 化学式解析 + 摩尔质量 + 元素明细（支持嵌套括号、结晶水 `CuSO4·5H2O`、离子电荷 `Fe3+`/`SO4^2-`、前置系数 `2H2O`） | `"H2SO4"` → 98.072 g/mol |
| `mass_moles` | 质量(g) ↔ 物质的量(mol) ↔ 粒子数（阿伏伽德罗常数 6.02214076e23），六种模式 | 98 g H2SO4 → 0.99927 mol |
| `solution` | 摩尔浓度 / 配制称量 / 稀释 C1V1=C2V2 / 百分比浓度（w/v、v/v、质量%） | 2 M × 62.5 mL = 0.5 M × 250 mL |
| `ph` | pH/pOH/[H+]/[OH-] 换算 + 缓冲溶液（Henderson–Hasselbalch），25 °C、Kw=1e-14 | pKa 4.76、[A-]/[HA]=0.5 → pH 4.459 |

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
pnpm test      # vitest，102 测试
pnpm lint      # oxlint src test
```

数值锚点全部由独立 Python oracle（`test/anchors-chem.py`，自建元素表+自建分词器）生成，并与公开教科书值交叉验证（H2O 18.015、H2SO4 98.08、NaCl 58.44、Ca(OH)2 74.09、CuSO4·5H2O 249.68、C6H12O6 180.16、KMnO4 158.03、NH3 17.03）。

## 许可证

MIT
