# Chromatin VUS Workflow

输入突变列表，获得可追溯的多数据库注释、图表和后续实验建议。流程面向人类表观修饰与染色质调控因子的功能研究。

它把论文中的“位点筛选 → 内源性编辑 → 表观修饰检测 → 功能表型验证”整理成可复用的研究入口。同一结构域的不同替换可能产生不同效应，因此工具保留每个位点的证据，而不把多个预测分数合并成一个“致病性总分”。

## 首版功能

- 支持 hg38/GRCh38 的 VCF、染色体位置加 REF/ALT 的 CSV/TSV，以及带版本号的转录本编码区 HGVS。
- 自动查询 REVEL、AlphaMissense、SIFT、PolyPhen、CADD，以及 phyloP、GERP++、phastCons 保守性指标；补充 ClinVar 与可取得的 gnomAD 频率。
- 输出可筛选的 HTML 报告、CSV/JSON 表格、PNG/SVG 图和逐位点实验建议。
- 保留查询时间、版本信息、参数、原始响应与校验值；可用缓存离线重放。
- 内置 14 个基因的起始面板，也可以提供自己审核的面板。首版只分析错义 SNV，其他类型明确提示。

## 使用

需要 Python 3.10 及以上。下载仓库后进入目录，执行：

```bash
python -m venv .venv
```

Windows PowerShell 激活环境：`.venv\Scripts\Activate.ps1`；macOS/Linux：`source .venv/bin/activate`。

```bash
python -m pip install -e .
chromatin-vus demo --out outputs/demo
chromatin-vus annotate examples/kmt2c_hgvs.csv --out outputs/kmt2c
```

打开输出目录里的 `report.html` 即可查看。`demo` 完全使用虚构数值，用于演示版式，不需要联网；KMT2C 示例联网查询公开变异，不包含课题组未发表的实验数据。

需要上传文件的本地界面时：

```bash
python -m pip install -e ".[ui]"
streamlit run app.py --server.address 127.0.0.1
```

## 输入示例

```csv
input_id,chrom,pos,ref,alt,gene,transcript,assembly
V4834A,7,152139219,A,G,KMT2C,NM_170606.3,GRCh38
```

也可使用：

```csv
input_id,gene,hgvs
V4834A,KMT2C,NM_170606.3:c.14501T>C
```

基因组坐标为 1-based、正链 REF/ALT。KMT2C 位于负链，所以转录本 `T>C` 对应基因组 `A>G`。流程核查这一转换并核对转录本版本。只有蛋白替换如 `V4834A` 不足以唯一确定基因组等位基因，首版要求补齐上述输入之一。

VCF 多个 ALT 分开处理，不解读样本基因型，也不按 FILTER 列筛选质量；请预先完成质控。暂不读取 `.vcf.gz`，也不自动转换 hg19 坐标。

## 如何理解结果

REVEL 等是模型预测，保守性指标反映序列约束，ClinVar 是数据库提交的临床解释，实验测量则是另一层证据。报告分别展示这些内容。缺失评分保留为 NA；缺失人群频率不当作 0。不同模型彼此相关，不能用“多数投票”累加为独立致病证据。

默认 REVEL ≥0.5 和 AF <0.001 是可修改的研究筛选参数，不是经过本项目临床验证的判断标准。AF 字段取返回的外显子组和基因组全体频率的较大值，不是祖源亚群 popmax。图中的 SIFT 越低代表预测影响越大，其他三个 0–1 模型通常相反。颜色表达原始数值，不代表统一的风险等级。

实验建议由公开规则生成，并给出触发原因及来源。例如 KMT2C 推荐以 H3K4me1 为起点，结合蛋白表达、染色质分布与细胞表型；DNA 修饰因子和染色质重塑因子会收到不同建议。建议不能替代研究者对基因机制和细胞模型的判断。

## 可复现与扩展

- 一次运行生成注释表、实验建议表、图、`provenance.json` 和 `raw/` 响应，保留完整目录便于复核。
- `--offline --cache 先前输出/raw` 可离线重放同一批查询；换用编码 HGVS 与基因组 HGVS 会形成不同请求，需先缓存对应请求。
- 复制 `src/chromatin_vus/data/panel.csv` 并修改，以 `--panel 自定义.csv` 替换内置面板。
- 每次指定新的输出目录，避免覆盖旧分析。默认最多 200 个等位基因。
- 方法边界、数据来源及引用见 [METHODS](docs/METHODS.md) 和 [DATA_SOURCES](docs/DATA_SOURCES.md)。

在线模式把变异查询发送给 Ensembl 和 MyVariant.info。代码采用 MIT 许可，但上游数据库和预测分数各有使用条款，不能把本项目的代码许可套用于第三方数据。
