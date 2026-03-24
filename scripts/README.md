# Scripts README

本目录包含用于数据检查与预处理的辅助脚本。

## 运行环境

- Python 3.8+（仅标准库，无第三方依赖）
- 在仓库根目录执行命令示例中的 `python3 scripts/...`

## 脚本说明

### `compare_starts.py`

用途：
- 对比原始 C 实现 TXT 输出与 Rust CSV 输出中的起始坐标是否一致。
- TXT 第 2 列是 0-based，脚本会自动转换为 1-based 后再比较。

输入：
- `txt`：C 实现输出的 TXT 文件
- `csv`：Rust 生成的 CSV 文件

参数：
- `--limit`：最多打印多少条 mismatch（默认 `20`）

示例：

```bash
python3 scripts/compare_starts.py ./old_impl.txt ./rust_impl.csv
python3 scripts/compare_starts.py ./old_impl.txt ./rust_impl.csv --limit 50
```

输出：
- 打印两侧总行数、mismatch 数量。
- 对每条 mismatch 打印当前行和前后邻居行，方便定位偏移或分组差异。

---

### `ncbi_region_split.py`

用途：
- 批量读取 NCBI assembly 下的 `*_genomic.gff(.gz)` 与 `*_genomic.fna(.gz)`。
- 提取 GFF 中 `feature=region` 的记录导出 CSV。
- 按 `seqid` 分割 FNA，输出 `seqid.fna`。

输入目录约定：
- `--input-root` 下每个子目录视为一个 assembly。
- 每个 assembly 目录中应各有且仅有 1 个：
  - `*_genomic.gff` 或 `*_genomic.gff.gz`
  - `*_genomic.fna` 或 `*_genomic.fna.gz`

参数：
- `--input-root`：assembly 根目录（必填）
- `--output-root`：输出根目录（必填）
- `--on-mismatch`：`warn_skip`（默认）或 `strict`
- `--summary-name`：总汇总 CSV 文件名（默认 `regions.all.csv`）

示例：

```bash
python3 scripts/ncbi_region_split.py \
  --input-root /path/to/refseq \
  --output-root /path/to/output
```

严格模式（遇到缺失/长度不一致返回非 0）：

```bash
python3 scripts/ncbi_region_split.py \
  --input-root /path/to/refseq \
  --output-root /path/to/output \
  --on-mismatch strict
```

输出：
- `output-root/<assembly_id>/regions.csv`
- `output-root/<assembly_id>/fna_split/<seqid>.fna`
- `output-root/<summary-name>`

`regions.csv`/汇总表主要字段：
- `assembly_id, seqid, region_start, region_end, region_length, source, strand`
- `genome, mol_type, is_circular, name, strain, plasmid_name, taxon_id`
- `fna_header, fna_length, length_match, status, attributes_json`

`status` 含义：
- `ok`：FNA 中存在该 `seqid` 且长度匹配，已分割输出
- `missing_in_fna`：GFF 中有但 FNA 中没有该序列
- `length_mismatch`：`region_length` 与 FNA 实际长度不一致
