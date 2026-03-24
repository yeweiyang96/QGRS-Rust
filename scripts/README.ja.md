# Scripts README

このディレクトリには、データ確認と前処理のための補助スクリプトがあります。

## 動作環境

- Python 3.8+（標準ライブラリのみ、外部依存なし）
- リポジトリのルートで `python3 scripts/...` 形式で実行

## スクリプト一覧

### `compare_starts.py`

目的:
- 旧 C 実装の TXT 出力と Rust 実装の CSV 出力で、開始座標が一致するかを比較します。
- TXT の第2列は 0-based のため、比較前に 1-based へ変換します。

入力:
- `txt`: C 実装が出力した TXT ファイル
- `csv`: Rust 実装が出力した CSV ファイル

引数:
- `--limit`: 表示する不一致件数の上限（デフォルト `20`）

例:

```bash
python3 scripts/compare_starts.py ./old_impl.txt ./rust_impl.csv
python3 scripts/compare_starts.py ./old_impl.txt ./rust_impl.csv --limit 50
```

出力:
- 両ファイルの行数と不一致件数を表示します。
- 各不一致について、該当行と前後の隣接行を表示し、ズレやグルーピング差分の特定を助けます。

---

### `ncbi_region_split.py`

目的:
- NCBI の `*_genomic.gff(.gz)` と `*_genomic.fna(.gz)` をバッチ処理します。
- GFF から `feature=region` を抽出して CSV に出力します。
- FNA を `seqid` ごとに分割し、`seqid.fna` を出力します。

入力ディレクトリ規約:
- `--input-root` 配下の1階層目サブディレクトリを 1 assembly とみなします。
- 各 assembly ディレクトリには、次のファイルがそれぞれ1つ必要です:
  - `*_genomic.gff` または `*_genomic.gff.gz`
  - `*_genomic.fna` または `*_genomic.fna.gz`

引数:
- `--input-root`: assembly ルートディレクトリ（必須）
- `--output-root`: 出力ルートディレクトリ（必須）
- `--on-mismatch`: `warn_skip`（デフォルト）または `strict`
- `--summary-name`: 集計 CSV のファイル名（デフォルト `regions.all.csv`）

例:

```bash
python3 scripts/ncbi_region_split.py \
  --input-root /path/to/refseq \
  --output-root /path/to/output
```

Strict モード（欠損/長さ不一致で非0終了）:

```bash
python3 scripts/ncbi_region_split.py \
  --input-root /path/to/refseq \
  --output-root /path/to/output \
  --on-mismatch strict
```

出力:
- `output-root/<assembly_id>/regions.csv`
- `output-root/<assembly_id>/fna_split/<seqid>.fna`
- `output-root/<summary-name>`

`regions.csv` / 集計 CSV の主な列:
- `assembly_id, seqid, region_start, region_end, region_length, source, strand`
- `genome, mol_type, is_circular, name, strain, plasmid_name, taxon_id`
- `fna_header, fna_length, length_match, status, attributes_json`

`status` の意味:
- `ok`: FNA に `seqid` が存在し長さ一致。分割ファイルを出力済み
- `missing_in_fna`: GFF にはあるが FNA には該当 `seqid` がない
- `length_mismatch`: `region_length` と FNA 実配列長が一致しない
