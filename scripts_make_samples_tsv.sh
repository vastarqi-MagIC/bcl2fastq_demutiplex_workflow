#!/usr/bin/env bash
set -euo pipefail

FASTQ_DIR="/mnt/nfs1_external/rscc/vastar/demultiplex/result_s0819a0901"
OUT="samples.tsv"

# 常见命名模式：SampleID_S1_L001_R1_001.fastq.gz
# 如果你的命名不同，请把 sed 那行贴给我再改。
echo -e "sample_id\tfq1\tfq2" > "$OUT"

# 用关联数组避免同一 sample 重复写入
declare -A SEEN

shopt -s nullglob
for r1 in "${FASTQ_DIR}"/*_R1_*.fastq.gz; do
  base=$(basename "$r1")
  # 提取 sample_id
  sample=$(echo "$base" | sed -E 's/_S[0-9]+_L[0-9]+_R1_001\.fastq\.gz$//')
  # 如果没有 _Sxx_Lxxx 结构，可尝试： sed 's/_R1_001\.fastq\.gz$//'
  [[ -z "$sample" ]] && { echo "WARN: 无法解析样本名: $base" >&2; continue; }

  # 找对应 R2
  r2="${r1/_R1_/_R2_}"
  if [[ ! -f "$r2" ]]; then
    echo "ERROR: 未找到匹配的 R2: $r2 (对应 $r1)" >&2
    exit 1
  fi

  # 只写一次
  if [[ -z "${SEEN[$sample]:-}" ]]; then
    echo -e "${sample}\t${r1}\t${r2}" >> "$OUT"
    SEEN[$sample]=1
  fi
done

echo "生成完成: $OUT"
wc -l "$OUT"