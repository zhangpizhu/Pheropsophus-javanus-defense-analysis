#!/bin/bash
PROTEOME=${1:-PJA.fasta}   # P. javanus predicted protein set (FASTA)
# Extract significant CAFE families and member gene lists/sequences.
# Usage: bash extract_significant.sh [proteome.fasta]   (default: PJA.fasta in this folder)

#提取显著扩张或收缩的orthogroupsID
cat Base_family_results.txt |grep "y"|cut -f1 > p0.05.significant 
#显著扩张/收缩的基因家族在每个节点的收缩与扩增数目
grep -f p0.05.significant Base_change.tab > Base_p0.05change.tab
#PJA显著扩张的orthogroupsID
cat Base_p0.05change.tab |awk '$15 > 0 {print $1}' > PJA_significant.expand
cat Base_p0.05change.tab | cut -f1,15 | grep "-" | cut -f1 > PJA_significant.contracted

#提取PJA显著扩张的基因
grep -f PJA_significant.expand Orthogroups.txt |sed "s/ /\n/g" | grep "PJA" | sort | uniq  > PJA_significant.expand.genes
grep -f PJA_significant.contracted Orthogroups.txt |sed "s/ /\n/g" | grep "PJA" | sort | uniq  > PJA_significant.contracted.genes

## 提取显著扩张的基因序列
seqkit grep -f PJA_significant.expand.genes "${PROTEOME:-PJA.fasta}" >PJA_significant.expand.pep.fa
## 提取显著扩张的基因序列
seqkit grep -f PJA_significant.contracted.genes "${PROTEOME:-PJA.fasta}" >PJA_significant.contracted.pep.fa
