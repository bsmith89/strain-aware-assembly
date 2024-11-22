gff=$1
depth=$2

# Setup data streams
## Gene
gene_info_pipe=$$.gene.pipe
mkfifo $gene_info_pipe
cat ${gff} | awk -v OFS='\t' '$0!~/^#/{print $1,$4 - 1,$5,$7}' > $gene_info_pipe &

## Depth
depth_pipe=$$.depth.pipe
mkfifo $depth_pipe
cat ${depth} | tqdm --unit-scale > $depth_pipe &

# Run SQL
cat <<EOF | sqlite3 -separator '	' :memory:

CREATE TABLE gene (
    seqid TEXT
  , left INT
  , right INT
  , strand VARCHAR(1)
  , PRIMARY KEY (seqid, left, right, strand)
);

CREATE TABLE depth (
    seqid TEXT
  , pos INT
  , depth INT
  , PRIMARY KEY (seqid, pos)
);

.import $gene_info_pipe gene
.import $depth_pipe depth

SELECT
    seqid || "[" || left || "-" || right || "]" || strand AS gene_name
  , (1.0 * SUM(depth)) / (right - left) AS depth
FROM depth
JOIN gene USING (seqid)
WHERE (depth.pos > gene.left) AND (depth.pos < gene.right)
GROUP BY seqid, left, right, strand
EOF

wait
rm $gene_info_pipe $depth_pipe
