for f in *001.fastq.gz; do
  zcat $f | sed '1,4d' | gzip > ${f//FB_IGO/CH_IGO}
done
rm *_FB_IGO*