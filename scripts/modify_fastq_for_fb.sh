for f in *001.fastq.gz; do
  gzcat $f | sed '1,4d' | gzip > ${f//FB/CH}
done
rm *_FB_IGO*