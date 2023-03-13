`(jupyter) cfan@stout:~/genomes/NOD_Ly49fix$ nohup STAR --runThreadN 10 --runMode genomeGenerate --genomeDir STAR_noAnnot --genomeFastaFiles NOD_Ly49fix.fa 1>star.gen.log 2>&1 &`

2022-12-01: back up `NOD_Ly49fix.fa` using mv NOD_Ly49fix.fa bk.NOD_Ly49fix.fa`, before generating a new version with the Ly49 region with repeats in lower cases.
```
mv Ly49_NOD.1.fna bk.Ly49_NOD.1.fna
```
