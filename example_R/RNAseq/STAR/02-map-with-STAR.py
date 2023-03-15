## Map cleaned reads against Long-read assembly + newly released annotation
import os
from glob import glob

outf = open("02-map-with-STAR.sh", 'w')

cmd = "STAR --genomeLoad LoadAndExit --genomeDir ../Reference/staridx\n"
outf.write(cmd)

for f in glob("../01-cleaned/*.gz"):
    s = f.split('/')[-1].replace('_SE.fastq.gz', '')
    print(s)
    cmd = "STAR --genomeDir ../Reference/staridx"
    cmd += " --runThreadN 30 --readFilesIn " +  f 
    cmd += " --outFileNamePrefix " + s
    cmd += " --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outBAMsortingBinsN 50"
    cmd += " --readFilesCommand zcat --quantMode GeneCounts --outSAMattributes Standard\n"
    #os.system(cmd)
    outf.write(cmd)

cmd = "STAR --genomeDir ../Reference/staridx --genomeLoad Remove\n"
outf.write(cmd)
outf.close()
