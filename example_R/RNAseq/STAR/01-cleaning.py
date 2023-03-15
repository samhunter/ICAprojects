from glob import glob
import os
import re

cleaning = open("01-cleaning_commands.sh", 'w')

# These are 3' QuantSeq libraries, and were only sequenced as SingleEnd 100bp reads
os.system('mkdir -p 01-cleaned')

for f in glob("./00-RawData/*.fastq.gz"):
    s = re.split('_S.*_L..._R', f.strip().split('/')[2])[0]
    print(s)
    log = f"./01-cleaned/{s}_stats.log"
    cmd = f"hts_Stats -L {log} -U {f} -N 'Raw read statistics' | " # record stats for raw reads
    cmd += f"hts_SeqScreener -N 'Screen PhiX' -k 12 -A {log} | "  # screen for phix  
    #cmd += "hts_SuperDeduper -e 100000 -S -O -A -L {log}  | " # remove PCR duplicates, record saturation
    cmd += f"hts_SeqScreener -r -s ./Ref/ixodes-scapularis_rRNA.fasta -N 'Count rRNA' -A {log}  | " # count rRNA reads
    cmd += f"hts_SeqScreener -N 'Screen Adapters' -A {log} -k 15 -x .01 --seq adapters.fa | "
    cmd += f"hts_NTrimmer -A {log} -N 'Trim Ns' | "
    cmd += f"hts_QWindowTrim -N 'Quality trim reads' -A {log} | "
    #cmd += "hts_PolyATTrim -l -m 75 -A {log}  "  # right trim
    #cmd += "hts_AdapterTrimmer -S -m 25 -A -L {log}  -O | " # trim adapters, throw out short fragements
    cmd += f"hts_LengthFilter -m 75 -N 'Discard reads less than 75bp' | "
    cmd += f"hts_Stats -A {log} -f ./01-cleaned/{s}" 
    #cmd += "hts_CutTrim -a 40 -m 90 -S -A -L " + log + " -O | "  # Cut off probe
    #cmd = "hts_NTrimmer -n -m 100 -L " + log + " -1 " + r1 + " -2 " + r2 + " -O |"
    #cmd += "hts_Overlapper -m 300 -n -S -o 15 -e 0.2 -A -L " + log + " -O |"
    #cmd += "hts_SeqScreener -g -S -A -L " + log + " -k 15 -x .01 --seq adapters.fa -f -p ./01-cleaned/" + s  # adapters
    cleaning.write(cmd+'\n')

cleaning.close()

