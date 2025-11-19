import sys
import tabix
import numpy as np
import re
import functools
import subprocess
from multiprocessing.dummy import Pool as ThreadPool

def load_chrlen(fai):
    M_chr_len = {}
    with open(fai, "r") as f_fai:
        for line in f_fai:
            line = line.strip()
            f = line.split()
            M_chr_len[f[0]] = f[1]
    return M_chr_len

fai = "/data/Analysis/huboqiang/Database_Meth/mm9/mm9_lambda.fa.fai"
M_chrlen = load_chrlen(fai)

def parse_region(region):
    chrom = re.split(":|-",region)[0]
    beg = int(re.split(":|-",region)[1])-1
    end = int(re.split(":|-",region)[2])+1
    return chrom, beg, end

func_cmp  = lambda x,TF: 1 if x.upper() == TF.upper() else 0
#pool = ThreadPool(8)

def parse_record(record, TF):
    l_motif = []
    M_motif_lv = {}
    for rec in record:
        l_info = rec[3].split("_")
        is_TF = sum(map(functools.partial(func_cmp, TF), l_info))
        if is_TF:
            if l_info[-1] not in M_motif_lv:
                M_motif_lv[l_info[-1]] = []
            
            M_motif_lv[l_info[-1]].append(int(rec[1]))
    
    for lv in ['t5000', 't20000', 't100000', 't250000']:
        if lv in M_motif_lv:
            l_motif = M_motif_lv[lv]
            break
    return l_motif

ucsc_dir = "/data/Analysis/huboqiang/software/UCSC"
def get_bwInfo(chrom, beg_out, end_out):
    shell_info = "%s/bigWigSummary %s %s %d %d %d" % (ucsc_dir, bigwigfile, chrom, beg_out, end_out, cnt_bin)
    p=subprocess.Popen(shell_info, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    if p.stdout:
        line = p.stdout.readline()
        line = line.strip()
        if line[0:2].upper() != "NO":
            line = line.replace('n/a', 'nan')
            f = line.split()
            chrpos = "%s:%d-%d" % (chrom, beg_out, end_out)
            out = "%s\t%s" % (chrpos, line)
            return out

def parse_line(line):
    line = line.strip()
        
    beg = int(f[1]) - 1
    end = int(f[2]) + 1
    try:
        record = tb_db.query(chrom, beg, end)
        l_motif = parse_record(record, TF)
        if len(l_motif) > 0:
            dist_ChIP_center = abs(np.array(l_motif) - (beg+end)/2)
            motif_pos = l_motif[dist_ChIP_center.argsort()[0]]
#            print l_motif, dist_ChIP_center, motif_pos
            beg_out = motif_pos - ext
            end_out = motif_pos + ext
            return get_bwInfo(chrom, beg_out, end_out)
    except:
        return None
        
    
infile = sys.argv[1]    #/data/Analysis/huboqiang/tmp/lilin_20151026/merge_mES_150913_col_1234.CTCF
motif_db = sys.argv[2] # /data/Analysis/huboqiang/Database_Meth/motif/mm9_motifs_split.bed.gz
bigwigfile = sys.argv[3] # /datd/lilin/Project/NOM-Seq/mouse_NOME_seq/02.SingleC/merge_mES_150913_col_1234/singleC/all.GCA.GCC.GCT.bw
TF = sys.argv[4]       # CTCF
ext = 1000
cnt_bin = 100

tb_db = tabix.open(motif_db)
with open(infile, "r") as f_infile:
    for line in f_infile:
        f = line.split()
        chrom = f[0]
        if not chrom in M_chrlen:
            continue
        
        out = parse_line(line)
        if out != None:
            print out

#    while 1:
#        l_lines = f_infile.readlines(10000)
#        if len(l_lines) == 0:
#            break
#        l_out = map(parse_line, l_lines)
#        print "\n".join(filter(lambda x: 0 if x == None else 1, l_out))
        
#pool.close()
#pool.join()
