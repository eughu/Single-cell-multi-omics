#-*- coding:utf-8 -*-
from __future__ import division
import re
import sys
import os
import subprocess
import cPickle as pickle
import time

import MethGC.utils.module_running_jobs as m_jobs
import MethGC.settings.projpath         as m_proj

class Scripts(m_proj.ProjInfo):
    def __init__(self):
        super( Scripts,self ).__init__()
        
    def define_files(self, ref):
        dir_db = "%s/%s" % (self.Database, ref)
        self.ref           = ref
        self.genome_gtf    = "%s/refGene.gtf"           % (dir_db)
        self.intragenic_bed= "%s/region.Intragenic.bed" % (dir_db)
        self.rmsk_bed      = "%s/chrom.sort.bed"        % (dir_db)

    def db_01_DownloadRef(self):
        l_sh_info = []
        l_sh_info.append("ref=$1")
        l_sh_info.append("dir_database=%s/$ref" % (self.Database))
        l_sh_info.append("dir_path=%s"          % (self.path))
        l_sh_info.append("""
cd $dir_database

wget http://hgdownload.soe.ucsc.edu/goldenPath/${ref}/bigZips/chromFa.tar.gz

tar -zxvf $dir_database/chromFa.tar.gz

for i in {1..22} M X Y
do
    cat $dir_database/chr$i.fa
done  >$dir_database/${ref}_lambda.fa && rm $dir_database/chr*fa          &&\\

cat $dir_path/database/lambda.fa >>$dir_database/${ref}_lambda.fa
        """)
        return l_sh_info

    def db_02_BuildRefIndex(self):
        l_sh_info = []
        l_sh_info.append("ref=$1")
        l_sh_info.append("dir_database=%s/$ref" % (self.Database))
        l_sh_info.append("bowtie_dir=%s"   % (self.sftw_bowtie_dir))
        l_sh_info.append("bismark_dir=%s"  % (self.sftw_bismark ))
        l_sh_info.append("samtools_exe=%s" % (self.sftw_samtools))
        l_sh_info.append("""
$bismark_dir/bismark_genome_preparation                                     \\
    --path_to_bowtie $bowtie_dir --verbose  $dir_database

$samtools_exe faidx $dir_database/${ref}_lambda.fa
        """)
        return l_sh_info
    
    def db_03_RefGene(self):
        l_sh_info = []
        l_sh_info.append("ref=$1")
        l_sh_info.append("dir_database=%s/$ref" % (self.Database))
        l_sh_info.append("bedtools_exe=%s" % (self.sftw_bedtools))
        l_sh_info.append("ucsc_dir=%s"     % (self.sftw_ucsc_dir))
        l_sh_info.append("bin=%s"          % (self.bin))
        l_sh_info.append("dir_path=%s"     % (self.path))
        l_sh_info.append("""
cd $dir_database
wget http://hgdownload.soe.ucsc.edu/goldenPath/${ref}/database/refGene.txt.gz

### remove chromosome fragments(unassembled).
for i in {1..22} M X Y
do
    zcat $dir_database/refGene.txt.gz | grep -w chr$i
done >$dir_database/tmp
gzip $dir_database/tmp
mv $dir_database/tmp.gz $dir_database/refGene.txt.gz

# ref.sort.txt
# For CIRCExplorer(circular RNA pipeline)
zcat $dir_database/refGene.txt.gz                                          |\\
awk '{
    OFS="\\t";
    print $13,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11
}' /dev/stdin                                                              |\\
sort /dev/stdin >$dir_database/ref.sort.txt                               &&\\

# refGene.bed
zcat $dir_database/refGene.txt.gz                                          |\\
awk '{
    tag="noncoding";
    if($4~/^NM/){tag="protein_coding"};
    OFS="\\t";
    print $3,$5,$6,$2,$4,$10,$11,tag,$13
}' /dev/stdin                                                              |\\
python $bin/s03_genePred2bed.py /dev/stdin                                 |\\
$bedtools_exe sort -i /dev/stdin >$dir_database/refGene.bed               &&\\

# region.Intragenic.bed
# For novo lncRNA detection
$bin/find_ExonIntronIntergenic/find_ExonIntronIntergenic                    \\
    $dir_database/refGene.bed                                               \\
    $dir_database/${ref}_lambda.fa.fai >$dir_database/pos.bed             &&\\

grep -v "Intergenic" $dir_database/pos.bed                                 |\\
    awk '{OFS="\t";print $1,$2,$3,"Intragenic"}' /dev/stdin                 \\
    >$dir_database/region.Intragenic.bed                                  &&\\

# refGene.gtf
# For mapping
zcat $dir_database/refGene.txt.gz                                          |\\
cut -f 2-                                                                  |\\
$ucsc_dir/genePredToGtf file stdin /dev/stdout                             |\\
grep -w exon                                                               |\\
$bedtools_exe sort -i /dev/stdin >$dir_database/refGene.gtf
        """)
        return l_sh_info

    def db_04_rmsk(self):
        l_sh_info = []
        l_sh_info.append("ref=$1")
        l_sh_info.append("dir_database=%s/$ref" % (self.Database))
        l_sh_info.append("bedtools_exe=%s" % (self.sftw_bedtools))
        l_sh_info.append("ucsc_dir=%s"     % (self.sftw_ucsc_dir))
        l_sh_info.append("bin=%s"          % (self.bin))
        l_sh_info.append("dir_path=%s"     % (self.path))
        l_sh_info.append("""
mkdir -p $dir_database/ChromOut
cd $dir_database/ChromOut

wget -O  $dir_database/ChromOut/chromOut.tar.gz                             \\
    http://hgdownload.soe.ucsc.edu/goldenPath/${ref}/bigZips/chromOut.tar.gz

tar -zxvf $dir_database/ChromOut/chromOut.tar.gz
for i in {1..22} M X Y 
do
    tail -n +4 $i/chr$i.fa.out | \\
    awk '{
        if ($9=="C"){$9="-"};
        OFS="\t";print $5,$6,$7,$1,".",".",".",$8,$9,$10,$11,$12,$13,$14,$15
    }'
done >$dir_database/chrom.bed

bedtools sort -i $dir_database/chrom.bed >$dir_database/chrom.sort.bed
        """)
        return l_sh_info
    
    def db_05_Region(self, query = "GCA.GCC.GCT"):
        l_sh_info     = []
        l_sh_info.append("ref=$1")
        l_sh_info.append("query=$2")
        l_sh_info.append("enzymeSite_exe=%s/get_enzymeSites/get_enzymeSite" % self.bin)
        l_sh_info.append("dir_database=%s/$ref" % (self.Database))
        l_sh_info.append("bgzip_exe=%s" % (self.sftw_bgzip))
        l_sh_info.append("tabix_exe=%s" % (self.sftw_tabix))
        l_sh_info.append("""
$enzymeSite_exe -s $query $dir_database/${ref}_lambda.fa

mkdir -p $dir_database/${ref}_lambda.fa.$query

$bgzip_exe -fc $dir_database/${ref}_lambda.fa.$query.bed >$dir_database/${ref}_lambda.fa.$query.bed.gz
$tabix_exe -f -p bed -s 1 -b 2 -e 3 $dir_database/${ref}_lambda.fa.$query.bed.gz

for i in `cut -f 1 $dir_database/${ref}_lambda.fa.fai`
do
    tabix $dir_database/${ref}_lambda.fa.$query.bed.gz $i >$dir_database/${ref}_lambda.fa.$query/$i.bed
done
        """)
        return l_sh_info





    def s_01_trim(self):
        l_sh_info     = []
        l_sh_info.append("in_fq1=$1")
        l_sh_info.append("in_fq2=$2")
        l_sh_info.append("meth_name=$3")
        l_sh_info.append("trim_dir=%s" % (self.dir_trim_data))
        l_sh_info.append("trim_galore=%s"  % (self.sftw_trim))
        l_sh_info.append("""
trim_sam_dir=$trim_dir/$meth_name

$trim_galore                                                                \\
    --quality 20 --stringency 3 --length 50 --clip_R1 9 --clip_R2 9         \\
    --paired --trim1 --phred33 --gzip                                       \\
    --output_dir $trim_sam_dir                                              \\
    $in_fq1 $in_fq2
        """)
        return l_sh_info

    
    def s_02_bismark(self):
        l_sh_info = []
        l_sh_info.append("meth_name=$1")
#        l_sh_info.append("trim_fq1=$2")
#        l_sh_info.append("trim_fq2=$3")
        l_sh_info.append("ref=$2")
        l_sh_info.append("perl_exe=%s" % (self.sftw_pl))
        l_sh_info.append("bismark_exe=%s/bismark.origin" % (self.sftw_bismark))
        l_sh_info.append("bowtie_dir=%s" % (self.sftw_bowtie_dir))
        l_sh_info.append("changeID_pl=%s/ChangeReadID.pl" % (self.bin))
        l_sh_info.append("trim_fq_dir=%s" % (self.dir_trim_data))
        l_sh_info.append("bam_dir=%s"     % (self.dir_bismark))
        l_sh_info.append("database=%s/$ref" % (self.Database))
        l_sh_info.append("samtools_exe=%s" % (self.sftw_samtools))
        l_sh_info.append("""
trim_fq1=${meth_name}_R1_val_1.fq.gz
trim_fq2=${meth_name}_R2_val_2.fq.gz

trim_sam_dir=$trim_fq_dir/$meth_name
bam_sam_dir=$bam_dir/$meth_name

sam=$bam_sam_dir/${trim_fq1}_bismark_pe.sam
unmap_fq1=$bam_sam_dir/${trim_fq1}_unmapped_reads_1.txt
unmap_fq2=$bam_sam_dir/${trim_fq2}_unmapped_reads_2.txt

sam1_unmap=$bam_sam_dir/unmap1/${trim_fq1}_unmapped_reads_1.txt_bismark
sam2_unmap=$bam_sam_dir/unmap2/${trim_fq2}_unmapped_reads_2.txt_bismark


$bismark_exe      --fastq    --non_directional   --unmapped                 \\
    --phred33-quals --path_to_bowtie $bowtie_dir                            \\
    --output_dir $bam_sam_dir --temp_dir $bam_sam_dir $database             \\
    -1 $trim_sam_dir/$trim_fq1 -2 $trim_sam_dir/$trim_fq2                && \\
    $samtools_exe view -u -b -S -t $database/${ref}_lambda.fa $sam         |\\
    $samtools_exe sort -m 1000000000 - $sam.sort


$bismark_exe      --fastq    --non_directional   --unmapped                 \\
    --phred33-quals --path_to_bowtie $bowtie_dir                            \\
    --output_dir $bam_sam_dir/unmap1 --temp_dir $bam_sam_dir/unmap1         \\
    $database   $unmap_fq1                                               && \\
    $samtools_exe view -uSb -t $database/${ref}_lambda.fa $sam1_unmap.sam  |\\
    $samtools_exe sort -m 1000000000 - $sam1_unmap.sort

$bismark_exe      --fastq    --non_directional   --unmapped                 \\
    --phred33-quals --path_to_bowtie $bowtie_dir                            \\
    --output_dir $bam_sam_dir/unmap2 --temp_dir $bam_sam_dir/unmap2         \\
    $database   $unmap_fq2                                               && \\
    $samtools_exe view -uSb -t $database/${ref}_lambda.fa $sam2_unmap.sam  |\\
    $samtools_exe sort -m 1000000000 - $sam2_unmap.sort

$perl_exe $changeID_pl $sam.sort.bam       $sam.sort.ReID.bam            && \\
$samtools_exe rmdup    $sam.sort.ReID.bam                                   \\
                       $sam.sort.ReID.rmdup.bam                          && \\

$perl_exe $changeID_pl $sam1_unmap.sort.bam $sam1_unmap.sort.ReID.bam    && \\
$samtools_exe rmdup -s $sam1_unmap.sort.ReID.bam                            \\
                       $sam1_unmap.sort.ReID.rmdup.bam                   && \\

$perl_exe $changeID_pl $sam2_unmap.sort.bam $sam2_unmap.sort.ReID.bam    && \\
$samtools_exe rmdup -s $sam2_unmap.sort.ReID.bam                            \\
                       $sam2_unmap.sort.ReID.rmdup.bam                   && \\

$samtools_exe merge -f $bam_sam_dir/$meth_name.rmdup.bam                    \\
                       $sam.sort.ReID.rmdup.bam                             \\
                       $sam1_unmap.sort.ReID.rmdup.bam                      \\
                       $sam2_unmap.sort.ReID.rmdup.bam                   && \\
    
$samtools_exe sort -m 1000000000  $bam_sam_dir/$meth_name.rmdup.bam          \\
    $bam_sam_dir/$meth_name.sort.rmdup                                   && \\

$samtools_exe index $bam_sam_dir/$meth_name.sort.rmdup.bam


#rm  $sam1_unmap.sort.bam      $sam2_unmap.sort.bam                          \\
#    $sam1_unmap.sort.ReID.bam $sam2_unmap.sort.ReID.bam                     \\
#    $bam_sam_dir/*reads*txt   $bam_sam_dir/$meth_name.rmdup.bam             \\
#    $sam  $sam.sort.bam  $sam.sort.ReID.bam  $sam.sort.ReID.rmdup.bam       \\
#    $bam_sam_dir/unmap1/*fq.gz_unmapped_reads_1.txt_unmapped_reads.txt      \\
#    $bam_sam_dir/unmap2/*fq.gz_unmapped_reads_2.txt_unmapped_reads.txt      \\
#    $bam_sam_dir/unmap1/*.gz_unmapped_reads_1.txt_bismark.sort.ReID.rmdup.bam \\
#    $bam_sam_dir/unmap2/*.fq.gz_unmapped_reads_2.txt_bismark.sort.ReID.rmdup.bam \\
#    ${sam1_unmap}*sam ${sam2_unmap}*sam                     
#rm $trim_sam_dir/$trim_fq1 $trim_sam_dir/$trim_fq2

rm $bam_sam_dir/*reads*txt $bam_sam_dir/unmap1/*fq.gz_unmapped_reads_1.txt_unmapped_reads.txt $bam_sam_dir/unmap2/*fq.gz_unmapped_reads_2.txt_unmapped_reads.txt $sam ${sam1_unmap}*sam ${sam2_unmap}*sam

        """)
        return l_sh_info



    def s_03_Bam2SingleC(self):
        l_sh_info = []
        l_sh_info.append("meth_name=$1")
        l_sh_info.append("chrom=$2")
        l_sh_info.append("ref=$3")
        l_sh_info.append("bam_dir=%s"       % (self.dir_bismark))
        l_sh_info.append("singleC_dir=%s"   % (self.dir_singleC))
        l_sh_info.append("database=%s/$ref" % (self.Database))
        l_sh_info.append("samtools_exe=%s"  % (self.sftw_samtools))
        l_sh_info.append("bedtools_exe=%s"  % (self.sftw_bedtools))
        l_sh_info.append("bgzip_exe=%s"     % (self.sftw_bgzip))
        l_sh_info.append("tabix_exe=%s"     % (self.sftw_tabix))
        l_sh_info.append("bam2singleC_exe=%s/bam2singleC_V2_non_direct/pileup2singleC" % (self.bin))
        l_sh_info.append("""
shift
shift
shift

bam_sam_dir=$bam_dir/$meth_name
singleC_sam_dir=$singleC_dir/$meth_name
sites_ref=$database/${ref}_lambda.fa.$query/$chrom.bed

$samtools_exe view -h $bam_sam_dir/$meth_name.sort.rmdup.bam $chrom        |\\
    $samtools_exe view -Sb /dev/stdin                                       \\
    >$singleC_sam_dir/bam/$meth_name.sort.rmdup.$chrom.bam                &&\\
$samtools_exe index $singleC_sam_dir/bam/$meth_name.sort.rmdup.$chrom.bam

for query in $@
do
    awk '{OFS="\\t";print $1,$2,$2+1}'                                          \\
        $database/${ref}_lambda.fa.${query}/$chrom.bed                         |\\
    $samtools_exe mpileup  -O  -f $database/${ref}_lambda.fa -l /dev/stdin      \\
            $singleC_sam_dir/bam/$meth_name.sort.rmdup.$chrom.bam              |\\
    awk '{OFS="\\t";print $1,$2,$2,$3,$4,$5,$6,$7}'                            |\\
    $bedtools_exe intersect -sorted -loj                                        \\
        -a $database/${ref}_lambda.fa.${query}/$chrom.rev.bed -b /dev/stdin    |\\
    $bam2singleC_exe $database/${ref}_lambda.fa /dev/stdin                      \\
        $singleC_sam_dir/singleC/$chrom.${query}.bed                           &&\\
    $bgzip_exe -f $singleC_sam_dir/singleC/$chrom.${query}.bed                 &&\\
    $tabix_exe -f -p bed -s 1 -b 2 -e 3 $singleC_sam_dir/singleC/$chrom.${query}.bed.gz
done
        """)
        return l_sh_info


    def s_04_lambdaMeth(self):
        l_sh_info = []
        l_sh_info.append("meth_name=$1")
        l_sh_info.append("ref=$2")
        l_sh_info.append("singleC_dir=%s"   % (self.dir_singleC))
        l_sh_info.append("database=%s/$ref" % (self.Database))
        l_sh_info.append("samtools_exe=%s"  % (self.sftw_samtools))
        l_sh_info.append("perl_exe=%s" % (self.sftw_pl))
        l_sh_info.append("mergeSC_py=%s/merge_chrSC.py" %(self.bin))
        l_sh_info.append("""
singleC_sam_dir=$singleC_dir/$meth_name
sites_ref=$database/${ref}_lambda.fa.fai

$samtools_exe mpileup -O -f $database/${ref}_lambda.fa $singleC_sam_dir/bam/$meth_name.sort.rmdup.lambda.bam > $singleC_sam_dir/bam/$meth_name.sort.rmdup.lambda.bam.pileup
$perl_exe /date/lilin/bin/MethGC/bin/singleC_metLevel.${ref}_GCH_WCG.lambda.pl $singleC_sam_dir/bam/$meth_name.sort.rmdup.lambda.bam.pileup >$singleC_sam_dir/bam/$meth_name.sort.rmdup.lambda.single5mC
rm $singleC_sam_dir/bam/$meth_name.sort.rmdup.lambda.bam.pileup

$samtools_exe mpileup -O -f $database/${ref}_lambda.fa $singleC_sam_dir/bam/$meth_name.sort.rmdup.chrM.bam  >$singleC_sam_dir/bam/$meth_name.sort.rmdup.chrM.bam.pileup
$perl_exe /date/lilin/bin/MethGC/bin/singleC_metLevel.${ref}_HCH.chrM.pl       $singleC_sam_dir/bam/$meth_name.sort.rmdup.chrM.bam.pileup   >$singleC_sam_dir/bam/$meth_name.sort.rmdup.chrM.single5mC
rm $singleC_sam_dir/bam/$meth_name.sort.rmdup.chrM.bam.pileup
        """)
        return l_sh_info

                
    def s_05_statMeth(self):
        l_sh_info = []
        l_sh_info.append("meth_name=$1")
        l_sh_info.append("ref=$2")
        l_sh_info.append("depth=$3")
        l_sh_info.append("query=$4")
        l_sh_info.append("sc=$5")
        l_sh_info.append("singleC_dir=%s"   % (self.dir_singleC))
        l_sh_info.append("database=%s/$ref" % (self.Database))
        l_sh_info.append("stat_depth_exe=%s/depth_stat/depth_stat" % (self.bin))
        l_sh_info.append("stat_3x_depth_exe=%s/depth_stat/high_depth/depth_stat" % (self.bin))
        l_sh_info.append("sc_stat_depth_exe=%s/for_single_cell/depth_stat/depth_stat" % (self.bin))
        l_sh_info.append("sc_3x_stat_depth_exe=%s/for_single_cell/depth_stat/high_depth/depth_stat" % (self.bin))
        l_sh_info.append("""
singleC_sam_prefix=$singleC_dir/$meth_name/singleC
sites_ref=$database/${ref}_lambda.fa.fai

if [ $sc == "0" ]
#   then $stat_depth_exe -p 8 -d $depth $sites_ref $singleC_sam_prefix $query 
   then $stat_3x_depth_exe -p 8 -d 3 $sites_ref $singleC_sam_prefix $query
fi
if [ $sc == "1" ]
   then $sc_stat_depth_exe -p 8 -d $depth $sites_ref $singleC_sam_prefix $query 
   $sc_3x_stat_depth_exe -p 8 -d 3 $sites_ref $singleC_sam_prefix $query
fi

grep all $singleC_sam_prefix/$query.log|cut -f 2,3 >$singleC_sam_prefix/$query.genome.CT_total.xls
grep all $singleC_sam_prefix/$query.3x.log|cut -f 2,3 >$singleC_sam_prefix/$query.genome.3x.CT_total.xls
        """)
        return l_sh_info


    def s_06_NDR(self):
        l_sh_info = []
        l_sh_info.append("meth_name=$1")
        l_sh_info.append("chrom=$2")
        l_sh_info.append("query=GCA.GCC.GCT")
        l_sh_info.append("depth=$3")
        l_sh_info.append("sc=$4")
        l_sh_info.append("singleC_dir=%s"   % (self.dir_singleC))
        l_sh_info.append("bgzip_exe=%s"     % (self.sftw_bgzip))
        l_sh_info.append("tabix_exe=%s"     % (self.sftw_tabix))
        l_sh_info.append("python_exe=%s"    % (self.sftw_py   ))
##        l_sh_info.append("NDR_py=%s/NDR_detect/bulk_test2.py" % (self.bin))
        l_sh_info.append("NDR_py=%s/NDR_detect/bulk_test2_v2p.py" % (self.bin))
        l_sh_info.append("scNDR_py=%s/for_single_cell/NDR_detect/sc_test2.py" % (self.bin))
        l_sh_info.append("NDR_tabix_5GC=%s/NDR_detect/NDR_tabix-5GCH.pl" % (self.bin))
        l_sh_info.append("""
singleC_sam_prefix=$singleC_dir/$meth_name/singleC

if [ $sc == "0" ]
   then $python_exe $NDR_py $singleC_sam_prefix/$chrom.$query.bed.gz $singleC_sam_prefix/$query.genome.3x.CT_total.xls $depth > $singleC_sam_prefix/$chrom.$query.NDR.bed.tmp && \\
   perl $NDR_tabix_5GC $singleC_sam_prefix/$chrom.$query.NDR.bed.tmp $singleC_sam_prefix/all.GCA.GCC.GCT.3x.bedGraph.gz >$singleC_sam_prefix/$chrom.$query.3x.p10.NDR.5GC.bed && \\
   $bgzip_exe -f $singleC_sam_prefix/$chrom.$query.3x.p10.NDR.5GC.bed     && \\
   $tabix_exe -f -p bed -s 1 -b 2 -e 3 $singleC_sam_prefix/$chrom.$query.3x.p10.NDR.5GC.bed.gz 
   rm $singleC_sam_prefix/$chrom.$query.NDR.bed.tmp
fi
if [ $sc == "1" ]
   then $python_exe $scNDR_py $singleC_sam_prefix/$chrom.$query.bed.gz $singleC_sam_prefix/$query.genome.CT_total.xls $depth > $singleC_sam_prefix/$chrom.$query.NDR.bed.tmp && \\
   perl $NDR_tabix_5GC $singleC_sam_prefix/$chrom.$query.NDR.bed.tmp $singleC_sam_prefix/all.GCA.GCC.GCT.bedGraph.gz >$singleC_sam_prefix/$chrom.$query.NDR.5GC.bed && \\
   $bgzip_exe -f $singleC_sam_prefix/$chrom.$query.NDR.5GC.bed          && \\
   $tabix_exe -f -p bed -s 1 -b 2 -e 3 $singleC_sam_prefix/$chrom.$query.NDR.5GC.bed.gz
   rm $singleC_sam_prefix/$chrom.$query.NDR.bed.tmp
fi
        """)
        return l_sh_info

    def s_07_NUC(self):
        l_sh_info = []
        l_sh_info.append("meth_name=$1")
        l_sh_info.append("chrom=$2")
        l_sh_info.append("sc=$3")
        l_sh_info.append("singleC_dir=%s"   % (self.dir_singleC))
        l_sh_info.append("bgzip_exe=%s"     % (self.sftw_bgzip))
        l_sh_info.append("tabix_exe=%s"     % (self.sftw_tabix))
        l_sh_info.append("python_exe=%s"    % (self.sftw_py   ))
##        l_sh_info.append("NUC_py=%s/NDR_detect/bulk_nucleosome_detect_v2.py" % (self.bin))
        l_sh_info.append("NUC_py=%s/NDR_detect/bulk_nucleosome_detect_v2_v2p.py" % (self.bin))
        l_sh_info.append("scNUC_py=%s/for_single_cell/NDR_detect/sc_nucleosome_detect_v2.py" % (self.bin))
        l_sh_info.append("NUC_tabix_3GC=%s/NDR_detect/NUC_tabix_3GC.pl" % (self.bin))
        l_sh_info.append("""
singleC_sam_prefix=$singleC_dir/$meth_name/singleC

if [ $sc == "0" ]
   then  $python_exe $NUC_py $singleC_sam_prefix/$chrom.GCA.GCC.GCT.bed.gz $singleC_sam_prefix/GCA.GCC.GCT.genome.3x.CT_total.xls 4 >$singleC_sam_prefix/$chrom.GCA.GCC.GCT.nucleosome.tmp     && \\
   perl $NUC_tabix_3GC $singleC_sam_prefix/$chrom.GCA.GCC.GCT.nucleosome.tmp $singleC_sam_prefix/all.GCA.GCC.GCT.3x.bedGraph.gz  > $singleC_sam_prefix/$chrom.GCA.GCC.GCT.3x.nucleosome.3GC.bed      && \\
   $bgzip_exe -f $singleC_sam_prefix/$chrom.GCA.GCC.GCT.3x.nucleosome.3GC.bed    &&\\
   $tabix_exe -f -p bed -s 1 -b 2 -e 3 $singleC_sam_prefix/$chrom.GCA.GCC.GCT.3x.nucleosome.3GC.bed.gz    && \\
   rm $singleC_sam_prefix/$chrom.GCA.GCC.GCT.nucleosome.tmp
fi
if [ $sc == "1" ]
   then  $python_exe $scNUC_py $singleC_sam_prefix/$chrom.GCA.GCC.GCT.bed.gz $singleC_sam_prefix/GCA.GCC.GCT.genome.3x.CT_total.xls 3 >$singleC_sam_prefix/$chrom.GCA.GCC.GCT.nucleosome.tmp     && \\
   perl $NUC_tabix_3GC $singleC_sam_prefix/$chrom.GCA.GCC.GCT.nucleosome.tmp $singleC_sam_prefix/all.GCA.GCC.GCT.bedGraph.gz  > $singleC_sam_prefix/$chrom.GCA.GCC.GCT.nucleosome.3GC.bed      && \\
   $bgzip_exe -f $singleC_sam_prefix/$chrom.GCA.GCC.GCT.nucleosome.3GC.bed    &&\\
   $tabix_exe -f -p bed -s 1 -b 2 -e 3 $singleC_sam_prefix/$chrom.GCA.GCC.GCT.nucleosome.3GC.bed.gz    && \\
   rm $singleC_sam_prefix/$chrom.GCA.GCC.GCT.nucleosome.tmp
fi
        """)
        return l_sh_info


    def s_08_NDR_IGV(self):
        l_sh_info = []
        l_sh_info.append("meth_name=$1")
        l_sh_info.append("query=$2")
        l_sh_info.append("ref=$3")
        l_sh_info.append("singleC_dir=%s"   % (self.dir_singleC))
        l_sh_info.append("database=%s/$ref" % (self.Database))
        l_sh_info.append("NDR_dir=%s"          % (self.dir_NDR))
#        l_sh_info.append("TDFs_dir=%s"         % (self.dir_TDF))
        l_sh_info.append("mergeSC_py=%s/merge_chrSC.py" %(self.bin)) 
#        l_sh_info.append("igvtools_exe=%s" % (self.sftw_igvtools))
        l_sh_info.append("python_exe=%s"    % (self.sftw_py   ))
        l_sh_info.append("ucsc_dir=%s"     % (self.sftw_ucsc_dir))
        l_sh_info.append("bgzip_exe=%s"     % (self.sftw_bgzip))
        l_sh_info.append("tabix_exe=%s"     % (self.sftw_tabix))
        l_sh_info.append("""
singleC_sam_prefix=$singleC_dir/$meth_name/singleC
sites_ref=$database/${ref}_lambda.fa.fai

$python_exe $mergeSC_py -r $sites_ref -c $query $singleC_sam_prefix

$bgzip_exe -c $singleC_sam_prefix/all.$query.bed  > $singleC_sam_prefix/all.$query.bed.gz  &&\\
$tabix_exe -f -p bed -s 1 -b 2 -e 3  $singleC_sam_prefix/all.$query.bed.gz


awk '{OFS="\\t";if ($4+$5>=1) print $1,$2-1,$3,$5/($4+$5)}' $singleC_sam_prefix/all.${query}.bed >$singleC_sam_prefix/all.${query}.bedGraph &&\\
$bgzip_exe -c $singleC_sam_prefix/all.${query}.bedGraph > $singleC_sam_prefix/all.${query}.bedGraph.gz &&\\
$tabix_exe -f -p bed -s 1 -b 2 -e 3  $singleC_sam_prefix/all.${query}.bedGraph.gz

$ucsc_dir/bedGraphToBigWig                                             \\
    $singleC_sam_prefix/all.${query}.bedGraph				\\
    $sites_ref  $singleC_sam_prefix/all.${query}.bw

awk '{OFS="\\t";if ($4+$5>=3) print $1,$2-1,$3,$5/($4+$5)}' $singleC_sam_prefix/all.${query}.bed >$singleC_sam_prefix/all.${query}.3x.bedGraph &&\\
$bgzip_exe -c $singleC_sam_prefix/all.${query}.3x.bedGraph >$singleC_sam_prefix/all.${query}.3x.bedGraph.gz &&\\
$tabix_exe -f -p bed -s 1 -b 2 -e 3  $singleC_sam_prefix/all.${query}.3x.bedGraph.gz

$ucsc_dir/bedGraphToBigWig                                             \\
    $singleC_sam_prefix/all.${query}.3x.bedGraph			\\
    $sites_ref  $singleC_sam_prefix/all.${query}.3x.bw

rm $singleC_sam_prefix/all.$query.bed
rm $singleC_sam_prefix/all.${query}.bedGraph
rm $singleC_sam_prefix/all.${query}.3x.bedGraph
       """)
        return l_sh_info

    def s_09_mergeNUC_NDR_rm(self):
        l_sh_info = []
        l_sh_info.append("meth_name=$1")
        l_sh_info.append("query=GCA.GCC.GCT")
        l_sh_info.append("ref=$2")
        l_sh_info.append("sc=$3")
        l_sh_info.append("database=%s/$ref" % (self.Database))
        l_sh_info.append("singleC_dir=%s"   % (self.dir_singleC))
        l_sh_info.append("NDR_dir=%s"          % (self.dir_NDR))
        l_sh_info.append("mergeSC_py=%s/merge_chrSC.py" %(self.bin))
        l_sh_info.append("python_exe=%s" % (self.sftw_py))
        l_sh_info.append("""
singleC_sam_prefix=$singleC_dir/$meth_name/singleC
sites_ref=$database/${ref}_lambda.fa.fai

if [ $sc == "0" ]
   then $python_exe $mergeSC_py -r $sites_ref -c $query.nucleosome.3GC $singleC_sam_prefix &&\\
$python_exe $mergeSC_py -r $sites_ref -c $query.NDR.5GC $singleC_sam_prefix
#ln -s $singleC_sam_prefix/all.$query.NDR_v2.5GC.bed $NDR_dir/$meth_name.$query.NDR_v2.5GC.bed 
#ln -s $singleC_sam_prefix/all.$query.nucleosome_v2.3GC.bed  $NDR_dir/$meth_name.$query.nucleosome_v2.3GC.bed
fi

if [ $sc == "1" ]
   then $python_exe $mergeSC_py -r $sites_ref -c $query.nucleosome.3GC $singleC_sam_prefix &&\\
$python_exe $mergeSC_py -r $sites_ref -c $query.NDR.5GC $singleC_sam_prefix   

#for i in {2..22}
#do
#   rm $singleC_dir/$meth_name/bam/$meth_name.sort.rmdup.chr$i.bam*
#done

#for i in {1..22} X Y M
#do
#   rm $singleC_sam_prefix/chr$i.ACG.TCG.bed.gz*
#   rm $singleC_sam_prefix/chr$i.GCA.GCC.GCT.bed.gz*
#done
fi
       """)
        return l_sh_info
 
    def s_11_plotDist(self):
        l_sh_info = []
        l_sh_info.append("mergeRatio=$1")
        l_sh_info.append("db_file=$2")
        l_sh_info.append("outprefix=$3")
        l_sh_info.append("sam_file=$4")
        l_sh_info.append("python_path=%s" % (self.sftw_py))
#        l_sh_info.append("pltDist_py=%s/module_plotDist_geneTSS-body.py" % (self.bin))
        l_sh_info.append("pltDist_py=%s/module_plotDist_ALL.py" % (self.bin))
        l_sh_info.append("""
$python_path $pltDist_py $mergeRatio $db_file $outprefix $sam_file
        """)
        return l_sh_info


    def s_12_NDR_flanking(self):
        l_sh_info = []
        l_sh_info.append("meth_name=$1")
        l_sh_info.append("ref=$2")
        l_sh_info.append("database=%s" % (self.Database))
#        l_sh_info.append("database=%s" % (self.Database_lyn))
        l_sh_info.append("singleC_dir=%s"   % (self.dir_singleC))
        l_sh_info.append("NDR_dir=%s"          % (self.dir_NDR))
        l_sh_info.append("NDR_flank_dir=%s"    % (self.dir_NDR_flank)) 
        l_sh_info.append("peak_plot_exe=%s/Peak_plot/peak_plot"        % (self.bin))
        l_sh_info.append("cut_NDR_support_dep_exe=%s/cut_NDR_support_dep.py"        % (self.bin))
        l_sh_info.append("bedtools_exe=%s" % (self.sftw_bedtools))
        l_sh_info.append("python_exe=%s"    % (self.sftw_py   ))
        l_sh_info.append("""
singleC_sam_prefix=$singleC_dir/$meth_name/singleC

ln -s $singleC_sam_prefix/all.GCA.GCC.GCT.3x.p10.NDR.5GC.bed $NDR_dir/$meth_name.GCH.3x.p10.NDR.5GC.bed
$bedtools_exe intersect -v -sorted                                              \\
    -a $NDR_dir/$meth_name.GCH.3x.p10.NDR.5GC.bed                              \\
    -b $database/$ref/refGene.TSS_2000.bed                                      |\\
    $python_exe $cut_NDR_support_dep_exe /dev/stdin                             \\
    | uniq >$NDR_dir/distal_region/$meth_name.GCH.3x.p10.NDR.5GC.distal.cnt10.bed

$bedtools_exe intersect  -u  -sorted                                          \\
    -a $NDR_dir/$meth_name.GCH.3x.p10.NDR.5GC.bed                              \\
    -b $database/$ref/refGene.TSS_2000.bed                                      |\\
    $python_exe $cut_NDR_support_dep_exe /dev/stdin                 \\
    | uniq >$NDR_dir/TSS_region/$meth_name.GCH.3x.p10.NDR.5GC.TSS.cnt10.bed

#awk -v OFS="\\t" '{print $1,$2,$3,"+" }' $NDR_dir/$meth_name.GCA.GCC.GCT.NDR.bed |\\
#    $peak_plot_exe $singleC_sam_prefix  /dev/stdin  $NDR_flank_dir/$meth_name.all

#awk -v OFS="\\t" '{print $1,$2,$3,"+" }' $NDR_dir/distal_region/$meth_name.GCA.GCC.GCT.NDR.distal.cnt10.bed |\\
#    $peak_plot_exe $singleC_sam_prefix  /dev/stdin  $NDR_flank_dir/$meth_name.distal.cnt10

#awk -v OFS="\\t" '{print $1,$2,$3,"+" }' $NDR_dir/TSS_region/$meth_name.GCA.GCC.GCT.NDR.TSS.cnt10.bed |\\
#    $peak_plot_exe $singleC_sam_prefix  /dev/stdin  $NDR_flank_dir/$meth_name.TSS.cnt10
       """)
        return l_sh_info

    def s_13_NDR_motif(self):
        l_sh_info = []
        l_sh_info.append("meth_name=$1")
        l_sh_info.append("ref=$2")
        l_sh_info.append("database=%s/$ref" % (self.Database))
        l_sh_info.append("perl_exe=%s" % (self.sftw_pl))
        l_sh_info.append("sftw_hommer=%s" % (self.sftw_hommer))
        l_sh_info.append("singleC_dir=%s"   % (self.dir_singleC))
        l_sh_info.append("NDR_dir=%s"          % (self.dir_NDR))
#        l_sh_info.append("NDR_flank_dir=%s"    % (self.dir_NDR_flank))
        l_sh_info.append("bin=%s"          % (self.bin))
        l_sh_info.append("bedtools_exe=%s" % (self.sftw_bedtools))
        l_sh_info.append("""
singleC_sam_prefix=$singleC_dir/$meth_name/singleC

#$perl_path $sftw_hommer                                                     \\
#    $NDR_dir/$meth_name.GCA.GCC.GCT.NDR.bed $ref                            \\
#    $NDR_dir/$meth_name.GCA.GCC.GCT.NDR.given -size 2000 -len 8 -S 100      \\
#    -mknown /data/Analysis/huboqiang/Database_Meth/mm9/motifs.motif

#$perl_path $sftw_hommer
#    $NDR_dir/distal_region/$meth_name.GCA.GCC.GCT.NDR.distal.cnt10.bed  $ref     \\
#    $NDR_dir/$meth_name.GCA.GCC.GCT.NDR.distal.cnt10.given -size 2000 -len 8 -S 100 \\
#    -mknown /data/Analysis/huboqiang/Database_Meth/mm9/motifs.motif

#$perl_path $sftw_hommer                                                     \\
#    $NDR_dir/TSS_region/$meth_name.GCA.GCC.GCT.NDR.TSS.cnt10.bed    $ref    \\
#    $NDR_dir/$meth_name.GCA.GCC.GCT.NDR.TSS.cnt10.given -size 2000 -len 8 -S 100  \\
#    -mknown /data/Analysis/huboqiang/Database_Meth/mm9/motifs.motif

#$perl_path $sftw_hommer                                                     \\
#    $NDR_dir/$meth_name.GCA.GCC.GCT.NDR.bed $ref                            \\
#    $NDR_dir/$meth_name.GCA.GCC.GCT.NDR -size 2000 -len 8 -S 100

$perl_exe $sftw_hommer                                                     \\
    $NDR_dir/distal_region/$meth_name.GCH.3x.p10.NDR.5GC.distal.cnt10.bed  $ref\\
    $NDR_dir/$meth_name.GCH.NDR.5GC.distal -size 2000 -len 8 -S 100

$perl_exe $sftw_hommer                                                     \\
    $NDR_dir/TSS_region/$meth_name.GCH.3x.p10.NDR.5GC.TSS.cnt10.bed  $ref      \\
    $NDR_dir/$meth_name.GCH.NDR.5GC.TSS -size 2000 -len 8 -S 100
       """)
        return l_sh_info

