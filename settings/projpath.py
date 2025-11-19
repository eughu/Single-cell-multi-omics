#-*- coding:utf-8 -*-
from __future__ import division
import re
import sys
import os
import subprocess
import numpy   as np
import cPickle as pickle
import pandas  as pd

class LoadSamp(object):
    def __init__(self):
        super( LoadSamp,self ).__init__()
        
    def load_RNA_samInfo(self,RNA_infile):
        self.samInfo_RNA_infile  = RNA_infile
        self.samInfo_pd_RNA  = pd.read_csv(self.samInfo_RNA_infile , sep="\t")
        
      
class DirSystem(object):
    def __init__(self):
        super( DirSystem,self ).__init__()
        self.home_dir               = os.path.abspath('./')
        path_file                   = os.path.realpath(__file__)
        self.dir_raw_data           = "%s/00.0.raw_data"               % (self.home_dir)
        self.dir_trim_data          = "%s/00.1.trim_data"              % (self.home_dir)
        self.dir_bismark            = "%s/01.bam"                      % (self.home_dir)
        self.dir_singleC            = "%s/02.SingleC"                  % (self.home_dir)
        self.dir_NDR                = "%s/03.NDR"                      % (self.home_dir)
 
        self.dir_StatInfo           = "%s/StatInfo"                    % (self.home_dir)
        self.path                   = "%s"     % ("/".join(path_file.split('/')[:-2]))
        self.bin                    = "%s/bin" % ("/".join(path_file.split('/')[:-2]))

        if not os.path.exists( self.dir_StatInfo ):
            os.mkdir( self.dir_StatInfo )
        #####################################
        ######## Revise this path!!! ########
        #####################################
        self.Database               = "/date/lilin/Database_Meth"
#        self.Database_lyn           = "/Share/BP/lilin/Database/MethGC"

    def define_scripts(self, s_idx):
        dir_script = "%s/scripts"      % (self.home_dir)
        self.scripts = "%s/scripts/%s" % (self.home_dir, s_idx)

        if not os.path.exists(dir_script):
            os.mkdir(dir_script)

        if not os.path.exists(self.scripts):
            os.mkdir(self.scripts)

        

class UsedSoftware(object):
    def __init__(self):
        super( UsedSoftware,self ).__init__()
        
        ##############################################
        ######## Revise the following path!!! ########
        ##############################################
        self.sftw_py        = "/date/zhangshu/software/anaconda2/bin/python"
	#"/date/lilin/software/anaconda2/bin/python"
        self.sftw_pl        = "/usr/bin/perl"
        
        self.sftw_samtools  = "/date/lilin/software/samtools-0.1.18/samtools"
        self.sftw_bedtools  = "/date/lilin/software/bedtools/bin/bedtools"
        self.sftw_bgzip     = "/date/lilin/software/tabix-0.2.6/bgzip"
        self.sftw_tabix     = "/date/lilin/software/tabix-0.2.6/tabix"
        self.sftw_ucsc_dir  = "/datf/xiehaoling/software/UCSC"
        
        self.sftw_trim      = "/date/lilin/software/trim_galore"
        self.sftw_bismark   = "/date/lilin/software/bismark_v0.7.6"
        self.sftw_changeID  = "/date/lilin/software/ChangeReadID.pl"
        self.sftw_bowtie_dir= "/date/lilin/software/bowtie-1.0.0"
#        self.sftw_igvtools  = "/WPS1/lilin/Software/IGVTools/igvtools" 
        self.sftw_hommer    = "/date/lilin/software/hommer/bin/findMotifsGenome.pl"

class ProjInfo(LoadSamp,DirSystem,UsedSoftware):
    def __init__(self):
        super( ProjInfo,self ).__init__()

