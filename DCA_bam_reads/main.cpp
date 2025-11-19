#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>
#include <numeric>
#include <inttypes.h>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>


#include "bgzf.h"
#include "tabix.h"
#include "knetfile.h"

//#include "split.h"
//#include "Fasta.h"

using namespace std;

const int g_mindep = 3;

int cnt_ovlp_tabix(tabix_t *t, string &region, vector<int> &V_index,  vector<string> &V_strand){
    int i, len;
    ti_iter_t iter;
    const char *s;
    const ti_conf_t *idxconf;

    if (ti_lazy_index_load(t) < 0 ) {
        fprintf(stderr,"[tabix] failed to load the index file.\n");
        return 1;
    }
    idxconf = ti_get_conf(t->idx);

    int32_t tid, beg, end;
    
    vector<double> V_value;
    
    int ovlp_cnt = 0;
    if (ti_parse_region(t->idx, region.c_str(), &tid, &beg, &end) == 0) {
        iter = ti_queryi(t, tid, beg, end);
        while ((s = ti_read(t, iter, &len)) != 0) {
            string lineCons = s;
            vector<string> V_Cons;
            boost::split(V_Cons, lineCons, boost::is_any_of(" \t\n"), boost::token_compress_on);
            int32_t pos = boost::lexical_cast<int32_t>(V_Cons[1]);
//            cout << pos << "\t" << beg << "\t" << pos-beg << endl;
            V_index.push_back(pos - beg );
            V_strand.push_back(V_Cons[3]);
            ++ovlp_cnt;
        }
        ti_iter_destroy(iter);
    }
    return ovlp_cnt;
}


double ratio_pass(string seq, int seq_isRev, vector<int> &V_index, vector<string> &V_strand, int &is_meth_cnt){
    vector<int> V_is_Meth;
    for (int i = 0;i < V_index.size(); ++i){
        int pos = V_index[i];
        char pos_seq = seq[pos];
        string pos_str = V_strand[i];
        int is_meth = -1;
        
        if (pos_str == "+" and not seq_isRev){
            if (pos_seq == 'C'){
                is_meth = 1;
            }
            if (pos_seq == 'T'){
                is_meth = 0;
            }
        }
        else{
            if (pos_str == "-" and seq_isRev){
                if (pos_seq == 'g'){
                    is_meth = 1;
                }
                if (pos_seq == 'a'){
                    is_meth = 0;
                }
            }
        }
//        cout << i << "\t"<< pos_seq << "\t"  << pos_str << "\t" << seq_isRev << "\t"<< pos << "\t" << is_meth << endl;
        if(is_meth >= 0){
               V_is_Meth.push_back(is_meth);
        }
    }
    double mean = -1;
    if (V_is_Meth.size() >= g_mindep){
        double sum = std::accumulate(V_is_Meth.begin(), V_is_Meth.end(), 0.0);
        mean = sum / V_is_Meth.size();
    }
    is_meth_cnt = V_is_Meth.size();
    return mean;
}
/*
    if (strand == "+" and type){
        if (base == 'C'){
            cnt_met += 1;
        }
        if (base == 'T'){
            cnt_umt += 1;
        }
    }
    if (strand == "-" and not type){
        if (base == 'G'){
            cnt_met += 1;
        }
        if (base == 'A'){
            cnt_umt += 1;
        }
    }
*/







int read_bam_interval( string sam_file, string bed_file_CG, string bed_file_GC){
    tabix_t *t_CG;
    tabix_t *t_GC;

    if ((t_CG = ti_open(bed_file_CG.c_str(), 0)) == 0) {
		fprintf(stderr, "[main] fail to open the data file.\n");
		return 1;
	}
    if ((t_GC = ti_open(bed_file_GC.c_str(), 0)) == 0) {
		fprintf(stderr, "[main] fail to open the data file.\n");
		return 1;
	}
    
    ifstream infile;
    infile.open(sam_file.c_str());
    if ( ! infile ){
        cerr << "Cannot find interval file" << sam_file << endl;
        exit(0);
    }
    
    string lineStr;
    while(getline( infile, lineStr, '\n' ))
 	{
        if(lineStr[0] == ' ' || lineStr[0] == '\n' || lineStr[0] == '@'){
            continue;
        }
        
        vector<string> lineVec;
        boost::split(lineVec,lineStr,boost::is_any_of("\t"), boost::token_compress_on);
	    
        string seq = lineVec[9];
        int flag = boost::lexical_cast<int>(lineVec[1]);
        int32_t begin_read = boost::lexical_cast<int32_t>(lineVec[3]);
        int32_t end_read   = begin_read + seq.size();
        string region = lineVec[2] + ":" + lineVec[3] + "-" + boost::lexical_cast<string>(end_read);
        
//        cout << region << "\t" << seq << endl;
        
        vector<int>    V_index_CG;
        vector<int>    V_index_GC;
        vector<string> V_strand_CG;
        vector<string> V_strand_GC;

        int seq_isRev = flag & 1<< 4;

      
        int ovlp_CG = cnt_ovlp_tabix(t_CG, region, V_index_CG,  V_strand_CG);
        int ovlp_GC = cnt_ovlp_tabix(t_GC, region, V_index_GC,  V_strand_GC);

        if(ovlp_CG < g_mindep or ovlp_GC < g_mindep){
            continue;
        }
        
        int is_meth_cnt_CG = 0;
        int is_meth_cnt_GC = 0;
        
        double ratio1 = ratio_pass(seq, seq_isRev, V_index_CG, V_strand_CG, is_meth_cnt_CG);
        double ratio2 = ratio_pass(seq, seq_isRev, V_index_GC, V_strand_GC, is_meth_cnt_GC);
        
        if (ratio1 == -1 or ratio2 == -1){
            continue;
        }

        if ( (ratio1 >= 0.9 and ratio2 <= 0.1) or (ratio1 <= 0.1 and ratio2 >= 0.9)   ){
            cout << region << "\t" << seq << "\t" << is_meth_cnt_CG << "\t" << is_meth_cnt_GC << "\t" << ratio1 << "\t" << ratio2 << endl;
        }

    }
    infile.close();
    ti_close(t_CG);
    ti_close(t_GC);
    
}


void usage()
{  
   cerr << endl;
   cerr << "Usage  :   ./DCA_bam /datd/huboqiang/test_NOM/02.SingleC/mESC_gF28_1/bam/mESC_gF28_1.sort.rmdup.chr10.sam /data/Analysis/huboqiang/Database_Meth/mm9/mm9_lambda.fa.ACG.TCG.bed.gz  /data/Analysis/huboqiang/Database_Meth/mm9/mm9_lambda.fa.GCA.GCC.GCT.bed.gz  " << endl;
   cerr << endl;
   cerr << "Example:  samtools view /datd/huboqiang/test_NOM/02.SingleC/mESC_gF28_1/bam/mESC_gF28_1.sort.rmdup.chr10.bam | ./DCA_bam /dev/stdin /data/Analysis/huboqiang/Database_Meth/mm9/mm9_lambda.fa.ACG.TCG.bed.gz  /data/Analysis/huboqiang/Database_Meth/mm9/mm9_lambda.fa.GCA.GCC.GCT.bed.gz " << endl;
   cerr << "Options:  " << endl;
   cerr << "               -s  Enzyme sites. Now only supports for CviAII(CATG) and DpnII(GATC). [CATG]" << endl;
   cerr << "               -h  get help information"   << endl;
   exit (0);
}

int main(int argc, char *argv[])
{
    int c;
    while ( (c=getopt(argc,argv,"s:h")) != -1 ){
        switch(c)
        {
            case 'h' : usage();break;
            default : usage();
        }
    }
    
    if (argc < 3) usage();
    string sam_file = argv[optind++];
    string bed_file_CG = argv[optind++];
    string bed_file_GC = argv[optind++];

    read_bam_interval( sam_file, bed_file_CG, bed_file_GC );
}
