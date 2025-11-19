#include "bed_read.h"
#include "gzstream.h"

using namespace std;

void load_chrlen(string fa_fai,map< string,unsigned > &ChrLen){
   ifstream infile;
   string file_name = fa_fai;
   infile.open(file_name.c_str());
   if ( ! infile ){
      cerr << "fail to open input file" << file_name << endl;
      exit(0);
   }
   string lineStr;
   while (getline(infile,lineStr,'\n')){
      if (lineStr[0] == ' ' || lineStr[0] == '\n'){
         continue;
      }
      vector<string> lineVec;
      boost::split(lineVec,lineStr,boost::is_any_of(":, \t\n"), boost::token_compress_on);
      unsigned length =  boost::lexical_cast<unsigned>(lineVec[1]);
      ChrLen[lineVec[0]] = length;
   }
   infile.close();
}

void load_bins(map< string,unsigned > &ChrLen, map< string,vector<float> > &ChrBin_ratio,map< string,vector<int> > &ChrBin_Met,map< string,vector<int> > &ChrBin_UnMet,map< string,vector<int> > &ChrBin_site,int bin){
   map< string,unsigned >::iterator I_it;
   for ( I_it=ChrLen.begin(); I_it!=ChrLen.end();I_it++ ){
      string chr = I_it->first;
      unsigned    len = I_it->second;
      unsigned    bin_len = len/bin + 1;
      vector<int>   len_zeros_site( bin_len,0);
      vector<int>   len_zeros_Met( bin_len,0);
      vector<int>   len_zeros_UnMet( bin_len,0);
      vector<float> len_zeros_ratio(bin_len,0.0);
       ChrBin_site.insert (  map< string,vector<int>   >::value_type(chr,len_zeros_site  ) ); 
       ChrBin_Met.insert  (  map< string,vector<int>   >::value_type(chr,len_zeros_Met   ) );
       ChrBin_UnMet.insert(  map< string,vector<int>   >::value_type(chr,len_zeros_UnMet ) );
       ChrBin_ratio.insert(  map< string,vector<float> >::value_type(chr,len_zeros_ratio ) );
  }
}

void bed_bin(string in_bed,map< string,vector<float> > &ChrBin_ratio,map< string,vector<int> > &ChrBin_Met,map< string,vector<int> > &ChrBin_UnMet,map< string,vector<int> > &ChrBin_site,int bin){
   ifstream infile;
   infile.open(in_bed.c_str());
   if ( ! infile ){
      cerr << "fail to open input file" << in_bed << endl;
      exit(0);
   }
   string lineStr;
   while (getline(infile,lineStr,'\n')){
      if (lineStr[0] == ' ' || lineStr[0] == '\n'){
         continue;
      }
      vector<string> lineVec;
      boost::split(lineVec,lineStr,boost::is_any_of(":, \t\n"), boost::token_compress_on);
		// my @f = split(/\t/,$lineStr);
		
      string chr = lineVec[0];
      unsigned begin = boost::lexical_cast<unsigned>(lineVec[1]);
      int      Met   = boost::lexical_cast<int>(   lineVec[3] );
      int      UnMet = boost::lexical_cast<int>(   lineVec[2] );
      float    ratio = boost::lexical_cast<float>( lineVec[4] );
      unsigned bin_n = begin / bin;
      if ( begin%bin > bin/2 ){
         bin_n += 1;
      }
//      cout << bin_n << "\t" <<  ChrBin_ratio[chr].size() << endl;
      if (bin_n >= ChrBin_ratio[chr].size()){
//         cout << chr << "\t" << bin_n << endl;
          bin_n =  ChrBin_ratio[chr].size() - 1;
      }
      ChrBin_ratio[chr][bin_n] += ratio;
      ChrBin_Met  [chr][bin_n] += Met  ;
      ChrBin_UnMet[chr][bin_n] += UnMet;
      ChrBin_site[chr][bin_n]  += 1;
   }
   infile.close();
}

void report_bin(map< string,vector<float> > &ChrBin_ratio,map< string,vector<int> > &ChrBin_Met,map< string,vector<int> > &ChrBin_UnMet,map< string,vector<int> > &ChrBin_site, string out_bed,int bin){
	//ogzstream outFile;
	ofstream outFile;
   outFile.open(out_bed.c_str());
   if ( ! outFile )
   {
      cerr << "fail to open input file" << outFile << endl;
      exit(0);
   }
   map< string,vector<int>  >::iterator I_it;
   for ( I_it=ChrBin_site.begin(); I_it!=ChrBin_site.end();I_it++ ){
      string chr = I_it->first;
      for (unsigned i=0;i<ChrBin_site[chr].size();i++){
         unsigned begin = i*bin + 1;
         unsigned end   = (i+1)*bin; 
			float sum_ratio = ChrBin_ratio[chr][i];
			int   sum_Met   = ChrBin_Met[chr][i];
			int   sum_UnMet = ChrBin_UnMet[chr][i];
			int   sum_site  = ChrBin_site[chr][i];
         outFile << chr << "\t" << begin << "\t" << end << "\t" << sum_ratio << "\t"  << sum_site << "\t" << sum_Met<< "\t"  << sum_UnMet << endl;
      }
   }
   outFile.close();
}
