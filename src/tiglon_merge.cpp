#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<algorithm>
#include<vector>
#include<map>
#include "tiglon_merge.h"
using namespace std;
typedef vector<int> path_t;
ofstream outgtf;
int SampleSize;
double Coverage = 0;
double SEED = 0;
struct info
{
    int number;
    double seed;
    int normal_edge_number;
    int partial_edge_number;
    int mj_number;
    int Nmj_number;
    int seed_sample_number;
};
map<string, vector<double> > id_cov_map;
map<string, info > id_info_map;

void load_info(char*file)
{
    ifstream in(file);
    string s;
    istringstream istr;
    while(getline(in,s))
    {
        istr.str(s);
	string temp,id;
	double seed;
	int normal_edge,partial_edge,mj,Nmj,seed_sample_number;
	istr>>temp>>id>>seed>>normal_edge>>partial_edge>>temp>>temp>>temp>>mj>>Nmj>>seed_sample_number;
	istr.clear();
	int N = atoi(id.substr(id.length() - 3,1).c_str());
	info info_={N,seed,normal_edge,partial_edge,mj,Nmj,seed_sample_number};
	id_info_map[id] = info_;
    }
    //cout<<"id_info_map.size(): "<<id_info_map.size()<<endl;
}
void process(string tranid,vector<string> oneTrans)
{
    if(id_cov_map.find(tranid) != id_cov_map.end() )//&& id_cov_map.find(tran_id) != id_cov_map.end())
    {
	    double cov = 0;
	    for(size_t i=0;i<id_cov_map[tranid].size();i++)
		    cov += id_cov_map[tranid][i];
	    cov = cov/(1.0*SampleSize);
	    bool flag = false;
	    if(cov > SEED) flag = true;
		    

	    if(flag)  
	    {
		    for(size_t i=0;i<oneTrans.size();i++)
		    outgtf<<oneTrans[i]<<" cov \""<<cov<<"\";"<<endl;
	    }
     }
    return;
}
void get_final_results(char*file)
{

    ifstream in(file);
    istringstream istr;
    string s;

    string chr, strand, lable;
    string tranid;
    string temp, Cov_s;
    double Cov;
    int exon_l,exon_r;
    vector<int> vecExon;
    vector<string> oneTrans;
    getline(in,s);
    istr.str(s);
    while(istr>>temp)
	    if( temp == "transcript_id")
		    istr>>tranid;
    istr.clear();
    oneTrans.push_back(s);

    while(getline(in,s))
    {
	istr.str(s);
 	string current_id;
	string curr_chr, curr_strand;
	int curr_el, curr_er;

	istr>>curr_chr>>temp>>lable>>curr_el>>curr_er>>temp>>curr_strand;


	while(istr>>temp)
	{
	    if( temp == "transcript_id")
		istr>>current_id;
	}
	if(current_id ==  tranid)
	{
	  oneTrans.push_back(s);
	}
	else
	{
	  process(tranid,oneTrans);
	  tranid = current_id;
	  oneTrans.clear();
	  oneTrans.push_back(s);
	}

	istr.clear();
    }
    process(tranid,oneTrans);
    return;

}
//./exe a.gtf b.gtf .. N.gtf input.info input.gtf output.gtf SEED sampleSize //only remove repeat ones
//./exe a.gtf b.gtf .. N.gtf input.gtf output.gtf SEED sampleSize//BamMerging
int main(int argc,char* argv[])
{
    //load_transref(argv[1],intron_trans_map,Chr_Junc_map);
    //cout<<argc<<endl;
    string sample_S = argv[argc-1];
    SampleSize=atoi(sample_S.c_str());
    string S = argv[argc-2];
    SEED = atof(S.c_str());
    //cout<<"filter: "<<SEED<<endl;

    outgtf.open(argv[argc-3]);

    for(int i=1;i<=argc-5;i++){
	//cout<<"load: "<<argv[i]<<" "<<i<<" sample..."<<endl;
        load_transref(argv[i],id_cov_map);
    
    }

    get_final_results(argv[argc-4]);
    outgtf.close();
    return 0;
}
