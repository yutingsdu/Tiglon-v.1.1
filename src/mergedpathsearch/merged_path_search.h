#ifndef PATH_SEARCH_H
#define PATH_SEARCH_H


#include <stdlib.h>
#include<vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "merged_simplify_graph.h"
#include "merged_junction_graph.h"
using namespace std;


typedef vector<int> path_t;
extern ofstream out_graph;
extern double SEED_Filter;
extern string out_name;
extern int trans_id ;
extern int Path_length ;
extern int rg_index;

int graph_number;
class PathSearch
{
  public:
    vector<int> edge_out,edge_in,exon_l,exon_r;
    vector<double> edge_weight,exon_cov;

    vector<vector<int> > Vec_edges;
    vector<double> Vec_weights;

    vector<path_t > Unused_pair_paths;
    vector<double> Unused_pair_paths_cov;

    //int line, strand, XS_plus, XS_minus;
    string strand,chr;
    vector<string> Node_seq_1;
    vector<int>Node_seq_2,Node_seq_3;
    vector<double>Node_seq_4;
    vector< vector< pair<int,double> > >junctions_MappingInfo;
    vector< pair<int,int> > junctions;

  public:
  PathSearch(vector<int>& edge_out_, vector<int>& edge_in_, vector<double>& edge_weight_,
	     vector<path_t >& Unused_pair_paths_, vector<double>& Unused_pair_paths_cov_,
	     vector<string>& Node_seq_1_, vector<int>& Node_seq_2_, vector<int>& Node_seq_3_, vector<double>& Node_seq_4_,
	     string strand_,string chr_)
	    //int line_,int strand_,int XS_plus_,int XS_minus_)
  {

	edge_out = edge_out_;
	edge_in = edge_in_;
	edge_weight = edge_weight_;

	Unused_pair_paths = Unused_pair_paths_;
	Unused_pair_paths_cov = Unused_pair_paths_cov_;

	Node_seq_1 = Node_seq_1_;
	Node_seq_2 = Node_seq_2_; Node_seq_3 = Node_seq_3_;
	Node_seq_4 = Node_seq_4_;

	//line = line_; strand = strand_; XS_plus = XS_plus_; XS_minus = XS_minus_;
	strand = strand_;
	chr = chr_;
	for(size_t i=0;i<Node_seq_1_.size();i++)
	{
	    exon_l.push_back(Node_seq_2_[i]);
	    exon_r.push_back(Node_seq_3_[i]);
	    exon_cov.push_back(Node_seq_4_[i]);
	}

	Unused_pair_paths_.clear();
    	Unused_pair_paths_cov_.clear();
        Node_seq_1_.clear();
        Node_seq_2_.clear(); Node_seq_3_.clear();
        Node_seq_4_.clear();
        edge_out_.clear(); edge_in_.clear(); edge_weight_.clear();	
  }

  void path_search()
  {
    //cout<<"Graph: "<<rg_index<<endl;
    Vec_edges.clear(); Vec_weights.clear();

    vector<int> Single_nodes, invalid_nodes;

    Combing simplify;
    simplify.junctions_MappingInfo = junctions_MappingInfo;
    simplify.junctions = junctions;
    simplify.build_node_set_and_simplify(exon_l,exon_r,exon_cov,edge_out,edge_in,edge_weight,Single_nodes,Unused_pair_paths,invalid_nodes);
    simplify.get_graph(Vec_edges,Vec_weights,Unused_pair_paths,Unused_pair_paths_cov);
    simplify.get_triplet_map(Unused_pair_paths,Unused_pair_paths_cov);
    simplify.packing();
    simplify.get_packing_result(Vec_edges);
    for(size_t i=0;i<Node_seq_1.size();i++)
    {
	if(find(invalid_nodes.begin(),invalid_nodes.end(),i) != invalid_nodes.end()) Node_seq_4[i] = 0.00;
    }
    /*
    if(!Vec_edges.empty())
    {
      out_graph<<"Edges"<<endl;
      for(int i=0;i<Vec_edges.size();i++) {
	      out_graph<<Vec_edges[i][0]<<" -> "<<Vec_edges[i][1]<<" : "<<Vec_weights[i]<<" (";
	      int n1 = Vec_edges[i][0], n2 = Vec_edges[i][1];
	      out_graph<<"junction: "<<Node_seq_3[n1]<<" "<<Node_seq_2[n2]<<")"<<endl;
      }
      out_graph<<"Nodes"<<endl;
      for(int i=0;i<Node_seq_1.size();i++) out_graph<<" "<<strand<<" "<<Node_seq_1[i]<<" "<<Node_seq_2[i]<<" "<<Node_seq_3[i]<<" "<<Node_seq_4[i]<<i<<endl;
      out_graph<<"Pair"<<endl;
      for(size_t i=0;i<Unused_pair_paths.size();i++) {
         vector<int> pp = Unused_pair_paths[i];
         for(size_t j=0;j<pp.size();j++) out_graph<<pp[j]<<" ";
         out_graph<<": "<<Unused_pair_paths_cov[i]<<endl;
       }
       out_graph<<"#Graph "<<rg_index<<endl;
    }
    */
    //Get_Junction_Graph junction_graph(Vec_edges,Vec_weights,Unused_pair_paths);
    //junction_graph.Construct_Junction_Graph();
    //simplify.get_packing_result_new(Vec_edges,junction_graph.Edges_left,junction_graph.Edges_right,junction_graph.Weights);
    simplify.path_search(strand,chr);
  }
};


void load_graph(char* file)
{
    vector<int> edge_out,edge_in,exon_l,exon_r;
    vector<double> edge_weight,exon_cov;

    vector<vector<int> > Vec_edges;
    vector<double> Vec_weights;

    vector<vector<int> > Unused_pair_paths;
    vector<double> Unused_pair_paths_cov;

    int line, strand, XS_plus, XS_minus;
    string strand_,chr_;
    vector<string> Node_seq_1;
    vector<int>Node_seq_2,Node_seq_3;
    vector<double>Node_seq_4;


    cerr<<"Loading graph.."<<endl;
    ifstream in(file);
    istringstream istr;
    string s,temp;
    int i=0;
    bool edge_flag = false, node_flag = false, pair_flag = false;
    while(getline(in,s))
    {
        i++;
        if(i % 100000 == 0) cerr<<"Loading "<<i<<" lines"<<endl;
	if( s == "Edges"){ edge_flag = true; node_flag = false;pair_flag = false;continue;}
	if( s == "Nodes") { edge_flag = false; node_flag = true;pair_flag = false;continue;} 
	if( s == "Pair") {edge_flag = false; node_flag = false;pair_flag = true;continue;}
	if( s[0] == '#') {
	    istr.str(s);
	    istr>>temp>>rg_index;
	    if(rg_index != graph_number)
            {
                trans_id = 1;
                graph_number = rg_index;
            }
	    istr.clear();
	    PathSearch ps(edge_out, edge_in, edge_weight,
             		  Unused_pair_paths, Unused_pair_paths_cov,
             		  Node_seq_1, Node_seq_2, Node_seq_3,  Node_seq_4,
			  strand_,chr_);
            		  //line,strand,XS_plus,XS_minus);
	    ps.path_search();

	    edge_flag = false; node_flag = false;pair_flag = false;continue;
	}

        if(edge_flag)
        {
	    istr.str(s);
	    int out,in;
	    vector<int> v;
	    double cov;        
	    istr>>out>>temp>>in>>temp>>cov;
	    istr.clear();
	    edge_out.push_back(out); edge_in.push_back(in);
	    edge_weight.push_back(cov);

	    //v.push_back(out); v.push_back(in);
	    //Vec_edges.push_back(v);
	    //Vec_weights.push_back(cov);

        } 
        if(node_flag)
        {
	    istr.str(s);
	    string chr,strd;
	    int left, right;
	    double cov;
	    istr>>strd>>chr>>left>>right>>cov;
	    chr_ = chr;
	    istr.clear();
	    strand_ = strd;
	    //if(strd == "+") strand = 2;
	    //else strand = 1;
	    //exon_l.push_back(left); exon_r.push_back(right);exon_cov.push_back(cov);
	    Node_seq_1.push_back(chr); Node_seq_2.push_back(left); Node_seq_3.push_back(right); Node_seq_4.push_back(cov);

	}	
	if(pair_flag)
	{
	    istr.str(s);
	    double cov;
	    vector<int> v;
	    while(istr>>temp)
	    {
	        if(temp == ":"){
		    istr>>cov;
		    break;
	 	}
		else
		    v.push_back(atoi(temp.c_str()) );
	    }
	    istr.clear();
	    Unused_pair_paths.push_back(v);
	    Unused_pair_paths_cov.push_back(cov);    
	}
    }



    return;
}
#endif
