#ifndef _NETWORK_H_
#define _NETWORK_H_

#include <map>
#include <deque>
#include <set>
#include <vector>
#include <stack>
#include <queue>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "data_index.h"
//#include "lapackpp.h"

namespace network{
     using namespace std;
     using namespace pattern;
     typedef unsigned long uint64;
     typedef map<uint64, double> OutEdge;
     typedef map<uint64, OutEdge> Graph;
     typedef map<uint64, int> Partition;
     const double PI =3.1415926;
     const double MAX_DOUBLE = 1e20;
     const double MIN_DOUBLE = 1e-20;

	struct Edge
	{
	    Edge(){}
		Edge(uint64 a, uint64 b, double w)
			:v1(a),v2(b),wt(w)
		{}
		uint64 v1,v2;
		double wt;
		double realwt;
	};

	bool operator > (const Edge & ea, const Edge & eb);
	bool operator < (const Edge & ea, const Edge & eb);

	bool operator > (const pair<uint64,double> & ea, const pair<uint64,double> & eb);
	bool operator < (const pair<uint64,double> & ea, const pair<uint64,double> & eb);

	template< typename T > void normalize1(map<T,double> & data){
        double s = 0;
        for(typename map<T,double>::iterator i = data.begin(); i != data.end(); ++i) s += abs(i->second);
        for(typename map<T,double>::iterator i = data.begin(); i != data.end(); ++i) i->second /= s;
	}

     template<typename T1, typename T2>
     bool GreaterSecond(const pair<T1,T2> & a, const pair<T1,T2> & b){
          return a.second > b.second;
     }
     template<typename T1, typename T2>
     bool LessSecond(const pair<T1,T2> & a, const pair<T1,T2> & b){
          return a.second < b.second;
     }

     //random number
     double rand01();
     double powMean(double x1, double x2, double p);
     double edgeSim(OutEdge & a, OutEdge & b);

     //graph load and save
     void loadGraph(const string & name, Graph & G);
     void loadSimpleGraph(const string & name, Graph & G);
     void saveGraphNet(const string & name, const Graph & G);
     void saveGraph(const string & name, const Graph & G);

     //random graphs
     void randomGraph(Graph & G, int n, double d);
     void randomGeometricGraph(Graph & G, int n, double d);
     void randomPowLaw(Graph & G, double c, double beta);

     int edges(const Graph & G);
     void degreeDistribution(const Graph & G, map<int,int> & dgs);
     void powerlaw(const Graph & G, double & c, double & beta);

     //basic
     void deleteVertex(Graph & G, uint64 v);
     int vertices(const Graph & G);
     void vertexSet(const Graph & G, set<uint64> & V);
     void symmetricGraph(const Graph & T, Graph & ST);
     void addweight(Graph & G, uint64 i, uint64 j, double w);
     double haveEdge(const Graph & G, uint64 u, uint64 v);
     void BFSVertices(const Graph & G, uint64 root, int depth, std::set<uint64> & vs, double ratio = 1,double minw=0);
     int component(const Graph & G);
     void subgraph(const Graph & G, std::set<uint64> & S, Graph & sub);
     int triangle(const Graph & G, uint64 v);
     int triple(const Graph & G, uint64 v);
     double clusterCoefficient(const Graph & G, uint64 v);
     double clusterCoefficient(const Graph & G);
     void shortestPathBFS(const Graph & G, uint64 root, OutEdge & path);
     double diameterBFS(const Graph & G);
     void dijkstra(const Graph & G, uint64 root, OutEdge & path);
     void eigen(const Graph & G, vector< double > & egs);
     double heterogeneity(const Graph & G, uint64 v);

     //partition
     double edgecut(const Graph & G, Partition & P, double b = 5);
     void init_partition(const Graph & G, Partition & P, map<uint64,int> vw);
     void kl_partition(Graph & G, map<uint64,int> & T, map<uint64,int> & vw, double balance, bool refine);
     void coarse(Graph & G, map<uint64, set<uint64> > & groups, map<uint64,int> & vw);
	 void coarseN(Graph & G, map<uint64, set<uint64> > & group, map<uint64,int> & vw, int N);
     void kl_partition_loop(Graph & G, Partition & T, map<uint64,int> & vw, double b = 5);
     void multiLevelPartition(Graph & G, Partition & T, map<uint64,int> & vw, double b = 5);
	 Partition gapartition(Graph & G, vector< Partition > & pop);
	 void repartition(const Graph & G, Partition & P);
	 void cross_partition1(const Graph & G, Partition & A, Partition & B);
	 void cross_partition2(const Graph & G, Partition & A, Partition & B);

     void bisection(Graph & G, Graph & G1, Graph & G2);
     void partition(Graph & G, Partition & T, int depth);

     void kmeanGraph(const Graph & G, int K, Partition & P);
     void meanCluster(const Graph & G, Partition & P, double pw);
     void meanShift(Graph & G);

     //rank
     void pagerank(const Graph & G, map<uint64,double> & rank);
     void betweennessCentrality(const Graph & G, map<uint64, double> & cb);
     void indegreerank(const Graph & G, map<uint64,double> & rank);

     //community MCL
     void markovGraph(const Graph & G, Graph & T);
     Graph mcl_expansion2(Graph & T0, Graph & T1);
     void mcl_inflation(Graph & T, double r);
     void mcl_attractor(const Graph & T, vector<uint64> & A);
     void mcl(const Graph & G, Partition & P, double r);
     void mcl_part(const Graph & G, const Graph & T, Partition & P);

	 //layout
	 void layoutDim1(const Graph & G, map<uint64,double> & x);

	 //mst
	 void cutMST(const Graph & G, Partition & P);
	 void primMST(const Graph & G, Graph & T);
	 void primRandomMST(const Graph & G, Graph & T);
	 void randomST(const Graph & G, Graph & T);
	 void sampleByST(const Graph & G, Graph & S, int N = 3);
	 void sampleByRandomMST(const Graph & G, Graph & S, int N = 3);
	 void sampleRandomly(const Graph & G, Graph & S, double p);
     void relevance(const Graph & G, uint64 vq, map<uint64,double> & rel, int N);

    struct PointSimFunc{
        double operator ()(const pair<double,double> & a, const pair<double,double> & b)
        {
            double d = (a.first - b.first) * (a.first - b.first) + (a.second - b.second) * (a.second - b.second);
            d = sqrt(d);
            return exp(-0.01 * d);
        }
    };

    struct VectorSimFunc{
        double operator ()(const vector< double > & a, const vector< double > & b){
            double d = 0;
            int N = a.size();
            for(int i = 0; i < N; ++i){
                d += (a[i] - b[i]) * (a[i] - b[i]);
            }
            d /= (double)(N);
            return exp(-2 * d);
        }
    };

template< int K >void nearestK(const Graph & G, uint64 root, set<uint64> & vs){
		Graph::const_iterator p = G.find(root);
		OutEdge oe = p->second;
		for(OutEdge::const_iterator j = p->second.begin(); j != p->second.end(); ++j){
				Graph::const_iterator pj = G.find(j->first);
				for(OutEdge::const_iterator i = pj->second.begin(); i != pj->second.end(); ++i){
						if(oe.find(i->first) != oe.end()) continue;
						oe[i->first] = sqrt(j->second * i->second);
				}
		}
		vector< pair<uint64,double> > voe(oe.begin(), oe.end());
		sort(voe.begin(), voe.end(), GreaterSecond<uint64,double>);
		for(int i = 0; i < voe.size() && i < K; ++i){
			if(voe[i].second * 2 > voe[0].second)
				vs.insert(voe[i].first);
		}
}


    template < typename SampleType, typename SimFunc >
    void relationSampling2(vector< SampleType > & dataset, Graph & G)
    {
        ofstream out("build.txt");
        G.clear();
        int N = dataset.size();
        vector< int > p(N);
        for(int i = 0; i < N; ++i){
            p[i] = i;
            uint64 ii = (uint64)(i);
            G[ii][ii] = 1;
        }

        SimFunc sim;
        double NN = N;
        NN = 2 * sqrt(NN);
        double maxew = 0;
        double minew = 1;
        for(int g = 0; g < (int)(NN); ++g){
            random_shuffle(p.begin(), p.end());
            for(int i = 1; i < p.size(); ++i){
                double s = sim(dataset[p[i-1]], dataset[p[i]]);
                uint64 ha = (uint64)(p[i-1]);
                uint64 hb = (uint64)(p[i]);
                G[ha][hb] = s;
                G[hb][ha] = s;
                maxew = max<double>(maxew, s);
                minew = min<double>(minew, s);
            }
        }
        out << edges(G) << endl;
        while(true){
            bool newedge = false;
            for(Graph::iterator i = G.begin(); i != G.end(); ++i){
                set<uint64> oe;
                nearestK<8>(G,i->first,oe);
                for(set<uint64>::iterator a = oe.begin(); a != oe.end(); ++a){
                    for(set<uint64>::iterator b = a; b != oe.end(); ++b){
                    		if(a == b) continue;
                        if(haveEdge(G, *a, *b)) continue;
                        int aa = (int)(*a);
                        int bb = (int)(*b);
                        double s = sim(dataset[aa],dataset[bb]);
                        G[*a][*b] = s;
                        G[*b][*a] = s;
                        maxew = max<double>(maxew, s);
                				minew = min<double>(minew, s);
                        newedge = true;
                    }
                }
            }
            out << edges(G) << endl;
            if(newedge == false) break;
        }
        for(Graph::iterator i = G.begin(); i != G.end(); ++i){
            		int findte = 0;
            		for(OutEdge::iterator j = i->second.begin(); j != i->second.end(); ++j){
            				if(j->second > maxew * 0.7){
            					findte++;
            				}
            		}
            		if(findte < 2){
            				for(int j = 0; j < N; ++j){
            						uint64 hj = (uint64)(j);
            						if(i->first == hj) continue;
            						if(i->second.find(hj) != i->second.end()) continue;
            						double s = sim(dataset[i->first], dataset[j]);
            						i->second[hj] = s;
            						G[hj][i->first] = s;
            				}
            		}
            }
            out << edges(G) << endl;
        Graph Gt = G;
        sampleByRandomMST(Gt,G,5);
        meanShift(G);
        out.close();
    }

    template < typename SampleType, typename SimFunc >
    void relationSampling(vector< SampleType > & dataset, Graph & G)
    {
        Cell2D< pair<double, double> > cell;
        double K = dataset.size();
        K = sqrt(K / 4);
        cell.N = (int)(K) + 1;
        buildCell2D(dataset, cell);
        SimFunc sim;
        for(int i = 0; i < dataset.size(); ++i){
            set<int> rp;
            searchCell2D(cell, dataset[i], 1, rp);
            for(set<int>::iterator a = rp.begin(); a != rp.end(); ++a){
                for(set<int>::iterator b = a; b != rp.end(); ++b){
                    if(a == b) continue;
                    double s = sim(dataset[*a],dataset[*b]);
                    uint64 ha = *a;
                    uint64 hb = *b;
                    G[ha][hb] = s;
                    G[hb][ha] = s;
                }
            }
        }
        Graph Gt = G;
        sampleByRandomMST(Gt,G,5);
        /*
        for(int g = 0; g < 1; ++g){
            Gt.clear();
            for(Graph::iterator i = G.begin(); i != G.end(); ++i){
                OutEdge oei = i->second;
                oei[i->first] = 1;
                for(OutEdge::iterator j = i->second.begin(); j != i->second.end(); ++j){
                    if(i->first >= j->first) continue;
                    OutEdge oej = G[j->first];
                    oej[j->first] = 1;
                    double s = edgeSim(oei,oej);
                    Gt[i->first][j->first] = s;
                    Gt[j->first][i->first] = s;
                }
            }
            G = Gt;
        }
        */
        //meanShift(G);
    }
};

#endif
