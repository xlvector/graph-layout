#include "network.h"

namespace network{
    double rand01()
    {
          double r = (double)(rand() % 10000);
          return r / 10000.0;
     }

     double edgeSim(OutEdge & a, OutEdge & b){
        set<uint64> nodes;
        for(OutEdge::iterator i = a.begin(); i != a.end(); ++i){
            nodes.insert(i->first);
        }
        for(OutEdge::iterator i = b.begin(); i != b.end(); ++i){
            nodes.insert(i->first);
        }
        double d = 0;
        for(set<uint64>::iterator i = nodes.begin(); i != nodes.end(); ++i){
            double da = 0;
            double db = 0;
            if(a.find(*i) != a.end()) da = a[*i];
            if(b.find(*i) != b.end()) db = b[*i];
            d += abs(da - db);
        }
        d /= (double)(nodes.size());
        return exp(-3 * d);
     }
double powMean(double x1, double x2, double p){
    	return pow((pow(x1,p) + pow(x2,p)) * 0.5, 1/p);
    }
    bool operator > (const Edge & ea, const Edge & eb){
		return ea.wt > eb.wt;
    }

	bool operator < (const pair<uint64,double> & ea, const pair<uint64,double> & eb){
		return ea.second < eb.second;
	}

    bool operator > (const pair<uint64,double> & ea, const pair<uint64,double> & eb){
		return ea.second > eb.second;
    }

	bool operator < (const Edge & ea, const Edge & eb){
		return ea.wt < eb.wt;
	}

     void loadGraph(const string & name, Graph & G)
     {
          ifstream in(name.c_str());
          uint64 i,j;
          double w;
          while(in >> i >> j >> w)
          {
               G[i][j] = w;
               G[j][i] = w;
          }
          in.close();
     }

     void loadSimpleGraph(const string & name, Graph & G)
     {
          ifstream in(name.c_str());
          uint64 i,j;
          while(in >> i >> j)
          {
               G[i][j] = 1;
               G[j][i] = 1;
          }
          in.close();
     }

     void saveGraph(const string & name, const Graph & G)
     {
          ofstream out(name.c_str());
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i)
          {
               for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j)
               {
                    out << i->first << "\t" << j->first << "\t" << j->second << endl;
               }
          }
          out.close();
     }

     void saveGraphNet(const string & name, const Graph & G)
     {
          int n = 1;
          map<uint64,int> vindex;
          ofstream out(name.c_str());
          out << "vertices " << G.size() << endl;
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i)
          {
               vindex[i->first] = n;
               out << n << "  " << i->first << endl;
               n++;
          }
          out << "edges" << endl;
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i)
          {
               for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j)
               {
                    out << vindex[i->first] << "\t" << vindex[j->first] << "\t" << j->second << endl;
               }
          }
     }

     int edges(const Graph & G)
     {
          int ret = 0;
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
               ret += i->second.size();
          }
          return ret>>1;
     }

     void deleteVertex(Graph & G, uint64 v){
          G.erase(v);
          vector<uint64> ep;
          for(Graph::iterator i = G.begin(); i != G.end(); ++i){
               i->second.erase(v);
               if(i->second.empty()) ep.push_back(i->first);
          }
          for(int i = 0; i < ep.size(); ++i) G.erase(ep[i]);
     }

	int triangle(const Graph & G, uint64 v){
		 Graph::const_iterator p = G.find(v);
		 int ret = 0;
		 for(OutEdge::const_iterator i = p->second.begin(); i != p->second.end(); ++i){
			  for(OutEdge::const_iterator j = i; j != p->second.end(); ++j){
				   if(i == j) continue;
				   if(haveEdge(G,i->first,j->first) != 0) ret++;
			  }
		 }
		 return ret;
	}

	int triple(const Graph & G, uint64 v){
		 Graph::const_iterator p = G.find(v);
		 int ret = p->second.size();
		 return ((ret * ret) >> 1) - (ret >> 1);
	}

    double clusterCoefficient(const Graph & G, uint64 v)
    {
        return (double)(triangle(G,v)) / (double)(triple(G,v));
    }

    double clusterCoefficient(const Graph & G)
    {
        double n = 0;
        double ret = 0;
        for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
            if(i->second.size() < 2) continue;
            ret += clusterCoefficient(G,i->first);
            n += 1;
        }
        return ret / n;
    }

     void degreeDistribution(const Graph & G, map<int,int> & dgs){
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i) dgs[i->second.size()]++;
     }

     void powerlaw(const Graph & G, double & c, double & beta)
     {
          map<int,int> dgs;
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i) dgs[i->second.size()]++;
          double a11 = 0;
          double a12 = 0;
          double a22 = 0;
          double b1 = 0;
          double b2 = 0;
          for(map<int,int>::iterator i = dgs.begin(); i != dgs.end(); ++i){
               double k = (double)(i->first);
               double n = (double)(i->second);
               k = log10(k);
               n = log10(n);
               a11 += k * k;
               a12 += k;
               a22 += 1;
               b1 += k * n;
               b2 += n;
          }
          double x = (b1 * a22 - b2 * a12)/(a11 * a22 - a12 * a12);
          double y = (b2 * a11 - b1 * a12)/(a11 * a22 - a12 * a12);
          beta = -1 * x;
          c = pow(10,y);
     }

     void randomGraph(Graph & G, int n, double d)
     {
          srand(1985718);
          double p = d / (double)(n - 1);
          double e = 0;
          for(int i = 0; i < n; ++i){
               for(int j = i + 1; j < n; ++j){
                    if(rand01() < p){
                         G[i][j] = 1;
                         G[j][i] = 1;
                         ++e;
                    }
               }
          }
          cout << G.size() << " : " << e << endl;
     }

     void randomGeometricGraph(Graph & G, int n, double d)
     {
          srand(1985718);
          double t = d / (double)(n);
          t = t / PI;
          double x[n],y[n];
          for(int i = 0; i < n; ++i){
               x[i] = rand01();
               y[i] = rand01();
          }
          double e = 0;
          for(int i = 0; i < n; ++i){
               for(int j = 0; j < n; ++j){
                    double d = (x[i] - x[j]) * (x[i] - x[j]) + (y[i] - y[j]) * (y[i] - y[j]);
                    if(d > t) continue;
                    G[i][j] = 1;
                    G[j][i] = 1;
                    ++e;
               }
          }
          cout << G.size() << " : " << e << endl;
     }

     void randomPowLaw(Graph & G, double c, double beta)
     {
          vector< int > dgs;
          double kk;
          for(int k = 1;; ++k){
               kk = (double)(k);
               double kk2 = pow(kk,beta);
               int n = (int)(c/kk2);
               if(n < 1) break;
               for(int i = 0; i < n; ++i) dgs.push_back(k);
          }
          cout << "max degree : " << kk << endl;
          int N = dgs.size();
          for(uint64 i = 0; i < N; ++i){
               int d = dgs[i];
               for(uint64 j = 0; j < d; ++j){
                    uint64 k = (uint64)(rand() % N);
                    G[i][k] = 1;
                    G[k][i] = 1;
               }
          }
     }


     double haveEdge(const Graph & G, uint64 u, uint64 v)
     {
          Graph::const_iterator oe = G.find(u);
          if(oe == G.end()) return 0;
          if(oe->second.find(v) == oe->second.end()) return 0;
          return 1;
     }

     void addweight(Graph & G, uint64 i, uint64 j, double w)
     {
          if(w < MIN_DOUBLE) return;
          Graph::iterator pg = G.find(i);
          if(pg == G.end()){
               G[i][j] = w;
               return;
          }
          OutEdge::iterator po = pg->second.find(j);
          if(po == pg->second.end()){
               pg->second[j] = w;
               return;
          }
          else{
               pg->second[j] += w;
               return;
          }
     }

     void dijkstra(const Graph & G, uint64 root, OutEdge & path)
     {
          double T = 3;
          map<double, set<uint64> > pQ;
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
               double d = MAX_DOUBLE;
               if(i->first == root) d = 0;
               pQ[d].insert(i->first);
               path[i->first] = d;
          }
          while(!pQ.empty()){
               map<double, set<uint64> >::iterator top = pQ.begin();
               uint64 v = *(top->second.begin());
               double minw = top->first;
               if(minw >= T) break;
               //cout << pQ.size() << "\t" << v << "\t" << minw << endl;
               top->second.erase(v);
               if(top->second.empty()) pQ.erase(top);

               Graph::const_iterator oe = G.find(v);
               OutEdge oev = oe->second;
               for(OutEdge::iterator k = oev.begin(); k != oev.end(); ++k){
                    double pp = minw + k->second;
                    double pk = path[k->first];
                    if(pp < pk){
                         map<double, set<uint64> >::iterator ppk = pQ.find(pk);
                         ppk->second.erase(k->first);
                         if(ppk->second.empty()) pQ.erase(ppk);
                         pQ[pp].insert(k->first);
                         path[k->first] = pp;
                    }
               }
          }
          //cout << endl;
          OutEdge tmp = path;
          path.clear();
          for(OutEdge::iterator i = tmp.begin(); i != tmp.end(); ++i)
               if(i->second < T)
                    path[i->first] = i->second;
     }

     void shortestPathBFS(const Graph & G, uint64 root, OutEdge & path)
     {
          map<uint64,int> nodedepth;
          nodedepth[root] = 0;
          deque<uint64> visited;
          visited.push_back(root);
          while(!visited.empty())
          {
               uint64 v = visited.front();
               visited.pop_front();
               int dp = nodedepth[v];
               Graph::const_iterator oe = G.find(v);
               if(oe == G.end()) continue;
               for(OutEdge::const_iterator k = oe->second.begin(); k != oe->second.end(); ++k)
               {
                    map<uint64,int>::iterator pm = nodedepth.find(k->first);
                    if(pm == nodedepth.end())
                    {
                         nodedepth[k->first] = dp + 1;
                         visited.push_back(k->first);
                    }
               }
          }
          for(map<uint64,int>::iterator i = nodedepth.begin(); i != nodedepth.end(); ++i){
               if(i->second == 0) continue;
               path[i->first] = (double)(i->second);
          }
     }

    double diameterBFS(const Graph & G)
    {
        double ret = 0;
        for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
            OutEdge path;
            shortestPathBFS(G,i->first,path);
            for(OutEdge::iterator j = path.begin(); j != path.end(); ++j) ret = max<double>(ret,j->second);
        }
        return ret;
    }

     void BFSVertices(const Graph & G, uint64 root, int depth, std::set<uint64> & vs, double ratio,double minw)
     {
          vs.clear();
          int maxsize = (int)(G.size() * ratio);
          map<uint64,int> nodedepth;
          nodedepth[root] = 0;
          deque<uint64> visited;
          visited.push_back(root);
          while(!visited.empty() && vs.size() < maxsize)
          {
               uint64 v = visited[0];
               visited.pop_front();
               int dp = nodedepth[v];
               if(dp >= depth) break;
               vs.insert(v);
               Graph::const_iterator oe = G.find(v);
               if(oe == G.end()) continue;
               for(OutEdge::const_iterator k = oe->second.begin(); k != oe->second.end(); ++k)
               {
                    if(k->second < minw) continue;
                    map<uint64,int>::iterator pm = nodedepth.find(k->first);
                    if(pm == nodedepth.end())
                    {
                         nodedepth[k->first] = dp + 1;
                         visited.push_back(k->first);
                    }
               }
          }
     }

     int component(const Graph & G){
          int ret = 0;
          set<uint64> V;
          vertexSet(G,V);
          while(!V.empty()){
               deque<uint64> Q;
               set<uint64> visited;
               uint64 root = *(V.begin());
               Q.push_back(root);
               while(!Q.empty()){
                    uint64 v = Q.front();
                    Q.pop_front();
                    if(visited.find(v) != visited.end()) continue;
                    visited.insert(v);
                    V.erase(v);
                    Graph::const_iterator oe = G.find(v);
                    if(oe == G.end()) continue;
                    for(OutEdge::const_iterator j = oe->second.begin(); j != oe->second.end(); ++j){
                         Q.push_back(j->first);
                    }
               }
               ++ret;
          }
          return ret;
     }

    void meanCluster(const Graph & G, Partition & P, double pw){
        ofstream out("mean.txt");
        int comp = 0;
        set<uint64> V;
        vertexSet(G,V);
        map<uint64, double> dm;
        map<uint64, int> dg;
        for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
        		dg[i->first] = i->second.size();
            double d = 0;
            for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j){
                d += pow(j->second,pw);
            }
            double nn = (double)(i->second.size());
            if(nn == 0) d = 0;
            else d = pow((d / nn), 1/pw);
            dm[i->first] = d;
            out << i->first << "\t" << d << endl;
        }
        while(!V.empty()){
            deque<uint64> Q;
            set<uint64> visited;
            uint64 root = *(V.begin());
            Q.push_back(root);
            while(!Q.empty()){
                uint64 v = Q.front();
                Q.pop_front();
                if(visited.find(v) != visited.end()) continue;
                P[v] = comp;
                visited.insert(v);
                V.erase(v);
                Graph::const_iterator oe = G.find(v);
                if(oe == G.end()) continue;
                for(OutEdge::const_iterator j = oe->second.begin(); j != oe->second.end(); ++j){
                    if(dg[j->first] > 1 && (j->second < powMean(dm[v],dm[j->first], 2))) continue;
                    Q.push_back(j->first);
                }
            }
            ++comp;
        }
    }

     void subgraph(const Graph & G, std::set<uint64> & S, Graph & sub)
     {
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i)
          {
               if(S.find(i->first) == S.end()) continue;
               for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j)
               {
                    if(S.find(j->first) == S.end()) continue;
                    sub[i->first][j->first] = j->second;
                    //sub[j->first][i->first] = j->second;
                    //cout << i->first << "\t" << j->first << "\t" << j->second << endl;
               }
          }
     }

     double edgecut(const Graph & G, Partition & P, double cur_balance)
     {
          double s_opt = (double)(G.size()) / 2;
          double s_p = s_opt * (100 + cur_balance) / 100;
          double s[2] = {0,0};
          double ret = 0;
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
               int ci = P[i->first];
               s[ci]++;
               for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j){
                    int cj = P[j->first];
                    if(ci != cj) ret += j->second;
               }
          }
          if(s[0] > s_p || s[1] > s_p) return 1E20;
          return ret/2;
     }


     void init_partition(const Graph & G, Partition & P, map<uint64,int> vw)
     {
          vector<uint64> nodes;
          int all = 0;
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
               nodes.push_back(i->first);
               P[i->first] = -1;
               all += vw[i->first];
          }
          random_shuffle(nodes.begin(), nodes.end());
          P[nodes[0]] = 0;
          int s0 = 1;

          for(int g = 0; g < G.size(); ++g){
               for(int i = 0; i < nodes.size(); ++i){
                    uint64 v = nodes[i];
                    int vtw = vw[v];
                    if((s0 + vtw) * 2 > all && s0*2 < all) continue;
                    int c = P[v];
                    if(c >= 0) continue;
                    Graph::const_iterator oe = G.find(v);
                    double ga0 = 0;
                    for(OutEdge::const_iterator j = oe->second.begin(); j != oe->second.end(); ++j){
                         int cj = P[j->first];
                         if(cj == 0) ga0++;
                         if(ga0 > 0) continue;
                    }
                    if(ga0 > 0){
                         P[v] = 0;
                         s0 += vw[v];
                    }
                    if(s0 * 2 > all) break;
               }
               if(s0 * 2 >= all) break;
          }
          for(Partition::iterator i = P.begin(); i != P.end(); ++i){
               if(i->second < 0) i->second = 1;
          }
     }

     void kl_partition(Graph & G, map<uint64,int> & T, map<uint64,int> & vw, double balance, bool refine)
     {
          double s_opt = 0;
          for(map<uint64,int>::iterator i = vw.begin(); i != vw.end(); ++i) s_opt += (double)(i->second);
          s_opt = s_opt / 2;
          double s_p = s_opt * (100 + balance) / 100;
          if(!refine) init_partition(G,T,vw);
          map<uint64,double> gain;
          double s[2] = {0,0};
          for(Graph::iterator i = G.begin(); i != G.end(); ++i){
               double ga = 0;
               int ci = T[i->first];
               s[ci] += vw[i->first];
               for(OutEdge::iterator j = i->second.begin(); j != i->second.end(); ++j){
                    int cj = T[j->first];
                    if(ci == cj) ga -= j->second;
                    else ga += j->second;
               }
               gain[i->first] = ga;
          }

          for(int g = 0; g < 100; ++g){
               bool find = false;
               double tg = max<double>(0,(double)(8 - g));
               for(map<uint64,double>::iterator i = gain.begin(); i != gain.end(); ++i){
                    int ci = T[i->first];
                    double vtw = (double)(vw[i->first]);
                    if(i->second > tg){
                         if(s[1-ci] + vtw >= s_p) continue;
                         find = true;
                         OutEdge & oe = G[i->first];
                         for(OutEdge::iterator j = oe.begin(); j != oe.end(); ++j){
                              int cj = T[j->first];
                              if(ci == cj) gain[j->first] += 2 * j->second;
                              else gain[j->first] -= 2 * j->second;
                         }
                         T[i->first] = 1 - ci;
                         i->second *= -1;
                         s[ci] -= vtw;
                         s[1-ci] += vtw;
                    }
               }
               if(!find && tg == 0) break;
          }
     }
/*
	void matching(const Graph & G, vector< Edge > & mc)
	{
		map<uint64,int> nodes;
		Grapj::const_iterator p = G.begin();
		uint64 u = p->first;
		nodes[u] = 0;
		deque<uint64> Q;
		Q.push_back(u);
		while(!Q.empty()){

		}
	}
*/
     void coarse(Graph & G, map<uint64, set<uint64> > & groups, map<uint64,int> & vw)
     {
          vector< pair<uint64,double> > vs;
          for(map<uint64,int>::iterator i = vw.begin(); i != vw.end(); ++i){
               pair<uint64,double> vsp(i->first, (double)(i->second) + rand01());
               vs.push_back(vsp);
          }
          random_shuffle(vs.begin(), vs.end());
          //sort(vs.begin(), vs.end(), GreaterSecond<uint,double>);
          //reverse(vs.begin(), vs.end());
          map<uint64,uint64> nodemap;
          for(int i = 0; i < vs.size(); ++i)
          {
               uint64 v = vs[i].first;
               if(nodemap.find(v) != nodemap.end()) continue;
               nodemap[v] = v;
               groups[v].insert(v);
               map<uint64, set<uint64> >::iterator p = groups.find(v);

               OutEdge & oe = G[v];
               int n = 0;
               double maxe = 0;
               uint64 maxv;
               for(OutEdge::iterator j = oe.begin(); j != oe.end(); ++j){
                    if(nodemap.find(j->first) != nodemap.end()) continue;
                    if(maxe < j->second){
                         maxe = j->second;
                         maxv = j->first;
                    }
                    if(maxe == j->second && rand() % 2 == 0) maxv = j->first;
               }
               if(maxe > 0){
                    p->second.insert(maxv);
                    nodemap[maxv] = v;
               }
          }
          groups.clear();
          for(map<uint64,uint64>::iterator i = nodemap.begin(); i != nodemap.end(); ++i) groups[i->second].insert(i->first);
          Graph ret;
          map<uint64,int> tvw = vw;
          vw.clear();
          for(Graph::iterator i = G.begin(); i != G.end(); ++i){
               uint64 ii = nodemap[i->first];
               vw[ii] += tvw[i->first];
               Graph::iterator pr = ret.find(ii);
               if(pr == ret.end()){
                    OutEdge oe;
                    ret[ii] = oe;
                    pr = ret.find(ii);
               }
               for(OutEdge::iterator j = i->second.begin(); j != i->second.end(); ++j){
                    uint64 jj = nodemap[j->first];
                    if(ii == jj) continue;
                    if(pr->second.find(jj) == pr->second.end()) pr->second[jj] = j->second;
                    else pr->second[jj] += j->second;
               }
          }
          G = ret;
     }

		void coarseN(Graph & G, map<uint64, set<uint64> > & group, map<uint64,int> & vw, int N)
		{
			 int n = G.size();
			 n = n >> N;
			 while(G.size() > N){
				  cout << G.size() << endl;
				  map<uint64, set<uint64> > g;
				  int ng = G.size();
				  coarse(G,g,vw);
				  if(G.size() == ng) break;
				  for(map<uint64,set<uint64> >::iterator i = g.begin(); i != g.end(); ++i){
					   set<uint64> tmp;
					   for(set<uint64>::iterator j = i->second.begin(); j != i->second.end(); ++j){
							tmp.insert(*j);
							map<uint64, set<uint64> >::iterator pg = group.find(*j);
							if(pg == group.end()) continue;
							for(set<uint64>::iterator k = pg->second.begin(); k != pg->second.end(); ++k)
								 tmp.insert(*k);
					   }
					   i->second = tmp;
				  }
				  group = g;
			 }
		}

     void kl_partition_loop(Graph & G, Partition & T, map<uint64,int> & vw, double cur_balance)
     {
          double mincut = 1E20;
          for(int i = 0; i < 16; ++i){
               Partition tT;
               kl_partition(G,tT,vw,cur_balance,false);
               double cut = edgecut(G,tT);
               if(mincut > cut){
                    mincut = cut;
                    T = tT;
               }
          }
     }

     void multiLevelPartition(Graph & G, Partition & T, map<uint64,int> & vw, double cur_balance)
     {
          if(G.size() < 120 + rand() % 100){
               kl_partition_loop(G,T,vw);
               return;
          }
          map<uint64, set<uint64> > groups;
          map<uint64,int> vw1 = vw;
          Graph G1 = G;
          coarse(G1,groups,vw);
          map<uint64,int> T2;
          multiLevelPartition(G1,T2,vw);
          for(map<uint64,int>::iterator i = T2.begin(); i != T2.end(); ++i)
          {
               set<uint64> & S = groups[i->first];
               for(set<uint64>::iterator k = S.begin(); k != S.end(); ++k) T[*k] = i->second;
          }
          kl_partition(G,T,vw1,cur_balance,true);
     }

	 void repartition(const Graph & G, Partition & P)
     {
          map<uint64,int> vw;
          Graph G1 = G;
          for(Graph::iterator i = G1.begin(); i != G1.end(); ++i)
          {
               vw[i->first] = 1;
               int ci = P[i->first];
               for(OutEdge::iterator j = i->second.begin(); j != i->second.end(); ++j){
                    int cj = P[j->first];
                    if(ci != cj) j->second = 1 - rand01() / 10;
               }
          }
          Partition P1;
          multiLevelPartition(G1,P1,vw);
          double c = edgecut(G,P);
          if(edgecut(G,P1) < c) P = P1;;
     }

	 void cross_partition1(const Graph & G, Partition & A, Partition & B)
     {
          set<uint64> bs0,bs1;
          map<uint64,int> vw;
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
               vw[i->first] = 1;
               int cai = A[i->first];
               int cbi = B[i->first];
               for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j){
                    int caj = A[j->first];
                    int cbj = B[j->first];
                    if(cai != caj){
                         bs0.insert(i->first);
                         bs0.insert(j->first);
                    }
                    if(cbi != cbj){
                         bs1.insert(i->first);
                         bs1.insert(j->first);
                    }
               }
          }
          set<uint64> bs;
          for(set<uint64>::iterator i = bs1.begin(); i != bs1.end(); ++i){
               if(bs0.find(*i) != bs0.end()) bs.insert(*i);
          }
		  Graph G1 = G;
          for(Graph::iterator i = G1.begin(); i != G1.end(); ++i){
               if(bs.find(i->first) != bs.end()) continue;
               for(OutEdge::iterator j = i->second.begin(); j != i->second.end(); ++j){
                    if(bs.find(j->first) != bs.end()) continue;
                    j->second = 2 - rand01() / 10;
               }
          }
          Partition P;
          multiLevelPartition(G1,P,vw);
          double ca = edgecut(G,A);
          if(edgecut(G,P) < ca) A = P;
     }

     void cross_partition2(Graph & G, Partition & A, Partition & B)
     {
          double n = (double)(A.size());
          Partition tA,tB,P;
          for(Partition::iterator i = A.begin(); i != A.end(); ++i){
               tA[i->first] = 0; tB[i->first] = 0;
          }
          for(Partition::iterator i = A.begin(); i != A.end(); ++i){
               int c = B[i->first];
               if(i->second == 0 && c == 0){
                    tA[i->first] = 0;
                    continue;
               }
               if(i->second == 1 && c == 1){
                    tA[i->first] = 1;
                    continue;
               }
               if(i->second == 0 && c == 1){
                    tB[i->first] = 0;
                    continue;
               }
               if(i->second == 1 && c == 0){
                    tB[i->first] = 1;
                    continue;
               }
          }
          Graph G1 = G;
          map<uint64,int> vw;
		  for(Graph::iterator i = G1.begin(); i != G1.end(); ++i)
          {
               vw[i->first] = 1;
               int ci = tA[i->first];
               if(ci < 0) continue;
               for(OutEdge::iterator j = i->second.begin(); j != i->second.end(); ++j){
                    int cj = tA[j->first];
                    if(ci != cj) j->second = 1 - rand01() / 10;
               }
          }
          multiLevelPartition(G1,P,vw);
          double ca = edgecut(G,A);
          if(edgecut(G,P) < ca) A = P;

          G1 = G;
          for(Graph::iterator i = G1.begin(); i != G1.end(); ++i)
          {
               int ci = tB[i->first];
               if(ci < 0) continue;
               for(OutEdge::iterator j = i->second.begin(); j != i->second.end(); ++j){
                    int cj = tB[j->first];
                    if(ci != cj) j->second = 1 - rand01() / 10;
               }
          }
          P.clear();
          multiLevelPartition(G1,P,vw);
          double cb = edgecut(G,B);
          if(edgecut(G,P) < cb) B = P;

     }

	 Partition gapartition(Graph & G, vector< Partition > & pop)
     {
          cout << "begin genetic partition ..." << endl;
          Partition gP;
          double mincut = 1E20;
          vector< pair<int,double> > rank;
          for(int i = 0; i < pop.size(); ++i){
               double c = edgecut(G,pop[i]);
               rank.push_back(make_pair<int,double>(i,c));
               if(mincut > c){
                    mincut = c;
                    gP = pop[i];
               }
          }
		  for(int g = 0; g < 8; ++g){
               random_shuffle(rank.begin(), rank.end());
               for(int i = 0; i < rank.size() - 1; ++i){
                    int a = rank[i].first;
                    int b = rank[i+1].first;
                    double r = rand01();
                    if(r < 0.3)
                         cross_partition1(G,pop[a],pop[b]);
                    if(r < 0.6 && r >= 0.3)
                         cross_partition2(G,pop[a],pop[b]);
                    if(r < 0.8 && r >= 0.6)
                         repartition(G,pop[a]);
                    if(r >= 0.8)
                         repartition(G,pop[b]);
                    rank[i].second = edgecut(G,pop[a]);
                    rank[i+1].second = edgecut(G,pop[b]);
               }
               for(int i = 0; i < rank.size(); ++i){
                    if(mincut > rank[i].second){
                         mincut = rank[i].second;
                         gP = pop[rank[i].first];
                    }
               }
          }
          return gP;
     }

     void bisection(Graph & G, Graph & G1, Graph & G2)
     {
          Partition gP;
          double mincd = 1E20;
          for(int g = 0; g < 8; ++g){
               Partition P;
               map<uint64,int> vw;
               for(Graph::iterator i = G.begin(); i != G.end(); ++i) vw[i->first] = 1;
               map<uint64,int> vw1 = vw;
               multiLevelPartition(G,P,vw);
               double cw = edgecut(G,P);
               if(mincd > cw){
                    mincd = cw;
                    gP = P;
               }
          }
          for(Graph::iterator i = G.begin(); i != G.end(); ++i){
               int ci = gP[i->first];
               for(OutEdge::iterator j = i->second.begin(); j != i->second.end(); ++j){
                    int cj = gP[j->first];
                    if(ci == 0 && cj == 0) G1[i->first][j->first] = j->second;
                    if(ci == 1 && cj == 1) G2[i->first][j->first] = j->second;
               }
          }
     }

     void partition(Graph & G, Partition & T,int depth)
     {
          if(depth > 2) return;
          Graph S1,S2;
          bisection(G,S1,S2);
          cout << S1.size() << "\t" << S2.size() << endl;
          for(Graph::iterator i = S1.begin(); i != S1.end(); ++i)
          {
               int v = T[i->first];
               v = (v << 1);
               T[i->first] = v;
          }
          for(Graph::iterator i = S2.begin(); i != S2.end(); ++i)
          {
               int v = T[i->first];
               v = (v << 1) + 1;
               T[i->first] = v;
          }
          if(S1.size() == 0 || S2.size() == 0) return;
          partition(S1,T,depth+1);
          partition(S2,T,depth+1);
     }

     void loadGraph_graph(const string & name, Graph & G)
     {
          ifstream in(name.c_str());
          int v,e;
          in >> v >> e;
          cout << v << "  " << e << endl;
          string line;
          int es = 0;
          uint64 i = 0;
          while(getline(in,line)){
               istringstream iss(line);
               uint64 node;
               while(iss >> node){
                    G[i][node] = 1;
                    es++;
               }
               i++;
          }
     }

     void symmetricGraph(const Graph & T, Graph & ST)
     {
          for(Graph::const_iterator i = T.begin(); i != T.end(); ++i){
               for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j){
                    ST[j->first][i->first] = j->second;
               }
          }
     }

     void markovGraph(const Graph & G, Graph & T)
     {
          T.clear();
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
               double d = 0;
               for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j){
                    d += j->second;
               }
               if(d < MIN_DOUBLE) continue;
               for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j){
                    double w = j->second / d;
                    T[i->first][j->first] = w;
               }
          }
     }

     Graph mcl_expansion2(Graph & T0, Graph & T1)
     {
          Graph T2;
          for(Graph::iterator i = T0.begin(); i != T0.end(); ++i){
               uint64 p = i->first;
               for(OutEdge::iterator j = i->second.begin(); j != i->second.end(); ++j){
                    uint64 r = j->first;
                    Graph::iterator pr = T1.find(r);
                    if(pr == T1.end()) continue;
                    OutEdge & ro = pr->second;
                    for(OutEdge::iterator k = ro.begin(); k != ro.end(); ++k){
                         uint64 q = k->first;
                         double w = j->second * k ->second;
                         addweight(T2,p,q,w);
                    }
               }
          }
          return T2;
     }
     void mcl_inflation(Graph & T, double r)
     {
          double small = 0.1 / (double)(T.size());
          for(Graph::iterator i = T.begin(); i != T.end(); ++i){
               double d = 0;
               vector< pair<uint64,double> > oe(i->second.begin(), i->second.end());
               sort(oe.begin(), oe.end(), GreaterSecond<uint64,double>);
               OutEdge ooe;
               for(int j = 0; j < oe.size() && j < 8; ++j){
                    double wr = pow(oe[j].second,r);
                    d += wr;
                    ooe[oe[j].first] = wr;
               }

               for(OutEdge::iterator j = ooe.begin(); j != ooe.end(); ++j) j->second /= d;
               i->second = ooe;
          }
     }

     void mcl_attractor(const Graph & T, vector<uint64> & A)
     {
          for(Graph::const_iterator i = T.begin(); i != T.end(); ++i){
               if(i->second.find(i->first) == i->second.end()) continue;
               A.push_back(i->first);
          }
     }


     void mcl_part(const Graph & G, const Graph & T, Partition & P)
     {
          vector<uint64> A;
          mcl_attractor(T,A);
          Graph ST;
          symmetricGraph(T,ST);
          for(int p = 0; p < A.size(); ++p)
          {
               OutEdge & oe = ST[A[p]];
               for(OutEdge::iterator j = oe.begin(); j != oe.end(); ++j){
                    if(j->second > 0.499999) P[j->first] = p;
               }
          }
     }

     void mcl(const Graph & G, Partition & P, double r)
     {
          Graph T;
          markovGraph(G,T);
          cout << T.size() << endl;
          for(int k = 0; k < 30; ++k){
               cout << k << endl;
               T = mcl_expansion2(T,T);
               mcl_inflation(T,r);
          }
          mcl_part(G,T,P);
     }

     void indegreerank(const Graph & G, map<uint64,double> & rank){
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
               for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j){
                    if(rank.find(j->first) == rank.end()) rank[j->first] = 0;
                    rank[j->first] += j->second;
               }
          }
     }

     void pagerank(const Graph & G, map<uint64,double> & rank)
     {
          double N = vertices(G);
          double q = 0.15;
          map<uint64,double> ds;
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
               rank[i->first] = 1;
               double d = 0;
               for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j)
                    d += j->second;
               ds[i->first] = d;
          }
          for(int g = 0; g < 12; ++g){
               for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
                    map<uint64,double>::iterator p = rank.find(i->first);
                    p->second = q/N;
                    for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j){
                         Graph::const_iterator k = G.find(j->first);
                         if(k == G.end()) continue;
                         p->second += (1 - q) * rank[j->first] * j->second / (ds[k->first]);
                    }
               }
          }
     }

	void primMST(const Graph & G, Graph & T)
	{
		priority_queue<Edge> pQ;
		Graph::const_iterator p = G.begin();
		for(OutEdge::const_iterator i = p->second.begin(); i != p->second.end(); ++i){
			Edge e(p->first,i->first,i->second);
			pQ.push(e);
		}
		set<uint64> visited;
		T.clear();
		while(!pQ.empty()){
			Edge e = pQ.top();
			pQ.pop();
			//if(visited.find(e.v1) != visited.end()) continue;
			if(visited.find(e.v2) != visited.end()) continue;
			visited.insert(e.v1);
			visited.insert(e.v2);
			T[e.v1][e.v2] = e.wt;
			T[e.v2][e.v1] = e.wt;
			p = G.find(e.v2);
			for(OutEdge::const_iterator i = p->second.begin(); i != p->second.end(); ++i){
				Edge ee(p->first,i->first,i->second);
				//if(visited.find(p->first) != visited.end()) continue;
				//if(visited.find(i->first) != visited.end()) continue;
				pQ.push(ee);
			}
		}
	}

	void primRandomMST(const Graph & G, Graph & T)
	{
		priority_queue<Edge> pQ;
		Graph::const_iterator p = G.begin();
		for(OutEdge::const_iterator i = p->second.begin(); i != p->second.end(); ++i){
			Edge e(p->first,i->first,i->second * (1 + rand01()));
			e.realwt = i->second;
			pQ.push(e);
		}
		set<uint64> visited;
		T.clear();
		while(!pQ.empty()){
			Edge e = pQ.top();
			pQ.pop();
			//if(visited.find(e.v1) != visited.end()) continue;
			if(visited.find(e.v2) != visited.end()) continue;
			visited.insert(e.v1);
			visited.insert(e.v2);
			T[e.v1][e.v2] = e.realwt;
			T[e.v2][e.v1] = e.realwt;
			p = G.find(e.v2);
			for(OutEdge::const_iterator i = p->second.begin(); i != p->second.end(); ++i){
				Edge ee(p->first,i->first,i->second * (1 + rand01()));
				e.realwt = i->second;
				pQ.push(ee);
			}
		}
	}

	void randomST2(const Graph & G, Graph & T)
	{
		priority_queue<Edge> pQ;
		Graph::const_iterator p = G.begin();
		for(OutEdge::const_iterator i = p->second.begin(); i != p->second.end(); ++i){
			Edge e(p->first,i->first,rand01() * i->second);
			e.realwt = i->second;
			pQ.push(e);
		}
		set<uint64> visited;
		T.clear();
		while(!pQ.empty()){
			Edge e = pQ.top();
			pQ.pop();
			if(visited.find(e.v2) != visited.end()) continue;
			visited.insert(e.v1);
			visited.insert(e.v2);
			T[e.v1][e.v2] = e.realwt;
			T[e.v2][e.v1] = e.realwt;
			p = G.find(e.v2);
			for(OutEdge::const_iterator i = p->second.begin(); i != p->second.end(); ++i){
				Edge ee(p->first,i->first,rand01() * i->second);
				ee.realwt = i->second;
				pQ.push(ee);
			}
		}
	}

	void randomST(const Graph & G, Graph & T){
          T.clear();
          set<uint64> V;
          vertexSet(G,V);
          uint64 v1 = *(V.begin());
          V.erase(v1);
          uint64 v2 = v1;
          while(!V.empty()){
               Graph::const_iterator oe = G.find(v1);
               vector<uint64> voe;
               for(OutEdge::const_iterator j = oe->second.begin(); j != oe->second.end(); ++j) voe.push_back(j->first);
               v2 = voe[rand() % voe.size()];
               if(V.find(v2) != V.end()){
                    T[v1][v2] = 1;
                    T[v2][v1] = 1;
                    V.erase(v2);
               }
               v1 = v2;
          }
     }

	void sampleByST(const Graph & G, Graph & S, int N){
        Graph T;
        for(int g = 0; g < N; ++g){
            randomST(G,T);
            for(Graph::iterator i = T.begin(); i != T.end(); ++i){
                for(OutEdge::iterator j = i->second.begin(); j != i->second.end(); ++j){
                    S[i->first][j->first] = j->second;
                }
            }
        }
	}

	void sampleByRandomMST(const Graph & G, Graph & S, int N){
		S.clear();
        primMST(G,S);
        for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
            double maxe = 0;
            for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j) maxe = max<double>(maxe,j->second);
            for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j){
                if(j->second > maxe * 0.7) S[i->first][j->first] = j->second;
            }
        }
        return;
        Graph T;
		for(int g = 0; g < N; ++g){
			if(g == 0) primMST(G,T);
			else primRandomMST(G,T);
			for(Graph::iterator i = T.begin(); i != T.end(); ++i){
				for(OutEdge::iterator j = i->second.begin(); j != i->second.end(); ++j){
					S[i->first][j->first] = j->second;
				}
			}
		}
	}


     void sampleRandomly(const Graph & G, Graph & S, double p){
          primMST(G,S);
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
               OutEdge oei = i->second;
               vector< pair<uint64,double> > oe(i->second.begin(), i->second.end());
               for(int k = 0; k < oe.size(); ++k) oe[k].second *= rand01();
               sort(oe.begin(), oe.end(), GreaterSecond<uint64,double>);
               int d = (int)(oe.size() * p);
               for(int j = 0; j < d; ++j){
                    double w = oei[oe[j].first];
                    S[i->first][oe[j].first] = w;
                    S[oe[j].first][i->first] = w;
               }
          }
     }

void meanShift(Graph & G){
	vector< Edge > es;
	for(Graph::iterator i = G.begin(); i != G.end(); ++i){
		for(OutEdge::iterator j = i->second.begin(); j != i->second.end(); ++j){
			if(i->first >= j->first) continue;
			Edge e(i->first, j->first, j->second);
		}
	}
	Graph Gt = G;
	for(int i = 0; i < es.size(); ++i){
		uint64 v1 = es[i].v1;
		uint64 v2 = es[i].v2;
		OutEdge o1 = G[v1];
		OutEdge o2 = G[v2];
		double w = es[i].wt;
		double nn = 1;
		for(OutEdge::iterator k = o1.begin(); k != o1.end(); ++k){
			w += k->second * exp(-3 * abs(k->second - w));
			nn += exp(-3 * abs(k->second - w));
		}
		w /= nn;
		Gt[v1][v2] = w;
		Gt[v2][v1] = w;
	}
	G = Gt;
}

	void layoutDim1(const Graph & G, map<uint64,double> & x)
	{
		vector<uint64> vs;
		for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
			x[i->first] = rand01() * 100;
			vs.push_back(i->first);
		}
		for(int g = 0; g < G.size() * 10; ++g){
			uint64 v = vs[rand() % vs.size()];
			Graph::const_iterator oe = G.find(v);
			double dx = 0;
			double xv = x[v];
			for(OutEdge::const_iterator i = oe->second.begin(); i != oe->second.end(); ++i){
				double fx = x[i->first] - xv;
				double dis = abs(fx) + 1;
				dis = dis / 8;
				dis *= dis;
				dx += fx * dis;
			}
			for(int i = 0; i < vs.size(); ++i){
				double fx = x[vs[i]] - xv;
				double dis = abs(fx) + 1;
				dis = dis / 8;
				dis *= dis;
				dx -= fx * dis;
			}
			dx = dx / abs(dx);
			xv += 3 * dx;
			if(xv < 10000 && xv > -10000) x[v] = xv;
		}
		double xmin = 1E20;
		double xmax = -1E20;
		for(map<uint64,double>::iterator i = x.begin(); i != x.end(); ++i){
			xmin = min<double>(xmin,i->second);
			xmax = max<double>(xmax,i->second);
		}
		for(map<uint64,double>::iterator i = x.begin(); i != x.end(); ++i){
			i->second = (i->second - xmin) / (xmax - xmin);
		}
	}

    /*
     void eigen(const Graph & G, vector< double > & egs){
          map<uint64,int> vindex;
          map<int,uint64> iv;
          int n = 0;
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i)
          {
               vindex[i->first] = n;
               iv[n] = i->first;
               ++n;
          }
          int N = G.size();

          LaGenMatDouble GA(N,N);
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
               int a = vindex[i->first];
               double d = 0;
               for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j) d += j->second;
               for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j){
                    int b = vindex[j->first];
                    GA(a,b) = j->second * (-1);
                    //GA(b,a) = j->second * (-1);
               }
               GA(a,a) = d;
          }

          LaVectorDouble eigreal(N),eigimg(N);
          LaGenMatDouble VR(N,N);
          LaEigSolve(GA,eigreal,eigimg,VR);
          egs.clear();
          for(int i = 0; i < N; ++i) egs.push_back(eigreal(i));
          sort(egs.begin(), egs.end());
     }

     double cosineMatrix(LaGenMatDouble & M, int ri, int rj){
          double co = 0;
          double si = 0;
          double sj = 0;
          for(int i = 0 ; i < M.cols(); ++i){
               if(i == ri) continue;
               if(i == rj) continue;
               si += M(ri,i) * M(ri,i);
               sj += M(rj,i) * M(rj,i);
               co += M(ri,i) * M(rj,i);
          }
          return co / sqrt(si * sj);
     }

     void relevance(const Graph & G, uint64 vq, map<uint64,double> & rel, int N){
          priority_queue<Edge> pQ;
          map<uint64,int> depth;
          depth[vq] = 0;
          Graph::const_iterator p = G.find(vq);
          for(OutEdge::const_iterator i = p->second.begin(); i != p->second.end(); ++i){
               Edge e(p->first,i->first,i->second);
               depth[i->first] = 1;
               pQ.push(e);
          }
          set<uint64> visited;
          while(!pQ.empty()){
               Edge e = pQ.top();
               pQ.pop();
               if(visited.find(e.v2) != visited.end()) continue;
               visited.insert(e.v1);
               visited.insert(e.v2);
               p = G.find(e.v2);
               int d = depth[e.v2];
               for(OutEdge::const_iterator i = p->second.begin(); i != p->second.end(); ++i){
                    Edge ee(p->first,i->first,i->second / d);
                    depth[i->first] = d + 1;
                    pQ.push(ee);
               }
               if(visited.size() > N) break;
          }

          Graph S;
          subgraph(G,visited,S);

          map<uint64,int> vindex;
          map<int,uint64> iv;
          int n = 0;
          for(Graph::iterator i = S.begin(); i != S.end(); ++i)
          {
               vindex[i->first] = n;
               iv[n] = i->first;
               ++n;
          }
          int K = S.size();
          cout << K << endl;
          LaGenMatDouble GA(K,K),A(K,K),B(K,K);
          GA = GA.zeros(K);
          for(Graph::iterator i = S.begin(); i != S.end(); ++i){
               int a = vindex[i->first];
               double d = 0;
               for(OutEdge::iterator j = i->second.begin(); j != i->second.end(); ++j) d += j->second;
               for(OutEdge::iterator j = i->second.begin(); j != i->second.end(); ++j){
                    int b = vindex[j->first];
                    GA(a,b) = j->second / d;
               }
          }
          B = GA; A = GA;
          for(int i = 0; i < 2; ++i){
               B = A * B;
               B *= 0.5;
               GA = GA + B;
          }
          int q = vindex[vq];
          for(int i = 0; i < K; ++i){
               if(i == q) continue;
               rel[iv[i]] = cosineMatrix(GA,q,i);
          }
     }
    */
     void vertexSet(const Graph & G, set<uint64> & V){
          for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
               V.insert(i->first);
               for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j)
                    V.insert(j->first);
          }
     }

     int vertices(const Graph & G){
          set<uint64> V;
          vertexSet(G,V);
          return V.size();
     }

     void betweennessCentrality(const Graph & G, map<uint64, double> & cb)
     {
          set<uint64> V;
          vertexSet(G,V);
          OutEdge d, delta, h;
          map<uint64, vector<uint64> > P;
          for(set<uint64>::iterator i = V.begin(); i != V.end(); ++i){
               cb[*i] = 0;
               d[*i] = 0;
               delta[*i] = 0;
               h[*i] = 0;
               P[*i].clear();
          }
          for(set<uint64>::iterator k = V.begin(); k != V.end(); ++k){
               //cout << *k << endl;
               stack<uint64> S;
               for(map<uint64, vector<uint64> >::iterator i = P.begin(); i != P.end(); ++i) i->second.clear();

               for(OutEdge::iterator i = d.begin(); i != d.end(); ++i) i->second = -1;
               d[*k] = 0;
               for(OutEdge::iterator i = delta.begin(); i != delta.end(); ++i) i->second = 0;
               delta[*k] = 1;
               deque<uint64> Q;
               Q.push_back(*k);
               while(!Q.empty()){
                    uint64 v = Q.front();
                    Q.pop_front();
                    S.push(v);
                    Graph::const_iterator oe = G.find(v);
                    if(oe == G.end()) continue;
                    for(OutEdge::const_iterator w = oe->second.begin(); w != oe->second.end(); ++w){
                         if(d[w->first] < 0){
                              Q.push_back(w->first);
                              d[w->first] = d[v] + 1;
                         }
                         if(d[w->first] = d[v] + 1){
                              delta[w->first] += delta[v];
                              P[w->first].push_back(v);
                         }
                    }
               }
               for(OutEdge::iterator i = h.begin(); i != h.end(); ++i) i->second = 0;
               while(!S.empty()){
                    uint64 w = S.top();
                    S.pop();
                    vector<uint64> & Pw = P[w];
                    double deltaw = delta[w];
                    for(vector<uint64>::iterator v = Pw.begin(); v != Pw.end(); ++v){
                         h[*v] += delta[*v] * (1 + h[w]) / deltaw;
                    }
                    if(w != *k){
                         cb[w] += h[w];
                         //cout << hw << endl;
                    }
               }
          }
     }



    void kmeanGraph(const Graph & G, int K, Partition & P){
        ofstream out("kmean.txt");
        vector< uint64 > vs;
        vector< Edge > es;
        for(Graph::const_iterator i = G.begin(); i != G.end(); ++i){
            P[i->first] = -1;
            vs.push_back(i->first);
            for(OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j){
                if(i->first >= j->first) continue;
                Edge e(i->first,j->first,j->second);
                es.push_back(e);
            }
        }
        sort(es.begin(), es.end());
        set<uint64> center;
        for(int i = 0; i < es.size(); ++i){
            center.insert(es[i].v1);
            center.insert(es[i].v2);
            if(center.size() > K * 8) break;
        }
        vector<uint64> vcenter(center.begin(), center.end());
        random_shuffle(vcenter.begin(), vcenter.end());
        for(int i = 0; i < K; ++i){
            P[vcenter[i]] = i;
        }
        for(int g = 0; g < 100; ++g){
            out << g << endl;
            bool stop = true;
            Partition P2 = P;
            for(int i = 0; i < vs.size(); ++i){
                uint64 v = vs[i];
                int pnow = P[v];
                Graph::const_iterator oe = G.find(v);
                vector< pair<uint64,double> > voe(oe->second.begin(), oe->second.end());
                sort(voe.begin(), voe.end(), GreaterSecond<uint64,double>);
                vector< double > sim(K,0);
                vector< double > nn(K,0.1);
                for(int j = 0; j < voe.size(); ++j){
                    int p = P[voe[j].first];
                    if(p < 0) continue;
                    sim[p] += voe[j].second * voe[j].second;
                    nn[p] += 1;
                }
                for(int j = 0; j < K; ++j){
                    sim[j] = sim[j] / nn[j];
                    sim[j] = sqrt(sim[j]);
                }
                double maxsim = 0;
                int p = 0;
                for(int j = 0; j < K; ++j){
                    if(maxsim < sim[j]){
                        maxsim = sim[j];
                        p = j;
                    }
                }
                if(p == pnow) continue;
                stop = false;
                double rd = rand01();
                if(rd < maxsim){
                    out << "p : " << v << "\t" << p << "\t" << maxsim << endl;
                    P2[v] = p;
                }
            }
            P = P2;
            if(stop && g > 20) break;
        }
        out.close();
    }
};
