#include "network.h"
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <algorithm>
#include <cmath>
#include "CImg.h"
using namespace cimg_library;
using namespace std;

namespace layout{
     struct Group
     {
          int x;
          int y;
          int n;
          set<int> ns;
     };

     class Node
     {
     public:
          Node(){}
          Node(string & n, double xx, double yy)
               :name(n),x(xx),y(yy),part(-1),rank(-1){}
          string name;
          string info;
          double x,y;
          int part;
          int rank;
     };

     typedef unsigned long uint64;
     typedef CImg<double> Image;
     
     int width;
     int height;
     double posdis;
     double edgelen;

     vector<Node> nodes;
     vector< map<int,double> > graph;
     map< int, Group > grid;

     bool toMapGraph(network::Graph & G)
     {
          int n = graph.size();
          for(int i = 0; i < n; ++i){
               for(map<int,double>::iterator j = graph[i].begin(); j != graph[i].end(); ++j){
                    uint64 ii = i;
                    uint64 jj = j->first;
                    G[ii][jj] = j->second;
               }
          }
          return true;
     }

     bool toVectorGraph(const network::Graph & G)
     {
          int n = G.size();
          graph.clear();
          graph = vector< map<int,double> >(n);
          for(network::Graph::const_iterator i = G.begin(); i != G.end(); ++i){
               for(network::OutEdge::const_iterator j = i->second.begin(); j != i->second.end(); ++j){
                    int v1 = (int)(i->first);
                    int v2 = (int)(j->first);
                    graph[v1][v2] = j->second;
               }
          }
          return true;
     }

     bool MST()
     {
          network::Graph G,T;
          toMapGraph(G);
          graph.clear();
          network::primMST(G,T);
          toVectorGraph(T);
          return true;
     }
     /*
     int gridIndex(int x, int y)
     {
          if(x < -10000 || y < -10000) return 0;
          int xx = (int)((double)(x + 10000) / 48);
          int yy = (int)((double)(y + 10000) / 48);
          return xx * 20000 + yy;
     }
     
     void resetGrid()
     {
          grid.clear();
          for(int i=0; i < nodes.size();i++)
          {
               int k = gridIndex(nodes[i].x, nodes[i].y);
               grid[k].ns.insert(i);
               grid[k].x += nodes[i].x;
               grid[k].y += nodes[i].y;
          }
     }
     */
     bool init()
     {
          double square = (double)(width*height);
          double vnum = (double)(nodes.size());
          posdis = 0.6 * sqrt(square/vnum);
          edgelen = posdis * 0.2;
          //resetGrid();
          return true;
     }
     /*
     void gridPos(int i, int & x, int & y)
     {
          int a = i % 20000;
          int b = (int)((double)(i) / 20000);
          y = (int)(a * 48 - 10000);
          x = (int)(b * 48 - 10000);
     }
     
     void gridNodes(int x, int y, set<int> & gd)
     {
          gd.clear();
          for(int i = -1; i <= 1; ++i){
               for(int j = -1; j <= 1; ++j){
                    double xx = (double)(x) + posdis * (double)(i);
                    double yy = (double)(y) + posdis * (double)(j);
                    int k = gridIndex((int)(xx), (int)(yy));
                    map<int, Group >::iterator p = grid.find(k);
                    if(p == grid.end()) continue;
                    for(set<int>::iterator i = p->second.ns.begin(); i != p->second.ns.end(); ++i) gd.insert(*i);
               }
          }
     }
     */
     
     bool initialLayout()
     {
          srand(time(0));
          int i;
          for(i=0; i < nodes.size();i++){
               nodes[i].x = (double)(rand() % width) / 1000;
               nodes[i].y = (double)(rand() % height) / 1000;
          }
          init();
          return true;
     }

     double distance2(int a, int b)
     {
          double dx = (double)(nodes[a].x - nodes[b].x);
          double dy = (double)(nodes[a].y - nodes[b].y);
          return dx*dx+dy*dy;
     }
     
     bool nextPosition(int i, int ss)
     {
          double x = nodes[i].x;
          double y = nodes[i].y;
          double dx,dy;
          dx = 0; dy = 0;
          for(map<int,double>::iterator p = graph[i].begin(); p != graph[i].end(); ++p)
          {
               int j = p->first;
               if(i == j) continue;
               double xj = nodes[j].x;
               double yj = nodes[j].y;
               double dis = (1 + distance2(i,j)) / (posdis*posdis);
               double fx = (x - xj);
               double fy = (y - yj);
               dx -= fx * dis;
               dy -= fy * dis;
          }
          for(int j = 0; j < graph.size(); ++j)
          {
               if(i == j) continue;
               double xj = nodes[j].x;
               double yj = nodes[j].y;
               double dis = (1 + distance2(i,j)) / (posdis*posdis);
               if(dis > 100) continue;
               double fx = (x - xj);
               double fy = (y - yj);
               dx += fx / dis;
               dy += fy / dis;
          }
          double e = sqrt(dx*dx + dy*dy);
          dx = dx / e;
          dy = dy / e;
          e = min<double>(2,sqrt(e));
          x += (int)(ss * dx * e);
          y += (int)(ss * dy * e);
          if(x < -10000 || y < -10000 || x > 10000 || y > 10000) return false;
          nodes[i].x = x;
          nodes[i].y = y;
          return false;
     }
     /*
     bool nextPositionGrid(int i, int ss)
     {
          double x = nodes[i].x;
          double y = nodes[i].y;
          double dx,dy;
          dx = 0; dy = 0;
          for(map<int,double>::iterator p = graph[i].begin(); p != graph[i].end(); ++p)
          {
               int j = p->first;
               if(i == j) continue;
               double xj = nodes[j].x;
               double yj = nodes[j].y;
               double dis = (1 + distance2(i,j)) / (posdis*posdis);
               double fx = (x - xj);
               double fy = (y - yj);
               dx -= fx * dis;
               dy -= fy * dis;
               dx += fx / dis;
               dy += fy / dis;
          }
          int k1 = gridIndex(nodes[i].x, nodes[i].y);
          
          for(map<int, Group >::iterator p = grid.begin(); p != grid.end(); ++p)
          {
               if(p->second.ns.empty()) continue;
               int xxj,yyj;
               gridPos(p->first, xxj, yyj);
               double xj = xxj; double yj = yyj;
               double dis = (x-xj)*(x-xj) + (y-yj)*(y-yj);
               dis = (1 + dis) / (posdis * posdis);
               double fx = (x - xj) * p->second.ns.size();
               double fy = (y - yj) * p->second.ns.size();
               dx += fx / dis;
               dy += fy / dis;
          }
          double e = sqrt(dx*dx + dy*dy);
          dx = dx / e;
          dy = dy / e;
          x += (int)((double)(ss) * dx);
          y += (int)((double)(ss) * dy);
          if(x < -10000 || y < -10000 || x > 10000 || y > 10000) return false;
          nodes[i].x = x; nodes[i].y = y;
          int k2 = gridIndex(nodes[i].x, nodes[i].y);
          if(k1 != k2){
               map<int,Group>::iterator pk1 = grid.find(k1);
               pk1->second.ns.erase(i);
               grid[k2].ns.insert(i);
               if(pk1->second.ns.empty())
                    grid.erase(pk1);
          }
          return false;
     }
     */
     bool littleMove()
     {
          int i = rand() % nodes.size();
          nodes[i].x += (double)(rand()%100-50);
          nodes[i].y += (double)(rand()%100-50);
          return true;
     }
     
     bool layout(int nn)
     {
          vector< pair<int,int> > vs;
          for(int i = 0; i < graph.size(); ++i) vs.push_back(make_pair<int,int>(i,graph[i].size()));
          sort(vs.begin(), vs.end(), network::GreaterSecond<int,int>);
          for(int k = 0; k < nn; k++)
          {
               if(k % 100 == 0) cout << "\n" << k << "\t";
               if(k % 10 == 0) cout << ".";
               cout.flush();
               for(int j=0;j<vs.size();j++)
               {
                    int i = vs[j].first;
                    nextPosition(i,4);
                    //if(rand() % 5 == 0 && k < nn - 30) littleMove();
               }
          }
     }

     bool loadNetFile(const string & filename)
     {
          width = 8000;
          height = 8000;
          graph.clear();
          nodes.clear();
          ifstream in(filename.c_str());
          if(in.fail()) return false;
          string buffer;
          int size,index;
          in>>buffer>>size;
          nodes = vector<Node>(size);
          string line;
          getline(in,line);
          for(int i = 0; i < size; i++){
               getline(in,line);
               istringstream iss(line);
               iss >> index >> buffer;
               nodes[index-1].x = 0;
               nodes[index-1].y = 0;
               nodes[index-1].name = buffer;
               nodes[index-1].info = line;
               nodes[index-1].rank = -1;
          }
          graph = vector< map<int,double> >(size);
          in>>buffer;
          int v1,v2;
          double value;
          while(!in.eof())
          {
               in>>v1>>v2>>value;
               graph[v1-1][v2-1] = value;
               graph[v2-1][v1-1] = value;
          }
          initialLayout();
          in.close();

          MST();
          cout << "load file [" << filename << "] ok!" << endl;
          return true;
     }

     bool savePNGImg(const string & filename)
     {
          double xmax = -1000000;
          double xmin = 1000000;
          double ymax = -1000000;
          double ymin = 1000000;
          for(int i = 0; i < nodes.size(); ++i){
               xmax = max<double>(xmax, nodes[i].x);
               xmin = min<double>(xmin, nodes[i].x);
               ymax = max<double>(ymax, nodes[i].y);
               ymin = min<double>(ymin, nodes[i].y);
          }
          for(int i = 0; i < nodes.size(); ++i){
               nodes[i].x = 100 + (double)(width - 200) * (nodes[i].x - xmin) / max<double>(xmax - xmin, ymax - ymin);
               nodes[i].y = 100 + (double)(height - 200) * (nodes[i].y - ymin) / max<double>(xmax - xmin, ymax - ymin);
          }
          Image img(width, height, 1, 3);
          double red_color[3] = {255,0,0};
          double green_color[3] = {0,255,0};
          double black_color[3] = {0,0,0};
          for(int i = 0; i < nodes.size(); ++i){
               img.draw_circle((int)(nodes[i].x), (int)(nodes[i].y), 2, red_color);
          }

          for(int i = 0; i < graph.size(); ++i){
               for(map<int,double>::iterator j = graph[i].begin(); j != graph[i].end(); ++j){
                    int x1 = (int)(nodes[i].x);
                    int y1 = (int)(nodes[i].y);
                    int x2 = (int)(nodes[j->first].x);
                    int y2 = (int)(nodes[j->first].y);
                    double color[3] = {200 * sqrt(j->second) + 55,200 * sqrt(j->second) + 55,0};
                    img.draw_line(x1, y1, x2, y2, color);
               }
          }

          for(int i = 0; i < nodes.size(); ++i){
               img.draw_text((int)(nodes[i].x + 5), (int)(nodes[i].y), nodes[i].name.c_str(), green_color, black_color);
          }
          img.save(filename.c_str());
          return true;
     }
};

int main(int argc, char ** argv){
     string src_file = argv[1];
     string out_file = "out.png";
     layout::loadNetFile(src_file);
     layout::layout(1000);
     layout::savePNGImg(out_file);
     cout << endl;
     return 0;
}
