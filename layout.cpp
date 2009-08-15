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
          set<int> nodes;
     };

     class Node
     {
     public:
          Node(){}
          Node(string & n, int x, int y, int w, int h)
               :name(n),x(x),y(y),width(w),height(h),part(-1),rank(-1){}
          string name;
          string info;
          int x,y;
          int width,height;
          int part;
          int rank;
     };

     typedef CImg<double> Image;
     
     int width;
     int height;
     double posdis;
     double edgelen;

     vector<Node> nodes;
     vector< map<int,double> > graph;
     map< int, Group > grid;

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
               grid[k].nodes.insert(i);
               grid[k].x += nodes[i].x;
               grid[k].y += nodes[i].y;
          }
     }

     bool init()
     {
          double square = (double)(width*height);
          double vnum = (double)(nodes.size());
          posdis = 0.6 * sqrt(square/vnum);
          edgelen = posdis * 0.2;
          resetGrid();
          return true;
     }

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
                    for(set<int>::iterator i = p->second.nodes.begin(); i != p->second.nodes.end(); ++i) gd.insert(*i);
               }
          }
     }

     
     bool initialLayout()
     {
          srand(time(0));
          int i;
          for(i=0; i < nodes.size();i++){
               nodes[i].x = rand()%width;
               nodes[i].y = rand()%height;
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
          int x = nodes[i].x;
          int y = nodes[i].y;
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
     
     bool nextPositionGrid(int i, int ss)
     {
          int x = nodes[i].x;
          int y = nodes[i].y;
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
               if(p->second.nodes.empty()) continue;
               int xxj,yyj;
               gridPos(p->first, xxj, yyj);
               double xj = xxj; double yj = yyj;
               double dis = (x-xj)*(x-xj) + (y-yj)*(y-yj);
               dis = (1 + dis) / (posdis * posdis);
               double fx = (x - xj) * p->second.nodes.size();
               double fy = (y - yj) * p->second.nodes.size();
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
               pk1->second.nodes.erase(i);
               grid[k2].nodes.insert(i);
               if(pk1->second.nodes.empty())
                    grid.erase(pk1);
          }
          return false;
     }
     
     bool littleMove()
     {
          for(unsigned int i=0; i < nodes.size();i++)
          {
               nodes[i].x += (rand()%100-50);
               nodes[i].y += (rand()%100-50);
          }
          return true;
     }
     
     bool layout(int nn)
     {
          vector< pair<int,int> > vs;
          for(int i = 0; i < graph.size(); ++i) vs.push_back(make_pair<int,int>(i,graph[i].size()));
          sort(vs.begin(), vs.end(), network::GreaterSecond<int,int>);
          for(int k=0;k<nn;k++)
          {
               for(int j=0;j<vs.size();j++)
               {
                    int i = vs[j].first;
                    if(k < 50) nextPositionGrid(i,4);
                    else nextPositionGrid(i,1);
               }
          }
     }

     bool loadNetFile(const string & filename)
     {
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
          cout << "001" << endl;
          graph = vector< map<int,double> >(size);
          in>>buffer;
          int v1,v2;
          double value;
          while(!in.eof())
          {
               in>>v1>>v2>>value;
               graph[v1-1][v2-1] = value;
          }
          initialLayout();
          in.close();
          return true;
     }

     bool savePNGImg(const string & filename)
     {
          Image img(5000, 5000, 1, 3);
          img.save(filename.c_str());
          return true;
     }
};

int main(int argc, char ** argv){
     string src_file = argv[1];
     layout::loadNetFile(src_file);
     return 0;
}
