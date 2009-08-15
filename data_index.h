#ifndef _DATA_INDEX_H_
#define _DATA_INDEX_H_

#include <vector>
#include <map>
#include <set>

namespace pattern{
    using namespace std;
    template <typename _Point> struct Cell2D{
        double xmin,ymin,xmax,ymax;
        int N;
        double xg,yg;
        map< int, set<int> > c;
        vector< _Point> points;
    };
    template <typename _Point> void buildCell2D(const vector< _Point > & ps, Cell2D< _Point > & cell){
        cell.points = ps;
        cell.xmin = ps[0].first;
        cell.xmax = cell.xmin;
        cell.ymin = ps[0].second;
        cell.ymax = cell.ymin;
        for(int i = 0; i < ps.size(); ++i){
            cell.xmin = min<double>(cell.xmin,ps[i].first);
            cell.xmax = max<double>(cell.xmax,ps[i].first);
            cell.ymin = min<double>(cell.ymin,ps[i].second);
            cell.ymax = max<double>(cell.ymax,ps[i].second);
        }
        cell.xg = (cell.xmax - cell.xmin)/((double)(cell.N));
        cell.yg = (cell.ymax - cell.ymin)/((double)(cell.N));
        for(int i = 0; i < ps.size(); ++i){
            int xn = (int)((ps[i].first - cell.xmin) / cell.xg);
            int yn = (int)((ps[i].second - cell.ymin) / cell.yg);
            int id = xn * cell.N + yn;
            cell.c[id].insert(i);
        }
    }
    template <typename _Point> void searchCell2D(const Cell2D< _Point > & cell, _Point & p, int K, set< int > & rp){
        int xn = (int)((p.first - cell.xmin) / cell.xg);
        int yn = (int)((p.second - cell.ymin) / cell.yg);
        for(int a = -1 * K; a <= K; ++a){
            for(int b = -1 * K; b <= K; ++b){
                int x = xn + a;
                int y = yn + b;
                int id = x * cell.N + y;
                map<int, set<int> >::const_iterator pp = cell.c.find(id);
                if(pp == cell.c.end()) continue;
                set< int > tp = pp->second;
                for(set<int>::iterator i = tp.begin(); i != tp.end(); ++i) rp.insert(*i);
            }
        }
    }

};

#endif
