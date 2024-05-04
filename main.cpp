#include <algorithm>
#include <bitset>
#include <cassert>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "fibonacci.hpp"

#define MAX_NUM_TERM 20

using cost = uint32_t;
struct Point {
   int32_t x, y, z;
};

auto
distance(const Point& p, const Point& q) -> cost {
   return abs(p.x - q.x) + abs(p.y - q.y) + abs(p.z - q.z);
}

class DijktraSteiner {
   using PointIdx = int32_t;
   using TermIdx = std::size_t;
   using Coord = int32_t;
   using CoordIdx = int32_t;

   struct PointRep {
      CoordIdx _x, _y, _z;

      PointRep() = default;
      PointRep(CoordIdx x, CoordIdx y, CoordIdx z) : _x(x), _y(y), _z(z){};

      bool
      operator==(const PointRep& other) const {
         return _x == other._x and _y == other._y and _z == other._z;
      }
      bool
      operator<(const PointRep& other) const {
         if (_x == other._x and _y == other._y)
            return _z < other._z;
         if (_x == other._x)
            return _y < other._y;
         return _x < other._y;
      }
      PointRep
      operator+(const PointRep& other) {
         return {_x + other._x, _y + other._y, _z + other._z};
      }
   };

   using TerminalSet = std::bitset<MAX_NUM_TERM>;
   using InstancePair = std::pair<PointRep, TerminalSet>;

   auto
   rep_to_idx(const PointRep& p) const -> PointIdx {
      return p._x + _x_coords.size() * (p._y + _y_coords.size() * p._z);
   }

   auto
   idx_to_rep(const PointIdx& p) const -> PointRep {
      return {p % (int)_x_coords.size(), (p / (int)_x_coords.size()) % (int)_y_coords.size(),
              p / ((int)_y_coords.size() * (int)_x_coords.size())};
   }

   auto
   rep_to_point(const PointRep& p) const -> Point {
      return {_x_coords[p._x], _y_coords[p._y], _z_coords[p._z]};
   }

   auto
   valid_point(const PointRep& p) const -> bool {
      return 0 <= p._x and p._x < (Coord)_x_coords.size() and 0 <= p._y and p._y < (Coord)_y_coords.size() and 0 <= p._z
             and p._z < (Coord)_z_coords.size();
   }

   auto
   lower_bound([[maybe_unused]] const InstancePair& vI) const -> cost {
      return 0;
   }

   auto
   length(const InstancePair& vI) const -> cost {
      if (_length[rep_to_idx(vI.first)].count(vI.second))
         return _length[rep_to_idx(vI.first)].at(vI.second);
      return std::numeric_limits<cost>::max() / 2;
   }

   auto
   update_length(const InstancePair& vI, cost new_cost) {
      if (_permanent_labels[rep_to_idx(vI.first)].count(vI.second)) {
         return;
      }
      // std::cerr << "Tentative subtree with root (" << _x_coords[vI.first._x] << ", " << _y_coords[vI.first._y] << ", "
      //           << _z_coords[vI.first._z] << ")\n";
      // std::cerr << "covering terminals " << vI.second << " and length " << new_cost << std::endl;
      if (_length[rep_to_idx(vI.first)].count(vI.second)) {
         _length[rep_to_idx(vI.first)][vI.second] = std::min(length(vI), new_cost);
         _pq.decrease_key(_fibheap_idx[rep_to_idx(vI.first)][vI.second]);
      } else {
         _length[rep_to_idx(vI.first)][vI.second] = new_cost;
         _fibheap_idx[rep_to_idx(vI.first)][vI.second] = _pq.emplace(vI);
      }
   }

   struct CompareInstancePairs {
      bool
      operator()(const InstancePair& a, const InstancePair& b) const {
         if (_DS.length(a) + _DS.lower_bound(a) == _DS.length(b) + _DS.lower_bound(b)) {
            if (a.first == b.first) {
               return a.second.to_ulong() < b.second.to_ulong();
            }
            return a.first < b.first;
         }
         return _DS.length(a) + _DS.lower_bound(a) < _DS.length(b) + _DS.lower_bound(b);
      }
      const DijktraSteiner& _DS;
   };

   // fields
   FibonacciHeap<InstancePair, CompareInstancePairs> _pq;
   std::vector<std::unordered_map<TerminalSet, cost>> _length;              // TODO: check if std::map is faster
   std::vector<std::unordered_map<TerminalSet, std::size_t>> _fibheap_idx;  // TODO: check if std::map is faster
   std::vector<std::unordered_set<TerminalSet>> _permanent_labels;
   PointIdx _num_terms;
   PointIdx _num_points;  // total number of points in the hanan grid
   std::vector<PointRep> _terminals;
   std::vector<Coord> _x_coords, _y_coords, _z_coords;

public:
   DijktraSteiner(std::vector<Point>&& terminals) : _pq(CompareInstancePairs{*this}), _num_terms(terminals.size()) {
      // extract hanan grid
      _x_coords.reserve(_num_terms);
      _y_coords.reserve(_num_terms);
      _z_coords.reserve(_num_terms);
      for (const auto& t : terminals) {
         _x_coords.push_back(t.x);
         _y_coords.push_back(t.y);
         _z_coords.push_back(t.z);
      }
      std::sort(_x_coords.begin(), _x_coords.end());
      std::sort(_y_coords.begin(), _y_coords.end());
      std::sort(_z_coords.begin(), _z_coords.end());
      _x_coords.erase(std::unique(_x_coords.begin(), _x_coords.end()), _x_coords.end());
      _y_coords.erase(std::unique(_y_coords.begin(), _y_coords.end()), _y_coords.end());
      _z_coords.erase(std::unique(_z_coords.begin(), _z_coords.end()), _z_coords.end());

      _terminals.reserve(_num_terms);
      for (const auto& t : terminals) {
         // get representation in hanan grid
         _terminals.emplace_back(std::find(_x_coords.begin(), _x_coords.end(), t.x) - _x_coords.begin(),
                                 std::find(_y_coords.begin(), _y_coords.end(), t.y) - _y_coords.begin(),
                                 std::find(_z_coords.begin(), _z_coords.end(), t.z) - _z_coords.begin());
         // TODO: replace with std::lower_bound
      }
      _num_points = _x_coords.size() * _y_coords.size() * _z_coords.size();
      _length.resize(_num_points);
      _fibheap_idx.resize(_num_points);
      _permanent_labels.resize(_num_points);
   }

   auto
   run() -> cost {
      TerminalSet remaining_bits(-1);
      for (PointIdx i = 1; i < _num_terms; ++i) {
         TerminalSet ts = 0;
         ts[i] = 1;
         _length[rep_to_idx(_terminals[i])][ts] = 0;
         _fibheap_idx[rep_to_idx(_terminals[i])][ts] = _pq.emplace({_terminals[i], ts});
         remaining_bits.flip(i);
      }
      while (not _pq.empty()) {
         PointRep v;
         TerminalSet I;
         std::tie(v, I) = _pq.delete_min();
         // std::cerr << "\nFound subtree with root (" << _x_coords[v._x] << ", " << _y_coords[v._y] << ", "
         //           << _z_coords[v._z] << ")\n";
         // std::cerr << "covering terminals " << I << " and length " << length({v, I}) << std::endl;
         if (v == _terminals[0] and (I | remaining_bits).all()) {
            return length({v, I});
         }
         std::vector<PointRep> dirs{
               {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1},
         };
         for (auto dir : dirs) {
            auto w = v + dir;
            if (not valid_point(w))
               continue;
            auto edgelength = distance(rep_to_point(v), rep_to_point(w));
            update_length({w, I}, length({v, I}) + edgelength);
         }
         for (auto J : _permanent_labels[rep_to_idx(v)]) {
            if ((J & (I | TerminalSet(1))).any())
               continue;
            update_length({v, I | J}, length({v, I}) + length({v, J}));
         }
         _permanent_labels[rep_to_idx(v)].insert(I);
      }
      return -1;
   }
};

int
main(int argc, char** argv) {
   assert(argc == 2);
   std::ifstream input(argv[1]);
   std::size_t n;
   input >> n;
   std::vector<Point> terminals(n);
   for (auto& p : terminals) {
      input >> p.x >> p.y >> p.z;
   }
   DijktraSteiner DS_inst(std::move(terminals));
   std::cout << DS_inst.run() << std::endl;
}
