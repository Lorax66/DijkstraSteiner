#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <set>
#include <vector>

template <typename T, class Compare>
class FibonacciHeap {
   using VertexIdx = int;
   struct Vertex {
      bool phi;
      T element;
      VertexIdx parent;
      std::size_t child_idx_at_parent;
      std::vector<VertexIdx> children;

      Vertex(T t, VertexIdx idx) {
         phi = false;
         element = t;
         parent = idx;
      }
   };

public:
   FibonacciHeap<T, Compare>(const Compare&& _comparator) : _num_removed(0), _comparator(_comparator){};

   // returns the index needed to use decrease_key
   auto
   emplace(T t) -> VertexIdx {
      auto idx = _vertices.size();
      _vertices.emplace_back(t, idx);
      plant(idx);
      return idx;
   }

   auto
   delete_min() -> T {
      auto min_idx = *std::min_element(_roots.begin(), _roots.end(), [this](auto v, auto w) -> bool {
         if (v == -1)
            return false;
         if (w == -1)
            return true;
         return _comparator(_vertices[v].element, _vertices[w].element);
      });
      assert(min_idx != -1);
      for (auto& r : _roots) {
         if (r == min_idx) {
            r = -1;
         }
      }
      for (auto c : _vertices[min_idx].children) {
         _vertices[c].parent = c;
         plant(c);
      }
      remove_parent_edge(min_idx);
      _num_removed++;
      return _vertices[min_idx].element;
   }

   void
   decrease_key(VertexIdx v) {
      std::vector<VertexIdx> p;
      while (_vertices[v].parent != v and _vertices[v].phi) {
         p.push_back(v);
         v = _vertices[v].parent;
         remove_parent_edge(p.back());
         _vertices[v].phi = not _vertices[v].phi;
      }
      for (auto u : p) {
         if (_vertices[u].children.empty()) {
            plant(u);
         }
      }
   }

   auto
   empty() -> bool {
      return _vertices.size() == _num_removed;
   }

private:
   void
   plant(VertexIdx v) {
      assert(v != -1);
      assert(_vertices[v].parent == v);
      _roots.resize(std::max(_roots.size(), _vertices[v].children.size() + 1), -1);
      auto r = _roots[_vertices[v].children.size()];
      if (r != -1 and _vertices[r].parent == r and r != v
          and _vertices[r].children.size() == _vertices[v].children.size()) {
         if (_comparator(_vertices[r].element, _vertices[v].element)) {
            add_edge(v, r);
            plant(r);
         } else {
            assert(_comparator(_vertices[v].element, _vertices[r].element));
            // _roots[_vertices[v].children.size()] = -1;
            add_edge(r, v);
            plant(v);
         }
      } else {
         _roots[_vertices[v].children.size()] = v;
      }
   }

   void
   add_edge(VertexIdx v, VertexIdx p) {
      assert(_vertices[v].parent == v);
      _vertices[v].parent = p;
      _vertices[v].child_idx_at_parent = _vertices[p].children.size();
      _vertices[p].children.push_back(v);
   }

   void
   remove_parent_edge(VertexIdx v) {
      if (_vertices[v].parent == v)
         return;
      _vertices[_vertices[_vertices[v].parent].children.back()].child_idx_at_parent = _vertices[v].child_idx_at_parent;
      _vertices[_vertices[v].parent].children[_vertices[v].child_idx_at_parent]
            = _vertices[_vertices[v].parent].children.back();
      _vertices[_vertices[v].parent].children.pop_back();
      for (auto c : _vertices[_vertices[v].parent].children) {
         assert(_vertices[_vertices[c].parent].children[_vertices[c].child_idx_at_parent] == c);
      }
      assert(_vertices[_vertices[v].parent].children.size()
             == std::set<VertexIdx>(_vertices[_vertices[v].parent].children.begin(),
                                    _vertices[_vertices[v].parent].children.end())
                      .size());
      _vertices[v].parent = v;
   }

   uint32_t _num_removed;
   std::vector<VertexIdx> _roots;
   std::vector<Vertex> _vertices;
   const Compare _comparator;
};
