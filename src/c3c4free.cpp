#include "c3c4free.hpp"

#include <algorithm>
#include <exception>
#include <iostream>
#include <set>

#include <boost/graph/undirected_dfs.hpp>

#include "io.hpp"

namespace graph {

C3c4free::C3c4free(std::unique_ptr<RawVertices> vertices)
    : initSize(vertices->size()) {
  matrix = std::make_shared<GraphMatrix>(vertices->size());
  solutionGraph = std::make_shared<GraphMatrix>(vertices->size());
  best = Solution(solutionGraph);
  auto j = vertices->begin();
  auto [vert_begin, vert_end] = boost::vertices(*matrix);
  for (VertexIt i = vert_begin; i != vert_end; ++i, ++j) {
    getVM(*i) = *j;
  }
  j = vertices->begin();
  auto [vert_begin_s, vert_end_s] = boost::vertices(*solutionGraph);
  for (VertexIt i = vert_begin_s; i != vert_end_s; ++i, ++j) {
    getVS(*i) = *j;
  }
  for (VertexIt i = vert_begin; i != vert_end; ++i) {
    for (VertexIt n = i + 1; n != vert_end; ++n) {
      boost::add_edge(*i, *n, *matrix);
    }
  }
  auto [edges_begin, edges_end] = boost::edges(*matrix);
  for (auto i = edges_begin; i != edges_end; ++i) {
    int weight = manhattanDistance(getVM(i->m_source), getVM(i->m_target));
    (*matrix)[*i].weight = weight;
    //    getVS(i->m_source).weights.insert({weight, i->m_target});
    //    getVS(i->m_target).weights.insert({weight, i->m_source});
  }
  std::cout << "Vertices count: " << boost::num_vertices(*matrix)
            << ", edges count: " << boost::num_edges(*matrix) << std::endl;
  std::cout << "Initialization done" << std::endl;
}

C3c4free::VertexIt C3c4free::getVItS(int index) {
  auto [vert_begin, vert_end] = boost::vertices(*solutionGraph);
  return vert_begin + index;
}

VertexImpl& C3c4free::getVS(int index) {
  return (*solutionGraph)[*getVItS(index)];
}

C3c4free::VertexIt C3c4free::getVItM(int index) {
  auto [vert_begin, vert_end] = boost::vertices(*matrix);
  return vert_begin + index;
}

VertexImpl& C3c4free::getVM(int index) {
  return (*matrix)[*getVItM(index)];
}

int C3c4free::getStartVertex() {
  switch (boost::num_vertices(*matrix)) {
    case 64:
      return 1;
    case 128:
      return 42;
    case 512:
      return 0;
    case 2048:
      return 1080;
    case 4096:
      return 3015;
    default:
      return 0;
  }
}

int initVertNum = 0;

void C3c4free::approximate() {
  //  const int graphSz = boost::num_vertices(*matrix);
  curWeight = 0;
  unsigned long no_changes = 0;
  long prevWeight[] = {0, 0, 0, 0, 0};
  const int method = 0;
  const unsigned long no_changes_max = 10000;
  while (true) {
    //  for (int i = 0; i < 1; ++i) {
    //    try {

    //    } catch (...) {
    //      // Если что-то пошло не так, просто скипаем решение
    //      continue;
    //    }

    bool changes = false;
    //    changes = addRandomEdges2();
    //    switch (method) {
    //      case 0:
    changes = addRandomEdges();
    //        no_changes_max = 10000;
    //        break;
    //      case 1:
    //    changes = add5Cycle();
    //        no_changes_max = 10000;
    //        break;
    //      case 2:
    //        changes = false;
    //        addBigEdges();
    //        no_changes_max = 0;
    //        break;
    //    }
    //    changes = add5Cycle();
    if (changes)
      no_changes = 0;
    else
      ++no_changes;

    if (no_changes > no_changes_max) {
      auto [edges_begin, edges_end] = boost::edges(*solutionGraph);
      if (curWeight > prevWeight[method]) {
        best = Solution(solutionGraph);
        for (auto i = edges_begin; i != edges_end; ++i) {
          (*solutionGraph)[*i].weight =
              (*matrix)[boost::edge(i->m_source, i->m_target, *matrix).first]
                  .weight;
          best.addEdge(i);
        }
        prevWeight[method] = curWeight;

        io::write(*this, "../solutions/" + std::to_string(initSize) + "_" +
                             std::to_string(best.getWeight()) + "_" +
                             std::to_string(method) + ".txt");
        std::cout << "wheight: " << best.getWeight() << " size: " << best.size()
                  << " mehtod: " << method << std::endl;
      }
      const long edge_avg = curWeight / boost::num_edges(*solutionGraph);
      //      boost::num_edges(*solutionGraph); std::cout << "num_edges: " <<
      //      boost::num_edges(*solutionGraph)
      //                << std::endl;
      //          unsigned long j = 0;

      unsigned long solution_edge_num = 0;
      for (auto i = edges_begin; i != edges_end; ++i) {
        //        if (j % 3 == 0) {
        long weight =
            (*matrix)[boost::edge(i->m_source, i->m_target, *matrix).first]
                .weight;
        if (weight < edge_avg) {
          curWeight -= weight;
          boost::remove_edge(*i, *solutionGraph);
          //          removeEdge(solution_edge_num);
        }
        ++solution_edge_num;
        //        ++j;
      }
    }
    //      std::cout << "removed:   " << k << std::endl;
    //      std::cout << initSize << ":" << curWeight << std::endl;
    //    curWeight = 0;
    //      method = (method + 1) % 2;
    //    std::cout << "wheight: " << curWeight << std::endl;
  }
}
void C3c4free::removeEdge(unsigned long solution_edge_num) {
  auto [solBegin, solEnd] = boost::edges(*matrix);
  auto e = *std::next(solBegin, solution_edge_num);
  boost::remove_edge(e.m_source, e.m_target, *solutionGraph);
  std::set<unsigned long> now_available;
  auto [vBegin, vEnd] = boost::vertices(*matrix);
  for (auto i = vBegin; i != vEnd; ++i) {
    if (*i != e.m_source && !hasC3C4(e.m_source, *i)) {
      auto [adjBegin, adjEnd] = boost::adjacent_vertices(*i, *matrix);
      if (adjBegin != adjEnd)
        now_available.insert(*i);
      boost::add_edge(e.m_source, *i, *matrix);
    }
    if (*i != e.m_target && !hasC3C4(e.m_target, *i)) {
      auto [adjBegin, adjEnd] = boost::adjacent_vertices(*i, *matrix);
      if (adjBegin != adjEnd)
        now_available.insert(*i);
      boost::add_edge(e.m_target, *i, *matrix);
    }
  }
  for (const auto& v : now_available) {
    for (auto i = vBegin; i != vEnd; ++i) {
      if (*i != v && !hasC3C4(v, *i)) {
        boost::add_edge(v, *i, *matrix);
      }
    }
  }
  auto matrix_e = boost::edge(e.m_source, e.m_target, *matrix);
  curWeight -= (*matrix)[matrix_e.first].weight;
}

void C3c4free::addEdge(unsigned long solution_edge_num) {
  auto [solBegin, solEnd] = boost::edges(*matrix);
  auto e = *std::next(solBegin, solution_edge_num);
  boost::add_edge(e.m_source, e.m_target, *solutionGraph);
  curWeight += (*matrix)[e].weight;
  std::set<unsigned long> forbidden_v;
  auto [sBegin, sEnd] = boost::adjacent_vertices(e.m_source, *solutionGraph);
  for (auto i = sBegin; i != sEnd; ++i)
    forbidden_v.insert(*i);
  auto [tBegin, tEnd] = boost::adjacent_vertices(e.m_target, *solutionGraph);
  for (auto i = tBegin; i != tEnd; ++i)
    forbidden_v.insert(*i);

  std::set<unsigned long> forbidden_v_2;
  for (const auto& j : forbidden_v) {
    auto [s2Begin, s2End] = boost::adjacent_vertices(j, *solutionGraph);
    for (auto i = s2Begin; i != s2End; ++i) {
      boost::clear_vertex(*i, *matrix);
    }
    boost::clear_vertex(j, *matrix);
  }
}

bool C3c4free::addRandomEdges2() {
  unsigned long edges_count = boost::num_edges(*matrix);
  if (edges_count == 0)
    return false;
  unsigned long e_to_add = std::rand() % edges_count;

  addEdge(e_to_add);
  return true;
}

bool C3c4free::addRandomEdges() {
  unsigned long u = 0, v = 0;
  while (u == v) {
    u = std::rand() % initSize;
    v = std::rand() % initSize;
  }
  auto edge = boost::edge(v, u, *solutionGraph);
  if (!edge.second && !hasC3C4(u, v)) {
    boost::add_edge(u, v, *solutionGraph);
    long weight = (*matrix)[boost::edge(v, u, *matrix).first].weight;
    curWeight += weight;
    return true;
  }
  return false;
}

void C3c4free::addBigEdges() {
  std::map<long, Edge> edges;
  auto [edges_begin, edges_end] = boost::edges(*matrix);
  for (auto i = edges_begin; i != edges_end; ++i) {
    edges.insert({(*matrix)[*i].weight, *i});
  }
  for (auto& e : edges) {
    if (!hasC3C4(e.second.m_source, e.second.m_target)) {
      boost::add_edge(e.second.m_source, e.second.m_target, *solutionGraph);
      curWeight += e.first;
    }
  }
}
// bool addRandomEdge() {}

bool C3c4free::add5Cycle() {
  std::set<unsigned long> vertices;
  while (vertices.size() < 5) {
    vertices.insert(std::rand() % initSize);
  }
  std::set<Edge> initEdges;
  //  std::set<Edge> biggestEdges;
  long initWeight = 0;
  for (auto v = vertices.begin(); v != vertices.end(); ++v) {
    getVS(*v).inCurrentSolution = false;
    for (auto u = std::next(v); u != vertices.end(); ++u) {
      auto edge = boost::edge(*v, *u, *solutionGraph);
      if (edge.second) {
        initEdges.insert(edge.first);
        long weight = (*matrix)[boost::edge(*v, *u, *matrix).first].weight;
        initWeight += weight;
        //        std::cout << "w: " << weight;
      }
    }
  }
  std::set<Edge> newEdges;
  auto v = vertices.begin();
  long newWeight = 0;
  getVS(*vertices.begin()).inCurrentSolution = true;
  //  for (int i = 0; i < 5; ++i) {
  while (true) {
    auto bestV = vertices.begin();
    long maxWeight = 0;
    for (auto u = vertices.begin(); u != vertices.end(); ++u) {
      if (!getVS(*u).inCurrentSolution && u != v && !hasC3C4(*u, *v)) {
        //        auto edge = boost::edge(*v, *u, *solutionGraph);
        //        long weight = (*matrix)[edge.first].weight;
        long weight = (*matrix)[boost::edge(*v, *u, *matrix).first].weight;
        if (weight >= maxWeight) {
          maxWeight = weight;
          bestV = u;
          //          std::cout << "w2: " << weight;
        }
      }
    }
    if (maxWeight == 0)
      break;
    auto newEdge = boost::edge(*v, *bestV, *solutionGraph).first;
    newWeight += maxWeight;
    getVS(*bestV).inCurrentSolution = true;
    newEdges.insert(newEdge);
    v = bestV;
  }
  auto lastEdge = boost::edge(*vertices.begin(), *v, *solutionGraph).first;
  newWeight += (*matrix)[lastEdge].weight;
  newEdges.insert(lastEdge);

  if (newWeight > initWeight) {
    for (auto& e : initEdges)
      boost::remove_edge(e.m_source, e.m_target, *solutionGraph);
    for (auto& e : newEdges)
      boost::add_edge(e.m_source, e.m_target, *solutionGraph);
    curWeight = curWeight + newWeight - initWeight;
    return true;
  }

  return false;
}

bool C3c4free::hasC3C4(unsigned long source, unsigned long target) {
  std::set<unsigned long> source_adj_1, target_adj_1, intersection;
  auto [sBegin, sEnd] = boost::adjacent_vertices(source, *solutionGraph);
  for (auto i = sBegin; i != sEnd; ++i)
    source_adj_1.insert(*i);
  auto [tBegin, tEnd] = boost::adjacent_vertices(target, *solutionGraph);
  for (auto i = tBegin; i != tEnd; ++i)
    target_adj_1.insert(*i);
  std::set_intersection(source_adj_1.begin(), source_adj_1.end(),
                        target_adj_1.begin(), target_adj_1.end(),
                        std::inserter(intersection, intersection.begin()));
  if (!intersection.empty()) {
    //    std::cout << "3 cycle\n";
    return true;
  }

  std::set<unsigned long> source_adj_2, target_adj_2;
  for (const auto& j : source_adj_1) {
    auto [s2Begin, s2End] = boost::adjacent_vertices(j, *solutionGraph);
    for (auto i = s2Begin; i != s2End; ++i) {
      size_t prev_size = source_adj_2.size();
      source_adj_2.insert(*i);
      if (prev_size == source_adj_2.size())
        return true;
    }
  }
  for (const auto& j : target_adj_1) {
    auto [s2Begin, s2End] = boost::adjacent_vertices(j, *solutionGraph);
    for (auto i = s2Begin; i != s2End; ++i) {
      size_t prev_size = target_adj_2.size();
      target_adj_2.insert(*i);
      if (prev_size == target_adj_2.size())
        return true;
    }
  }
  std::set_intersection(source_adj_2.begin(), source_adj_2.end(),
                        target_adj_2.begin(), target_adj_2.end(),
                        std::inserter(intersection, intersection.begin()));
  if (!intersection.empty()) {
    //    std::cout << "4 cycle\n";
    return true;
  }
  return false;
}

const Solution<GraphMatrix> C3c4free::getSolution() const {
  return best;
}

int C3c4free::getProblemSize() const {
  return initSize;
}

bool C3c4free::hasCycle(unsigned long startV, int parentV) {
  getVM(startV).inCurrentSolution = true;
  for (auto& i : best.edges)
    if (i->m_source == startV || i->m_target == startV) {
      int newV = i->m_target == startV ? i->m_source : i->m_target;
      if (!getVM(newV).inCurrentSolution) {
        if (hasCycle(newV, startV))
          return true;
      } else if (newV != parentV)
        return true;
    }
  return false;
}

void C3c4free::clearMarks() {
  auto [vBegin, vEnd] = boost::vertices(*solutionGraph);
  for (VertexIt i = vBegin; i != vEnd; ++i) {
    getVS(*i).inCurrentSolution = false;
  }
}

bool C3c4free::check() {
  if (hasCycle((*best.edges.begin())->m_source, -1))
    return false;
  return true;
}

bool C3c4free::isReachable(unsigned long source,
                           unsigned long target,
                           unsigned long prev) {
  if (source == target)
    return true;
  auto [vBegin, vEnd] = boost::adjacent_vertices(source, *solutionGraph);
  for (auto i = vBegin; i != vEnd; ++i) {
    if (*i != prev && isReachable(*i, target, source))
      return true;
  }
  return false;
}

}  // namespace graph
