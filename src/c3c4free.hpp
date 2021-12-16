#pragma once

#include <filesystem>
#include <fstream>
#include <set>
#include <utility>

#include "solution.hpp"
#include "types.hpp"

namespace graph {
  // k-minimum spanning tree
  class C3c4free {
    using Vertex = boost::graph_traits<GraphMatrix>::vertex_descriptor;
    using Edge = GraphMatrix::edge_descriptor;
    using EdgeIt = GraphMatrix::out_edge_iterator;
    using VertexIt = GraphMatrix::vertex_iterator;

    public:
    C3c4free(std::unique_ptr<RawVertices> vertices);

    // Найти приближённое решение
    void approximate();

    const Solution<GraphMatrix> getSolution() const;

    int getProblemSize() const;

    // Соответствует ли решение условиям задачи
    bool check();

    private:
    //    void removeBigEdges(int startShift);

    /* Получить индекс вершины, оптимальной для начала поиска решения.
     * Вершины изначально выбраны руками на основе полученных картинок (потом
     * сделал перебор по всем вершинам).
     * Для общего случая имеет смысл искать
     * скопления точек и в них выбирать центральную или с наименьшим ребом
     */
    int getStartVertex();

    // Вершина по индексу
    VertexImpl& getVS(int index);

    // Итератор вершины по индексу
    VertexIt getVItS(int index);

    // Вершина по индексу
    VertexImpl& getVM(int index);

    // Итератор вершины по индексу
    VertexIt getVItM(int index);

    bool hasCycle(unsigned long startV, int parentV);

    // Проверяет, дас ли добавление ребра {source, target} 3- или 4- цикл
    bool hasC3C4(unsigned long source, unsigned long target);

    // Добавлять рёбра в порядке убывания веса, если это возможно
    void addBigEdges();

    bool addRandomEdges();

    bool addRandomEdges2();
    bool addEdge(unsigned long e);
    void removeEdge(Edge e);

    // Создаём большой цикл, ходим по нему, соединяя вершины на расстаянии 5
    void cycleWithCycles();

    void clearMarks();

    void recalculateAvalibleEdges();

    // Добавить случайный цикл на 5 вершинах. Возвращает true, если были
    // изменены рёбра, иначае false
    bool add5Cycle();

    bool isReachable(unsigned long source,
                     unsigned long target,
                     unsigned long prev);

    std::set<unsigned long> available_e;
    const int initSize;
    long curWeight = 0;
    std::shared_ptr<GraphMatrix> matrix;
    std::shared_ptr<GraphMatrix> availabileEGraph;
    std::shared_ptr<GraphMatrix> solutionGraph;
    long bestWeight = 0;
    Solution<GraphMatrix> best;
  };

  }  // namespace graph
