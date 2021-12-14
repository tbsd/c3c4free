#include <ctime>
#include <filesystem>
#include <iostream>
#include <thread>

#include "c3c4free.hpp"
#include "io.hpp"

void solve(const std::filesystem::path& inFile,
           const std::filesystem::path& outFile) {
  std::cout << "Now solveing: " << inFile << std::endl;
  graph::C3c4free task{graph::io::read(inFile)};
  task.approximate();
  graph::io::write(task, outFile);
}

// Чтобы скомпилить, нужно положить библиотеку boost в ../lib/boost/
// python3 imgs.py рисует картинки по файам задач/решений
// Так было выяснено, что распределение вершин совсем не равномерное
int main() {
  std::time_t seed = std::time(nullptr);
  srand(seed);
  std::cout << "Seed: " << seed << std::endl;
  //  std::thread t2(solve, "../Taxicab_2048.txt", "../Kurbatov_2048.txt");
  //  std::thread t1(solve, "../Taxicab_512.txt", "../Kurbatov_512.txt");
  //  std::thread t1(solve, "../Taxicab_4096.txt", "../Kurbatov_4096.txt");
  //    std::thread t2(solve, "../Taxicab_128.txt", "../Kurbatov_128.txt");
  //    solve("../Taxicab_4096.txt", "../Kurbatov_4096.txt");
  solve("../Taxicab_64.txt", "../Kurbatov_64.txt");
  //  solve("../Taxicab_128.txt", "../Kurbatov_128.txt");
  //  solve("../Taxicab_512.txt", "../Kurbatov_512.txt");
  //  solve("../Taxicab_2048.txt", "../Kurbatov_2048.txt");
  //  solve("../Taxicab_4096.txt", "../Kurbatov_4096.txt");
  //  t2.join();
  //  t1.join();
  return 0;
}
