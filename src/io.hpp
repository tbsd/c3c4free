#pragma once

#include <filesystem>

#include "c3c4free.hpp"
#include "types.hpp"

namespace graph::io {

// Чтение из входного файла
std::unique_ptr<RawVertices> read(const std::filesystem::path& inPath);

// Вывод в файл
void write(const C3c4free& problem, const std::filesystem::path& outPath);

}  // namespace graph::io
