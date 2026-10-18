// Golden data for fdsp_blackman (0.6.0): fldigi 4.2.13's blackman() of src/include/misc.h
// as its reference build compiles it (Apple clang -O2, FP contraction on, arm64): the FSQ and
// IFKP receivers fill their 4096-point window with blackman(1.0 * i / 4096) (fsq.cxx:206-207,
// ifkp.cxx:162-163). Generate on that platform (make golden).
#include <cstdint>
#include <cstdio>
#include <vector>
#include "misc.h"

int main(int argc, char **argv)
{
    const char *dir = argc > 1 ? argv[1] : ".";
    std::vector<double> out;
    const int sizes[3] = {4096, 1024, 257};
    for (int k = 0; k < 3; k++)
        for (int i = 0; i < sizes[k]; i++)
            out.push_back(blackman(1.0 * i / sizes[k]));
    char path[1024];
    snprintf(path, sizeof path, "%s/blackman.f64", dir);
    FILE *f = fopen(path, "wb");
    if (!f) return 1;
    fwrite(out.data(), sizeof(double), out.size(), f);
    fclose(f);
    return 0;
}
