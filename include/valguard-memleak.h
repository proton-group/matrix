#include <sstream>
#ifdef __unix__
#include <valgrind/memcheck.h>

std::stringstream valgrind_leaks()
{
    std::stringstream out;
    unsigned long leaked = 0;
    unsigned long dubious = 0;
    unsigned long rechable = 0;
    unsigned long suppresed = 0;
    VALGRIND_DO_CHANGED_LEAK_CHECK;
    VALGRIND_COUNT_LEAKS(leaked, dubious, rechable, suppresed);
    unsigned long sum = leaked + dubious + suppresed;
    if(sum != 0)
    {
        out << "Memory leak of " << sum << " byte(s) detected.\n";
    }
    return out;
}

#elif __APPLE__
#include <valgrind/memcheck.h>

std::stringstream valgrind_leaks()
{
    std::stringstream out;
    unsigned long leaked = 0;
    unsigned long dubious = 0;
    unsigned long rechable = 0;
    unsigned long suppresed = 0;
    VALGRIND_DO_CHANGED_LEAK_CHECK;
    VALGRIND_COUNT_LEAKS(leaked, dubious, rechable, suppresed);
    unsigned long sum = leaked + dubious + suppresed;
    if(sum != 0)
    {
        out << "Memory leak of " << sum << " byte(s) detected.\n";
    }
    return out;
}

#elif WIN32
#include "gtest-memleak.h"
std::stringstream valgrind_leaks()
{
    std::stringstream out;
    return out;
}
#endif