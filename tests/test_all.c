#include <stdio.h>
#include <stdint.h>
#include <tests/rollo-i-128-tests/all_rollo-i-128_tests.h>
#include <time.h>


uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32u) | lo;
}

void printresults(const char* measured_function, clock_t total_time, uint64_t total_iterations, uint64_t total_cycles){
    double total_time_in_seconds = (double) total_time /
                                   (double) CLOCKS_PER_SEC;
    double operations_per_second = (double) total_iterations * (double) CLOCKS_PER_SEC/ (double) total_time;
    double cycles_per_op = (double) total_cycles/ (double) total_iterations;

    printf("\r%60s:  %12.2f cycles/op, (%3.2f s, %4.3e ops, %4.3e cycles, %3.2f op/s, %4.3f GHz)\n", measured_function,
           cycles_per_op, total_time_in_seconds, (double) total_iterations, (double) total_cycles, operations_per_second, cycles_per_op*operations_per_second/1000000000.);
}

int main(int argc, char **argv) {
    putchar('\n');

    int result = 0;
    int total_failed = 0;

    result = rollo_i_128_unit_tests();
    total_failed += result;

    printf("\nTotal failed tests: %d\n", total_failed);

    rollo_i_128_performance_tests();

    return 0;
}
