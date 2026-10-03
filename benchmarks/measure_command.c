/* Fresh, small native launcher: child RSS must not inherit a Python validator's
 * large address space. Timings cover exec through completed command exit. */
#define _GNU_SOURCE 1
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/resource.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

int main(int argc, char **argv) {
    if (argc < 3) return 2;
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);
    pid_t child = fork();
    if (child < 0) return 3;
    if (child == 0) {
        execv(argv[2], argv + 2);
        perror("execv");
        _exit(127);
    }
    struct rusage usage;
    int status = 0;
    while (wait4(child, &status, 0, &usage) < 0) {
        if (errno != EINTR) return 4;
    }
    clock_gettime(CLOCK_MONOTONIC, &end);
    int code = WIFEXITED(status) ? WEXITSTATUS(status) : 128 + WTERMSIG(status);
    double seconds = end.tv_sec - start.tv_sec + (end.tv_nsec - start.tv_nsec) * 1e-9;
    double cpu = usage.ru_utime.tv_sec + usage.ru_utime.tv_usec * 1e-6
               + usage.ru_stime.tv_sec + usage.ru_stime.tv_usec * 1e-6;
    double rss = usage.ru_maxrss;
#ifdef __APPLE__
    rss /= 1024.;
#endif
    FILE *output = fopen(argv[1], "w");
    if (!output) return 5;
    fprintf(output, "{\"seconds\":%.9f,\"cpu_seconds\":%.9f,\"peak_rss_kib\":%.0f,\"exit_code\":%d}\n",
            seconds, cpu, rss, code);
    if (fclose(output) != 0) return 6;
    return code;
}
