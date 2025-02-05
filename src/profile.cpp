
#include "profile.hpp"

#include <sstream>
#include <unistd.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <signal.h>


namespace cuttlefish
{
    void profile(const std::function<void()>& f, const std::string& record)
    {
        // Fork a child.

        const auto pid = getpid();
        const auto forked = fork();
        if(forked == 0) // This is the child process.
        {
            // Launch profiler through the child.

            const auto fd = open("/dev/null", O_RDWR);
            dup2(fd, 1);
            dup2(fd, 2);

            // Execute `perf` from the child and attach it to the parent.
            // For precise-sampling to deal with perf-skid: https://www.brendangregg.com/perf.html
            exit(execl("/usr/bin/perf", "perf", "record", "-o", record.c_str(), "-e", "cpu-cycles:pp,cache-misses", "-p", std::to_string(pid).c_str(), nullptr));
        }
        else    // This is the parent process.
        {
            // Execute `f`.

            f();

            // Kill the profiler.

            kill(forked, SIGINT);
            waitpid(forked, nullptr, 0);
        }
    }
}
