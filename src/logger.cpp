#include <sweep/logger.h>

#include <spdlog/sinks/stdout_color_sinks.h>

namespace sweep {

spdlog::logger& logger()
{
    static auto default_logger = spdlog::stdout_color_mt("sweep");
    return *default_logger;
}

} // namespace sweep
