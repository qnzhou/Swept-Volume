if (TARGET yaml-cpp::yaml-cpp)
    return()
endif()

message(STATUS "Third-party (external): creating target 'yaml-cpp::yaml-cpp'")

include(CPM)
CPMAddPackage(
    NAME yaml-cpp
    GITHUB_REPOSITORY jbeder/yaml-cpp
    GIT_TAG yaml-cpp-0.7.0
)
