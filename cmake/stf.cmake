if (TARGET stf::stf)
    return()
endif()

include(CPM)
CPMAddPackage(
    NAME stf
    GITHUB_REPOSITORY adobe-research/space-time-functions
    GIT_TAG 8a286acafc54edd5319eedb0e95d2ab17ad2512a
)

