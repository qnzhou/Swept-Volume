if (TARGET stf::stf)
    return()
endif()

include(CPM)
CPMAddPackage(
    NAME stf
    GITHUB_REPOSITORY adobe-research/space-time-functions
    GIT_TAG 80da61c0a693b7ef96a0150a8165cf2e6d546643
)
