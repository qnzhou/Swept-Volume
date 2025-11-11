if (TARGET stf::stf)
    return()
endif()

include(CPM)
CPMAddPackage(
    NAME stf
    GITHUB_REPOSITORY adobe-research/space-time-functions
    GIT_TAG 261207038fe50515a67b1231ea27bcdcc475abdf
)

