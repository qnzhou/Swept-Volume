if (TARGET stf::stf)
    return()
endif()

include(CPM)
CPMAddPackage(
    NAME stf
    GITHUB_REPOSITORY adobe-research/space-time-functions
    GIT_TAG 1e16987a3dfa23b9c90728f6f942cfb25a6505a5
)

