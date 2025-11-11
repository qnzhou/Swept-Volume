if(TARGET arrangement::arrangement)
    return()
endif()

include(CPM)
CPMAddPackage(
    NAME arrangement
    GITHUB_REPOSITORY adobe-research/arrangement-benchmark
    GIT_TAG 5933c167f2277a10c1342e3e4c6de72bdfb94b4a
)
