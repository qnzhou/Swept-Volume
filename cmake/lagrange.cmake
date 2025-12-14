if(TARGET lagrange::core)
    return()
endif()

include(CPM)
CPMAddPackage(
    NAME lagrange-open
    GITHUB_REPOSITORY adobe/lagrange
    GIT_TAG v6.39.0
)
lagrange_include_modules(io)

