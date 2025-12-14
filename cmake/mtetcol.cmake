if (TARGET mtetcol::mtetcol)
    return()
endif()

include(CPM)
CPMAddPackage(
  NAME mtetcol
  GITHUB_REPOSITORY adobe-research/mtetcol
  GIT_TAG 8901d44292bd3bbc9c5fe36654c780f5040a06c9
)

