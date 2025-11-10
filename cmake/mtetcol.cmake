if (TARGET mtetcol::mtetcol)
    return()
endif()

include(CPM)
CPMAddPackage(
  NAME mtetcol
  GITHUB_REPOSITORY adobe-research/mtetcol
  GIT_TAG 68f73c6630b2b1a2b724c88cc375b8c6350cb961
)

