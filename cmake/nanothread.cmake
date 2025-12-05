if(TARGET nanothread::nanothread)
    return()
endif()

include(CPM)
CPMAddPackage(
  NAME nanothread
  GITHUB_REPOSITORY mitsuba-renderer/nanothread
  GIT_TAG 5d354917ea574488d133d242eda1a365a70b80e9
)

add_library(nanothread::nanothread ALIAS nanothread)
