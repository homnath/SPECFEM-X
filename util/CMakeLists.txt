# make util
# REVISION:
#   HNG, Sep 06,2011;

# build exodus2sem
# new version
IF(BUILD_UTIL_EXODUSNEW2SEM)
  ADD_EXECUTABLE(exodusnew2sem
              trelis2mesh.c)
ENDIF(BUILD_UTIL_EXODUSNEW2SEM)

# old version
IF(BUILD_UTIL_EXODUSOLD2SEM)
  ADD_EXECUTABLE(exodusold2sem
              exodus2mesh.c)
ENDIF(BUILD_UTIL_EXODUSOLD2SEM)

# build write_sos
IF(BUILD_UTIL_WRITE_SOS)
  ADD_EXECUTABLE(write_sos
              write_sos.f90)
ENDIF(BUILD_UTIL_WRITE_SOS)

