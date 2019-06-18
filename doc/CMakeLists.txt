# build documentation
# REVISION:
#   HNG, Sep 06,2011;
IF(ISDOCUMENTATION)

  # copy all LaTeX source and required files 
  FILE(COPY ${SDIR}/doc DESTINATION ${CMAKE_BINARY_DIR})
  
  # compile documenent
  # pdflatex1
  ADD_CUSTOM_COMMAND( 
    OUTPUT    manual_SPECFEMX.aux
    DEPENDS   manual_SPECFEMX.tex
    COMMAND   ${PDFLATEX_COMPILER}
    ARGS      -interaction=batchmode manual_SPECFEMX    
    COMMENT   "pdflatex1"
  )

  # bibtex
  ADD_CUSTOM_COMMAND(
  OUTPUT    manual_SPECFEMX.bbl
  DEPENDS   manual_SPECFEMX.aux
  COMMAND   ${BIBTEX_COMPILER}
  ARGS      -terse manual_SPECFEMX
  COMMENT   "bibtex"
  )
  
  # pdflatex2
  ADD_CUSTOM_COMMAND( 
    OUTPUT    manual_SPECFEMX.pdf
    DEPENDS   manual_SPECFEMX.bbl
    COMMAND   ${PDFLATEX_COMPILER}
    ARGS      -interaction=batchmode manual_SPECFEMX
    
    COMMENT   "pdflatex2"
  )

  # pdflatex3
  ADD_CUSTOM_COMMAND( 
    OUTPUT    manual_SPECFEMX.log
    DEPENDS   manual_SPECFEMX.bbl
              manual_SPECFEMX.pdf
    COMMAND   ${PDFLATEX_COMPILER}
    ARGS      -interaction=batchmode manual_SPECFEMX
    
    COMMENT   "pdflatex3"
  )

  ADD_CUSTOM_TARGET(Documentation ALL echo
    DEPENDS   manual_SPECFEMX.log
  )
ENDIF(ISDOCUMENTATION)

