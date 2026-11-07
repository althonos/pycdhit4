

cdef extern from "cdhit-common.h" nogil:

    const char* CDHIT_VERSION
    
    const size_t MAX_SEQ
    const size_t MAX_AA
    const size_t MAX_NA
    const size_t MAX_UAA
    const size_t MAX_DIAG
    const size_t MAX_GAP
    const size_t MAX_DES
    const size_t MAX_LINE_SIZE
    const size_t MAX_FILE_NAME
    const size_t MAX_SEG
    const size_t MAX_BIN_SWAP
    const size_t MAX_TABLE_SIZE
    const size_t CLOCK_TICKS

    const char[MAX_AA] aa
    const char[26] aa2idx
    const char[26] na2idx
    
    enum:
        FAILED_FUNC
        OK_FUNC

    enum:
        IS_REP
        IS_REDUNDANT
        IS_PROCESSED
        IS_MINUS_STRAND


    