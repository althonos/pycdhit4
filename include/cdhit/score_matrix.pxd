cdef extern from "cdhit-common.h" nogil:

    const size_t MAX_AA
    const size_t MAX_NA

    cppclass ScoreMatrix:
        int[MAX_AA][MAX_AA] matrix
        int gap
        int ext_gap

        ScoreMatrix()
        void init()
        void set_gap(int gap1, int ext_gap1)
        void set_matrix(const int* mat1)
        void set_to_na()
        void set_match(int score)
        void set_mismatch(int score)