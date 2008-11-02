
def dumps(int imin, int imax, int rng=67108864):
    cdef char ilist[32]
    cdef int i = 0
    rng >>= 1
    cdef int mid = rng
    while rng > 1:
        rng >>=1
        if imin <= mid and imax <= mid:
            mid -= rng
            ilist[i] = '0' 
        elif imin > mid and imax > mid:
            mid += rng
            ilist[i] = '1' 
        else:
            break
        i += 1
    return ilist[:i]

def loads(hashed, int rng=67108864):
    cdef int i=0, imin = 0, imax = rng
    cdef char *chashed = hashed
    cdef char io = chashed[0]
    while io != '\0':
        i += 1
        rng >>= 1;
        if io == c'0': imax -= rng
        else: imin += rng
        io = chashed[i]
    return (imin, imax)
