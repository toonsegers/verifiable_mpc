import functools
   
    
@functools.lru_cache(maxsize=None)
def _recombination_vectors(field, xs, xr):
    """Compute and store recombination vectors.

    Recombination vectors depend on the field, the x-coordinates xs
    of the shares and the x-coordinates xr of the recombination points.
    """
    modulus = field.modulus
    xs = [x % modulus for x in xs]  # also for conversion from
    xr = [x % modulus for x in xr]  # int to type(modulus)
    d = [None] * len(xs)
    for i, x_i in enumerate(xs):
        q = 1
        for j, x_j in enumerate(xs):
            if i != j:
                q *= x_i - x_j
                q %= modulus
        d[i] = q
    matrix = [None] * len(xr)
    for r, x_r in enumerate(xr):
        matrix[r] = [None] * len(xs)
        p = 1
        for j, x_j in enumerate(xs):
            p *= x_r - x_j
            p %= modulus
        p = field(p)
        for i, x_i in enumerate(xs):
            matrix[r][i] = (p / field((x_r - x_i) * d[i])).value
    return matrix


def recombine(field, points, x_rs=0):  ##### ONLY for shares that are single numbers
    """Recombine shares given by points into secrets.

    Recombination is done for x-coordinates x_rs.
    """
    xs, shares = list(zip(*points))
    if not isinstance(x_rs, list):
        x_rs = (x_rs,)
    m = len(shares)
    width = len(x_rs)
    T_is_field = isinstance(shares[0], field)  # all elts assumed of same type
    vector = _recombination_vectors(field, xs, tuple(x_rs))
    sums = [0] * width
    for i in range(m):
        s = shares[i]
        if T_is_field:
            s = s.value
        # type(s) is int or gfpx.Polynomial
        for r in range(width):
            sums[r] += s * vector[r][i]
    for r in range(width):
        sums[r] = field(sums[r])
    if isinstance(x_rs, tuple):
        return sums[0]

    return sums