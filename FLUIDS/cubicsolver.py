import math


def cubicsolver(a: float, b: float, c: float) -> list:
    x0 = None
    x1 = None
    x2 = None

    q = (a * a - 3 * b)
    r = (2 * a * a * a - 9 * a * b + 27 * c)

    Q = q / 9

    R = r / 54

    Q3 = Q * Q * Q

    R2 = R * R

    CR2 = 729 * r * r

    CQ3 = 2916 * q * q * q

    if R == 0 and Q == 0:
        x0 = - a / 3
        x1 = - a / 3
        x2 = - a / 3
        return [3, x0, x1, x2]
    elif CR2 == CQ3:
        """ /* this test is actually R2 == Q3, written in a form suitable
         for exact computation with integers */

      /* Due to finite precision some double roots may be missed, and
         considered to be a pair of complex roots z = x +/- epsilon i
         close to the real axis. */"""
        sqrtQ = math.sqrt(Q)

        if R > 0:
            x0 = -2 * sqrtQ - a / 3
            x1 = sqrtQ - a / 3
            x2 = sqrtQ - a / 3
        else:
            x0 = - sqrtQ - a / 3
            x1 = - sqrtQ - a / 3
            x2 = 2 * sqrtQ - a / 3
        return [3, x0, x1, x2]
    elif R2 < Q3:
        # {
        # double
        sgnR = (1 if R >= 0 else -1)
        # double
        ratio = sgnR * math.sqrt(R2 / Q3);
        # double
        theta = math.acos(ratio);
        # double
        norm = -2 * math.sqrt(Q);
        x0 = norm * math.cos(theta / 3) - a / 3;
        x1 = norm * math.cos((theta + 2.0 * math.pi) / 3) - a / 3;
        x2 = norm * math.cos((theta - 2.0 * math.pi) / 3) - a / 3;

        """/* Sort *x0, *x1, *x2 into increasing order */"""

        if x0 > x1:
            x0, x1 = x1, x0

        if x1 > x2:
            x1, x2 = x2, x1
            if x0 > x1:
                x0, x1 = x1, x0
        return [3, x0, x1, x2]
    else:
        sgnR = (1 if R >= 0 else -1)
        A = -sgnR * math.pow(math.fabs(R) + math.sqrt(R2 - Q3), 1.0 / 3.0)
        B = Q / A
        x0 = A + B - a / 3
        return [1, x0, x1, x2]
