from . import write_util

__all__ = ['write_util']


def list_contains(List1, List2):
    check = False

    # Iterate in the 1st list
    for m in List1:

        # Iterate in the 2nd list
        for n in List2:

            # if there is a match
            if m == n:
                check = True
                return check

    return check

def findclosest(list, value):
    a = min((abs(x - value), x, i) for i, x in enumerate(list))
    return a[2], a[1]
