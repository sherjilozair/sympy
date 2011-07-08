def _slice_to_bounds(key, defmax): 
    """
        Takes slice or number and returns (min,max) for iteration
        Takes a default maxval to deal with the slice ':' which is (none, none)
    """
        lo, hi = 0, defmax
        if key.start is not None:
            if key.start >= 0:
                lo = key.start
            else:
                lo = defmax+key.start
        if key.stop is not None:
            if key.stop >= 0:
                hi = key.stop
            else:
                hi = defmax+key.stop
        return lo, hi
