#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
def bincount2d(arr, bins=None):
    if bins is None:
        bins = np.max(arr) + 1
    count = np.zeros(shape=[len(arr), bins], dtype=np.int64)
    indexing = np.arange(len(arr))
    for col in arr.T:
        count[indexing, col] += 1
    return count


t = np.array([[1,2,3],[4,5,6],[3,2,2]], dtype=np.int64)
print(bincount2d(t))