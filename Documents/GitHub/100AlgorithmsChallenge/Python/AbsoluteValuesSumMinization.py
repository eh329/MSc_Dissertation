def absolute_values_sum(arr):
    
    arr.sort()
    length = len(arr)
    median = length // 2

    if length % 2 == 0:
        return arr[median - 1]

    return arr[median]
