def longest_string(arr):
    length = max([len(x) for x in arr])
    res = [x for x in arr if len(x) == length]
    return res