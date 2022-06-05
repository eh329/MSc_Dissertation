def alternate(arr):
    group1 = sum([v for i, v in enumerate(arr) if  i % 2 == 0])
    group2 = sum([v for i, v in enumerate(arr) if i % 2 != 0])
    final = [group1, group2]
    return final