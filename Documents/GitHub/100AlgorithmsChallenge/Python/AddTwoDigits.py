def add_digits(n):
    digits = sum(list(map(int, [x for x in str(n)])))
    return digits