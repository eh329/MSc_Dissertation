def adjacent_elements_products(array):
    products = [array[x] * array[x + 1] for x in range(0, len(array) - 1)]
    return max(products))