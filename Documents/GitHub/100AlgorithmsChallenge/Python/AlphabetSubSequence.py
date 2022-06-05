from string import ascii_lowercase

def sub_sequence(string):

    alphabet = ascii_lowercase
    ordered_alphabet = {v: i for i, v in enumerate(alphabet)}
    order = [ordered_alphabet[x] for x in string.lower()]

    for i in range(0, len(order) - 1):
        if order[i + 1] > order[i]:
            continue
        
        else:
            return False
            
    return True