from string import ascii_lowercase
alpha1 = alphabet[0: 1]
alpha2 = alphabet[1:]
new_alphabet = alpha2 + alpha1

def shift(string):
    alphabet_order = {v: i for i, v in enumerate(alphabet)}
    old_order = [alphabet_order[x] for x in string.lower()]
    new_alphabet_order = {v: i for i, v in enumerate(new_alphabet)}
    new_order = [new_alphabet_order[x] for x in string.lower()]
    return "".join([new_alphabet[x] for x in old_order])