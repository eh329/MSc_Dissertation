def add_border(ls):
    new_ls = ["*" + x + "*" for x in ls]
    new_ls.insert(0, "*****")
    new_ls.insert(len(new_ls), "*****")
    return new_ls