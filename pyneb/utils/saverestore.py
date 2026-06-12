import pickle

def save(file_, *args, **kwargs):
    """
    Save the value of some data in a file.
    Usage: save('misdatos.pypic','a',b=b)

    """
    dict = kwargs
    for name in args:
        dict[name] = eval(name)
    with open(file_, "wb") as f:
        pickle.dump(dict, f, protocol=2)

def restore(file_):
    """
    Read data saved with save function.
    Usage: datos = restore('misdatos.pypic')

    """
    with open(file_, "rb") as f:
        result = pickle.load(f)
    return result
