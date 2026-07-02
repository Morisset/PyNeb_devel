import pickle
import sys

def save(file_, *args, **kwargs):
    """
    Save the value of some data in a file.
    Usage: save('misdatos.pypic','a',b=b)
    Positional arguments are variable names looked up in the caller's namespace.

    """
    dict_ = kwargs
    if args:
        caller = sys._getframe(1)
        for name in args:
            if name in caller.f_locals:
                dict_[name] = caller.f_locals[name]
            elif name in caller.f_globals:
                dict_[name] = caller.f_globals[name]
            else:
                raise NameError('save: variable {0} not found in the caller namespace'.format(name))
    with open(file_, "wb") as f:
        pickle.dump(dict_, f, protocol=2)

def restore(file_):
    """
    Read data saved with save function.
    Usage: datos = restore('misdatos.pypic')

    """
    with open(file_, "rb") as f:
        result = pickle.load(f)
    return result
