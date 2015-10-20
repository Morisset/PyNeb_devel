import sys
import re
import pdb
import os

IS_PYTHON3 = sys.version_info[0] >= 3

def _parser(tokens, version=None):
    # Parse the full edit descriptors
    eds = _parse_tokens(tokens, reversion=False, version=None)
    # Parse the edit descriptors used for the reversion of format control 
    # (F95 format 12.2.2)
    reversion_eds = _parse_tokens(tokens, reversion=True, version=None)
    return eds, reversion_eds

def _parse_tokens(tokens, reversion=False, version=None):
    # Remove outer parens is there are any
    tokens = _remove_outer_parens(tokens)
    # Get only the reversion tokens
    if reversion == True:
        tokens = _get_reversion_tokens(tokens)
    # First expand the parentheses
    tokens = _expand_parens(tokens)
    # Split on commas
    token_sets = _split_on_commas(tokens)
    # Split on ED9 (i.e. :)
    token_sets = _split_on_ed9(token_sets)
    # Split the ED10 (i.e. /)
    token_sets = _split_on_ed10(token_sets)
    # Split the ED8 (i.e. P edit descriptors)
    token_sets = _split_on_ed8(token_sets)
    # Process each set of edit descriptors
    eds = []
    for token_set in token_sets:
        # Assume first edit descriptor is the one to process
        ed_type = None
        ed_value = None
        for token in token_set:
            if token.type in ['ED1', 'ED2', 'ED3', 'ED4', 'ED5', 'ED6', 'ED7', 'ED8', 'ED9', 'ED10', 'QUOTED_STRING']:
                ed_type = token.type
                ed_value = token.value
                break
        # TODO: something more responsible here ...
        if ed_type is None:
            continue
        # If repeatable and first token is repeat number then cache
        repeat = None
        if ed_value in REPEATABLE_EDS and (token_set[0].type in ['NZUINT', 'UINT']):
            repeat = token_set[0].value
            token_set = token_set[1:]
        # Process the edit descriptor
        if ed_type == 'QUOTED_STRING':
            ed = _read_quoted_string(token_set)
        elif ed_type == 'ED1':
            ed = _read_ed1(token_set)
        elif ed_type == 'ED2':
            ed = _read_ed2(token_set)
        elif ed_type == 'ED3':
            ed = _read_ed3(token_set)
        elif ed_type == 'ED4':
            ed = _read_ed4(token_set)
        elif ed_type == 'ED5':
            ed = _read_ed5(token_set)
        elif ed_type == 'ED6':
            ed = _read_ed6(token_set)
        elif ed_type == 'ED7':
            ed = _read_ed7(token_set)
        elif ed_type == 'ED8':
            ed = _read_ed8(token_set)
        elif ed_type == 'ED9':
            ed = _read_ed9(token_set)
        elif ed_type == 'ED10':
            ed = _read_ed10(token_set)
        else:
            raise InvalidFormat('Could not identify edit descriptor in sequence $s' % str(token_set))
        # If there is a repeat number cached, then apply
        if repeat is not None:
            ed.repeat = repeat
        # Add to the list
        eds.append(ed)
    return eds


# Functions that group the tokens into sets

def _expand_parens(tokens):
    new_tokens = []
    get_tokens = iter(tokens)
    for t0 in get_tokens:
        if t0.type != 'LEFT_PARENS':
            new_tokens.append(t0)
        else:
            # Read in all tokens in subformat and recurse back to self
            paren_tokens = []
            nesting = 1
            while nesting > 0:
                try:
                    if IS_PYTHON3:
                        t1 = next(get_tokens)
                    else:
                        t1 = get_tokens.next()
                except StopIteration:
                    raise InvalidFormat('Open parens in format')
                if t1.type == 'LEFT_PARENS':
                    nesting = nesting + 1
                elif t1.type == 'RIGHT_PARENS':
                    nesting = nesting - 1
                paren_tokens.append(t1)
            # Remove last right paren
            paren_tokens = paren_tokens[:-1]
            # If repeated, then repeat the tokens accordingly
            if (len(new_tokens) > 0) and (new_tokens[-1].type in ['NZUINT', 'UINT']):
                repeat = new_tokens[-1].value
                # Remove the repeat NZUINT, UINT
                new_tokens = new_tokens[:-1]
                new_tokens.extend(repeat * (_expand_parens(paren_tokens) + [Token('COMMA', None)]))
            else:
                new_tokens.extend(_expand_parens(paren_tokens))
    return new_tokens


def _split_on_commas(tokens):
    token_sets = []
    set_buff = []
    for t0 in tokens:
        if t0.type == 'COMMA':
            token_sets.append(set_buff)
            set_buff = []
        else:
            set_buff.append(t0)
    token_sets.append(set_buff)
    return token_sets


def _split_on_ed9(token_sets):
    '''Splits on :'''
    new_token_sets = []
    for token_set in token_sets:
        if 'ED9' not in [t.type for t in token_set]:
            new_token_sets.append(token_set)
        else:
            buff = []
            for token in token_set:
                if token.type == 'ED9':
                    if len(buff) > 0:
                        new_token_sets.append(buff)
                        buff = []
                    new_token_sets.append([token])
                else:
                    buff.append(token)
            if len(buff) > 0:
                new_token_sets.append([token])
    return new_token_sets


def _split_on_ed10(token_sets):
    '''Splits on /'''
    new_token_sets = []
    for token_set in token_sets:
        # May have a repeat on the slash if preceded by a comma
        if (len(token_set) > 2) and ((token_set[0].type in ['UINT', 'NZUINT']) and (token_set[1].type == 'ED10')):
            new_token_sets.append(token_set[:2])
            token_set = token_set[2:]
        buff = []
        for token in token_set:
            if token.type == 'ED10':
                if len(buff) > 0:
                    new_token_sets.append(buff)
                    buff = []
                new_token_sets.append([token])
            else:
                buff.append(token)
        if len(buff) > 0:
            new_token_sets.append(buff)
    return new_token_sets


def _split_on_ed8(token_sets):
    '''Splits on ED8 (i.e. P edit descriptors)'''
    new_token_sets = []
    for token_set in token_sets:
        # Append to new list if no ED8
        if 'ED8' not in [t.type for t in token_set]:
            new_token_sets.append(token_set)
        # Otherwise split string on ED9
        elif (token_set[0].type in ['INT', 'UINT', 'NZUINT']) and (token_set[1].type == 'ED8'):
            new_token_sets.append(token_set[:2])
            new_token_sets.append(token_set[2:])
        else:
            raise InvalidFormat('P edit descriptor in invalid position')
    return new_token_sets

# Function to extract only the tokens for the reversion of control

def _get_reversion_tokens(tokens):
    reversion_tokens = []
    # Easier to work backwards
    nesting = None
    for token in tokens[::-1]:
        # End of loop condition
        if (nesting is not None) and (nesting < 1):
            # Parens may have a repeat number
            if token.type in ['UINT', 'NZUINT']:
                reversion_tokens.append(token)
            break
        # Read till the first right parens
        if token.type == 'RIGHT_PARENS':
            if nesting is None:
                nesting = 1
            else:
                nesting = nesting + 1
        elif token.type == 'LEFT_PARENS':
            if nesting is None:
                raise InvalidFormat('Unbalanced parens in format')
            else:
                nesting = nesting - 1
        reversion_tokens.append(token)
    # Tokens are in reverse order
    reversion_tokens.reverse()
    return reversion_tokens

# The functions that read particular edit descriptors sequences

def _read_quoted_string(tokens):
    # Of form X only
    type_string = ",".join([t.type for t in tokens])
    if type_string != "QUOTED_STRING":
        raise InvalidFormat('Token %s has invalid neighbouring token' % tokens[0])
    ed = QuotedString()  
    ed.char_string = tokens[0].value
    return ed

def _read_ed1(tokens):
    # Of form X only
    type_string = ",".join([t.type for t in tokens])
    if type_string != "ED1":
        raise InvalidFormat('Token %s has invalid neighbouring token' % tokens[0])
    ed = get_edit_descriptor_obj(tokens[0].value)   
    return ed

def _read_ed2(tokens):
    # Of form nX only
    type_string = ",".join([t.type for t in tokens])
    if type_string != "NZUINT,ED2":
        raise InvalidFormat('Token %s has invalid neighbouring token' % tokens[0])
    ed = get_edit_descriptor_obj(tokens[1].value)
    ed.num_chars = tokens[0].value
    return ed

def _read_ed3(tokens):
    # Of form Xn only
    type_string = ",".join([t.type for t in tokens])
    if type_string != "ED3,NZUINT":
        raise InvalidFormat('Token %s has invalid neighbouring token' % tokens[0])
    ed = get_edit_descriptor_obj(tokens[0].value)
    # L edit descriptor has a width rather than num_chars
    if hasattr(ed, 'width'):
        ed.width = tokens[1].value
    else:
        ed.num_chars = tokens[1].value
    return ed

def _read_ed4(tokens):
    # Of form X or Xn
    type_string = ",".join([t.type for t in tokens])
    if type_string in ["ED4", "ED4,NZUINT"] or \
        (ALLOW_ZERO_WIDTH_EDS and (type_string == "ED4,UINT")):
        ed = get_edit_descriptor_obj(tokens[0].value)
        if len(tokens) > 1:
            ed.width = tokens[1].value
    else:
        raise InvalidFormat('Token %s has invalid neighbouring token' % tokens[0])
    return ed

def _read_ed5(tokens):
    # Of form Xn.m only
    type_string = ",".join([t.type for t in tokens])
    if type_string in ["ED5,NZUINT,DOT,UINT", "ED5,NZUINT,DOT,NZUINT"] or \
      (ALLOW_ZERO_WIDTH_EDS and (type_string in \
      ["ED5,UINT,DOT,UINT", "ED5,UINT,DOT,NZUINT"])):
        ed = get_edit_descriptor_obj(tokens[0].value)
        ed.width = tokens[1].value
        ed.decimal_places = tokens[3].value
    else:
        raise InvalidFormat('%s has invalid neighbouring token' % tokens[0])
    return ed

def _read_ed6(tokens):
    # Of form Xn or Xn.m
    type_string = ",".join([t.type for t in tokens])
    if type_string == "ED6,NZUINT" or \
        (ALLOW_ZERO_WIDTH_EDS and (type_string == "ED6,UINT")):
        ed = get_edit_descriptor_obj(tokens[0].value)
        ed.width = tokens[1].value
        ed.min_digits = None
    elif type_string in ["ED6,NZUINT,DOT,UINT", "ED6,NZUINT,DOT,NZUINT"] or \
      (ALLOW_ZERO_WIDTH_EDS and (type_string in \
      ["ED6,UINT,DOT,UINT", "ED6,UINT,DOT,NZUINT"])):
        ed = get_edit_descriptor_obj(tokens[0].value)
        ed.width = tokens[1].value
        ed.min_digits = tokens[3].value
    else:
        raise InvalidFormat('%s has invalid neighbouring token' % tokens[0])
    return ed

def _read_ed7(tokens):
    # Of form Xn.m or Xn.mEe
    type_string = ",".join([t.type for t in tokens])
    if type_string in ["ED7,NZUINT,DOT,UINT", "ED7,NZUINT,DOT,NZUINT"] or \
      (ALLOW_ZERO_WIDTH_EDS and (type_string in \
      ["ED7,UINT,DOT,UINT", "ED7,UINT,DOT,NZUINT"])):
        ed = get_edit_descriptor_obj(tokens[0].value)
        ed.width = tokens[1].value
        ed.decimal_places = tokens[3].value
        ed.exponent = None
    elif type_string in ['ED7,NZUINT,DOT,NZUINT,ED7,NZUINT', \
            'ED7,NZUINT,DOT,NZUINT,ED7,UINT', \
            'ED7,NZUINT,DOT,NZUINT,ED7,INT', \
            'ED7,NZUINT,DOT,UINT,ED7,NZUINT', \
            'ED7,NZUINT,DOT,UINT,ED7,UINT', \
            'ED7,NZUINT,DOT,UINT,ED7,INT'] or \
        (ALLOW_ZERO_WIDTH_EDS and (type_string in \
            ['ED7,UINT,DOT,NZUINT,ED7,NZUINT', \
            'ED7,UINT,DOT,NZUINT,ED7,UINT', \
            'ED7,UINT,DOT,NZUINT,ED7,INT', \
            'ED7,UINT,DOT,UINT,ED7,NZUINT', \
            'ED7,UINT,DOT,UINT,ED7,UINT', \
            'ED7,UINT,DOT,UINT,ED7,INT'])):
        ed = get_edit_descriptor_obj(tokens[0].value)
        ed.width = tokens[1].value
        ed.decimal_places = tokens[3].value
        ed.exponent = tokens[5].value
    else:
        raise InvalidFormat('%s has invalid neighbouring token' % tokens[0])
    return ed

def _read_ed8(tokens):
    # Of form kX only, where k is a signed integer, may omit comma if followed
    # by Type 5 or 7 edit descriptor
    type_string = ",".join([t.type for t in tokens])
    if type_string in ["NZUINT,ED8", "UINT,ED8", "INT,ED8"]:
        ed = get_edit_descriptor_obj(tokens[1].value)
        ed.scale = tokens[0].value
    else:
        raise InvalidFormat('%s has invalid neighbouring token' % tokens[0])
    return ed

def _read_ed9(tokens):
    # Of form X only, may omit comma either side
    type_string = ",".join([t.type for t in tokens])
    if type_string == "ED9":
        ed = get_edit_descriptor_obj(tokens[0].value)
    else:
        raise InvalidFormat('%s has invalid neighbouring token' % tokens[0])
    return ed

def _read_ed10(tokens):
    # Of form X only, may omit following comma and leading comma if no repeat
    type_string = ",".join([t.type for t in tokens])
    if type_string == "ED10":
        ed = get_edit_descriptor_obj(tokens[0].value)
    else:
        raise InvalidFormat('%s has invalid neighbouring token' % tokens[0])
    return ed


# Functions that pre-process the token list

def _remove_outer_parens(tokens):
    # Finally, remove outer parens is there are any
    if (tokens[0].type == 'LEFT_PARENS') and (tokens[-1].type == 'RIGHT_PARENS'):
        tokens = tokens[1:-1]
    return tokens


# Run some tests if run as a script

#if __name__ == '__main__':
#    import doctest
#    import os
#    from _lexer import lexer
#    globs = {'lexer' : lexer, 'parser' : parser}
#    # Need to normalize whitespace since pasting into VIM converts tabs to
#    # spaces
#    doctest.testfile(os.path.join('tests', 'parser_test.txt'), \
#        globs=globs, optionflags=doctest.NORMALIZE_WHITESPACE)

# Some lexer tokens to look out for
DIGITS = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']
SIGNS = ['+', '-']
COMMA = [',']
DOT = ['.']
WHITESPACE = [' ', '\t', '\n']
QUOTE_CHARS = ['"', "'"]
DOUBLE_EDIT_DESCRIPTORS = ['EN', 'ES', 'TR', 'TL', 'BN', 'BZ', 'SP', 'SS']
SINGLE_EDIT_DESCRIPTORS = ['A', 'B', 'D', 'E', 'F', 'G', 'I', 'L', 'O', 'P', 'S', 'T', 'X', 'Z', ':', '/']
H_EDIT_DESCRIPTOR = ['H']
LEFT_PARENS = ['(']
RIGHT_PARENS = [')']
COLON = ':'
SLASH = '/'


def _lexer(format):
    '''Lex the FORTRAN format statement into tokens'''
    tokens = []
    s = -1
    h_chars = None
    while True:
        # Get the next set of characters
        s = s + 1
        c0, c1, c2 = _get_chars(format, s)
        # If at end of format, end it all - aieee!
        if c0 is None:
            break
        # Read in an H edit descriptor string
        elif h_chars is not None:
            buff = format[s:s+h_chars]
            tokens.append(Token('QUOTED_STRING', buff))
            s = s + (h_chars - 1)
            h_chars = None
        # Skip whitespace
        elif c0 in WHITESPACE:
            continue
        # Read in a quoted string
        elif c0 in QUOTE_CHARS:
            buff = ''
            delim = c0
            while True:
                s = s + 1
                c0, c1, c2 = _get_chars(format, s)
                # Check if an escaped delimiter
                if (c0 == delim) and (c1 == delim):
                    s = s + 1
                    buff = buff + delim
                elif (c0 == delim):
                    break
                elif c0 is None:
                    # Premature end of format
                    raise InvalidFormat('Premature end of quoted string in format')
                else:
                    buff = buff + c0
            tokens.append(Token('QUOTED_STRING', buff))
        # Read in an integer
        elif c0 in DIGITS + SIGNS:
            # Check sign followed by at least one digit
            if (c0 in SIGNS) and (c1 not in DIGITS):
                # TODO: Is whitesapce allowed between sign and digit?
                raise InvalidFormat("Orphaned sign '%s' with no digits at position %d" % (c0, s))
            buff = c0
            while True:
                s = s + 1
                c0, c1, c2 = _get_chars(format, s)
                if (c0 not in DIGITS) or (c0 is None):
                    break
                else:
                    buff = buff + c0
            s = s - 1
            val = int(buff)
            if buff[0] in SIGNS:
                tokens.append(Token('INT', val))
            elif val == 0:
                tokens.append(Token('UINT', val))
            else:
                tokens.append(Token('NZUINT', val))
        # Read in a comma
        elif c0 in COMMA:
            tokens.append(Token('COMMA', None))
        # Read in a dot
        elif c0 in DOT:
            tokens.append(Token('DOT', None))
        # Read in double lettered edit descriptors
        elif (c1 is not None) and ((c0 + c1).upper() in DOUBLE_EDIT_DESCRIPTORS):
            ed_type = _get_ed_type((c0 + c1).upper())
            tokens.append(Token(ed_type, (c0 + c1).upper()))
            s = s + 1
        # Read in an H edit descriptor
        elif c0.upper() in H_EDIT_DESCRIPTOR:
            if (len(tokens) > 0) and (tokens[-1].type in ('NZUINT', 'UINT')):
                h_chars = tokens[-1].value
                tokens = tokens[:-1]
            else:
                raise InvalidFormat("Missing H descriptor number argument at position %d" % s)
        # Read in single lettered edit descriptors
        elif c0.upper() in SINGLE_EDIT_DESCRIPTORS:
            ed_type = _get_ed_type(c0.upper())
            tokens.append(Token(ed_type, c0.upper()))
        # Read in left parens
        elif c0 in LEFT_PARENS:
            tokens.append(Token('LEFT_PARENS', None))
        # Read in right parens
        elif c0 in RIGHT_PARENS:
            tokens.append(Token('RIGHT_PARENS', None))
        else:
            raise InvalidFormat('Character %s not recognised at position %d' % (c0, s))
    return tokens

def _get_ed_type(ed_string):
    if ed_string in ED1:
        ed_type = 'ED1'
    elif ed_string in ED2:
        ed_type = 'ED2'
    elif ed_string in ED3:
        ed_type = 'ED3'
    elif ed_string in ED4:
        ed_type = 'ED4'
    elif ed_string in ED5:
        ed_type = 'ED5'
    elif ed_string in ED6:
        ed_type = 'ED6'
    elif ed_string in ED7:
        ed_type = 'ED7'
    elif ed_string in ED8:
        ed_type = 'ED8'
    elif ed_string in ED9:
        ed_type = 'ED9'
    elif ed_string in ED10:
        ed_type = 'ED10'
    else:
        ed_type = None
    return ed_type

def _get_chars(format, s):
    try:
        c0 = format[s]
    except IndexError:
        c0 = None
    try:
        c1 = format[s+1]
    except IndexError:
        c1 = None
    try:
        c2 = format[s+2]
    except IndexError:
        c2 = None
    return (c0, c1, c2)


class InvalidFormat(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)


class Token(object):
    def __init__(self, type, value):
        self.type = type
        self.value = value
    def __repr__(self):
        return "\n  Token: type=%s,\tvalue=%s" % (self.type, str(self.value))

# Do some testing when run as a module

#if __name__ == '__main__':
#    import doctest
#    import os
#    globs = {'lexer' : lexer}
#    # Need to normalize whitespace since pasting into VIM converts tabs to
#    # spaces
#    doctest.testfile(os.path.join('tests', 'lexer_test.txt'), \
#        globs=globs, optionflags=doctest.NORMALIZE_WHITESPACE)

class InvalidFormat(Exception):
    pass

def get_edit_descriptor_obj(name):
    '''Returns a new object instance from a string'''
    name = name.upper()
    if name == 'A':
        return A()
    elif name == 'B':
        return B()
    elif name == 'BN':
        return BN()
    elif name == 'BZ':
        return BZ()
    elif name == ':':
        return Colon()
    elif name == 'D':
        return D()
    elif name == 'E':
        return E()
    elif name == 'EN':
        return EN()
    elif name == 'ES':
        return ES()
    elif name == 'F':
        return F()
    elif name == 'G':
        return G()
    elif name == 'H':
        return H()
    elif name == 'I':
        return I()
    elif name == 'L':
        return L()
    elif name == 'O':
        return O()
    elif name == 'P':
        return P()
    elif name =='S':
        return S()
    elif name == '/':
        return Slash()
    elif name == 'SP':
        return SP()
    elif name == 'SS':
        return SS()
    elif name == 'T':
        return T()
    elif name == 'TL':
        return TL()
    elif name == 'TR':
        return TR()
    elif name == 'X':
        return X()
    elif name == 'Z':
        return Z()
    else:
        raise InvalidFormat('Expected an edit descriptor, got %s' % name)

# All the tokens defined in the F77 specification unless specified

class A(object):
    def __init__(self):
        self.repeat = None
        self.width = None
    def __repr__(self):
        return '<A repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + '>'

class QuotedString(object):
    def __init__(self, char_string=None):
        self.char_string = char_string
    def get_width(self):
        return len(self.char_string)
    width = property(get_width)
    def __repr__(self):
        return '<QuotedString char_string=' + str(self.char_string) + '>'

# Only in F95
class B(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.min_digits = None
    def __repr__(self):
        return '<B repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' min_digits=' + str(self.min_digits) + '>'
    
class BN(object):
    def __init__(self):
        pass
    def __repr__(self):
        return '<BN>'

class BZ(object):
    def __init__(self):
        pass
    def __repr__(self):
        return '<BZ>'

class Colon(object):
    def __init__(self):
        pass
    def __repr__(self):
        return '<Colon>'
    
class D(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.decimal_places = None
    def __repr__(self):
        return '<D repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' decimal_places=' + str(self.decimal_places) + '>'

class E(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.decimal_places = None
        self.exponent = None
    def __repr__(self):
        return '<E repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' decimal_places=' + str(self.decimal_places) + \
                ' exponent=' + str(self.exponent) + '>'
    
# Only in F95
class EN(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.decimal_places = None
        self.exponent = None
    def __repr__(self):
        return '<EN repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' decimal_places=' + str(self.decimal_places) + \
                ' exponent=' + str(self.exponent) + '>'

# Only in F95
class ES(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.decimal_places = None
        self.exponent = None
    def __repr__(self):
        return '<ES repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' decimal_places=' + str(self.decimal_places) + \
                ' exponent=' + str(self.exponent) + '>'

class F(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.decimal_places = None
    def __repr__(self):
        return '<F repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' decimal_places=' + str(self.decimal_places) + '>'
    
class FormatGroup(object):
    pass

class G(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.decimal_places = None
        self.exponent = None
    def __repr__(self):
        return '<G repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' decimal_places=' + str(self.decimal_places) + \
                ' exponent=' + str(self.exponent) + '>'

# Only in F77
class H(object):
    def __init__(self):
        self.num_chars = None
        self.char_string = None
    def __repr__(self):
        return '<H num_chars=' + str(self.num_chars) + \
                ' char_string=' + str(self.char_string) + '>'
    
class I(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.min_digits = None
    def __repr__(self):
        return '<I repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' min_digits=' + str(self.min_digits) + '>'
    
class L(object):
    def __init__(self):
        self.repeat = None
        self.width = None
    def __repr__(self):
        return '<L repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + '>'

# Only in F95
class O(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.min_digits = None
    def __repr__(self):
        return '<O repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' min_digits=' + str(self.min_digits) + '>'

class P(object):
    def __init__(self):
        self.scale = None
    def __repr__(self):
        return '<P scale=' + str(self.scale) + '>'
    
class S(object):
    def __init__(self):
        pass
    def __repr__(self):
        return '<S>'
    
class Slash(object):
    def __init__(self):
        self.repeat = None
        pass
    def __repr__(self):
        return '<Slash repeat=' + str(self.repeat) + '>'
    
class SP(object):
    def __init__(self):
        pass
    def __repr__(self):
        return '<SP>'
    
class SS(object):
    def __init__(self):
        pass
    def __repr__(self):
        return '<SS>'
    
class T(object):
    def __init__(self):
        self.num_chars = None
    def __repr__(self):
        return '<T num_chars=' + str(self.num_chars) + '>'
    
class TL(object):
    def __init__(self):
        self.num_chars = None
    def __repr__(self):
        return '<TL num_chars=' + str(self.num_chars) + '>'
    
class TR(object):
    def __init__(self):
        self.num_chars = None
    def __repr__(self):
        return '<TR num_chars=' + str(self.num_chars) + '>'

class X(object):
    def __init__(self):
        self.num_chars = None
    def __repr__(self):
        return '<X num_chars=' + str(self.num_chars) + '>'

# Only in F95
class Z(object):
    def __init__(self):
        self.repeat = None
        self.width = None
        self.min_digits = None
    def __repr__(self):
        return '<Z repeat=' + str(self.repeat) + \
                ' width=' + str(self.width) + \
                ' min_digits=' + str(self.min_digits) + '>'

# Categorise the edit descriptors depnding on how they should be parsed

ED1 = ['BN', 'BZ', 'SP', 'SS', 'S'] # Of form X only
ED2 = ['X'] # Of form nX only
ED3 = ['T', 'TR', 'TL', 'L'] # Of form Xn only
ED4 = ['A'] # Of form X or Xn
ED5 = ['D', 'F'] # Of form Xn.m only
ED6 = ['B', 'I', 'O', 'Z'] # Of form Xn or Xn.m
ED7 = ['E', 'EN', 'ES', 'G'] # Of form Xn.m or Xn.mEe
ED8 = ['P'] # Of form kX only, where k is a signed integer, may omit comma if followed by Type 5 or 7 edit descriptor
ED9 = [':'] # Of form X only, may omit comma either side
ED10 = ['/'] # Of form X only, may omit following comma and leading comma if no repeat
REPEATABLE_EDS = ['L', 'A', 'D', 'F', 'B', 'I', 'O', 'Z', 'E', 'EN', 'ES', 'G', '/']
OUTPUT_EDS = (L, A, D, F, B, I, O, Z, E, EN, ES, G)
CONTROL_EDS = (BN, BZ, P, SP, SS, S, X, T, TR, TL, Colon, Slash)
NON_REVERSION_EDS = (P, S, SP, SS, BN, BZ)
ALL_ED = ED1 + ED2 + ED3 + ED4 + ED5 + ED6 + ED7 + ED8 + ED9 + ED10


    
'''
Miscellaneous functions, classes etc
'''

class has_next_iterator(object):
  '''
  A wrapper class for iterators so that the .has_next() method is implemented

  See - http://stackoverflow.com/questions/1966591/hasnext-in-python-iterators
  '''
  def __init__(self, it):
    self.it = iter(it)
    self._has_next = None
  def __iter__(self): return self
  def __next__(self):
    if self._has_next:
      result = self._the_next
    else:
      if IS_PYTHON3:
        result = next(self.it)
      else:
        result = self.it.next()
    self._has_next = None
    return result
  def next(self):
    if self._has_next:
      result = self._the_next
    else:
      if IS_PYTHON3:
        result = next(self.it)
      else:
        result = self.it.next()
    self._has_next = None
    return result
  def has_next(self):
    if self._has_next is None:
      try: 
        if IS_PYTHON3:
          self._the_next = next(self.it)
        else:
          self._the_next = self.it.next()
      except StopIteration: self._has_next = False
      else: self._has_next = True
    return self._has_next


def expand_edit_descriptors(eds):
    expanded_eds = []
    for ed in eds:
        if hasattr(ed, 'repeat') and (ed.repeat is not None):
            expanded_eds.extend(ed.repeat * [ed])
        else:
            expanded_eds.append(ed)
    return expanded_eds


# Should all edit descriptor values be returned even if they were not
# written to?
RET_WRITTEN_VARS_ONLY = False
# Should 'None' values be returned when no record is available to read
# from or the FORTRAN 'default'?
RET_UNWRITTEN_VARS_NONE = True
# The order in which edit desciptors are tried by default when G edit
# descriptor encountered on input
G_INPUT_TRIAL_EDS = ['F', 'L', 'A']
# Contrary to specification, many compilers allow zero width edit
# descriptors
ALLOW_ZERO_WIDTH_EDS = True
# Set the characters that separate the records
RECORD_SEPARATOR = '\n'

# The maximum size for an integer
if sys.version_info[0] >= 3:
    PROC_MAXINT = sys.maxsize
else:
    PROC_MAXINT = sys.maxint
# Processor dependant default for including leading plus or not
PROC_INCL_PLUS = False 
# Option to allow signed binary, octal and hex on input (not a FORTRAN feature)
PROC_ALLOW_NEG_BOZ = False
# Prcessor dependant padding character
PROC_PAD_CHAR = ' '
# Interpret blanks or jsut a negative as a zero, as in ifort behaviour
PROC_NEG_AS_ZERO = True
# Show a sign for zero?
PROC_SIGN_ZERO = False
PROC_MIN_FIELD_WIDTH = 46
PROC_DECIMAL_CHAR = '.'
G0_NO_BLANKS = False
PROC_NO_LEADING_BLANK = False
# The default value if BN, BZ edit descriptors are not specified
PROC_BLANKS_AS_ZEROS = False

def reset():
    global RET_WRITTEN_VARS_ONLY, RET_UNWRITTEN_VARS_NONE, PROC_INCL_PLUS, \
        PROC_ALLOW_NEG_BOZ, PROC_PAD_CHAR, PROC_NEG_AS_ZERO, PROC_SIGN_ZERO, \
        PROC_MIN_FIELD_WIDTH, PROC_DECIMAL_CHAR, G0_NO_BLANKS, \
        PROC_NO_LEADING_BLANK, PROC_BLANKS_AS_ZEROS, PROC_MAXINT, G_INPUT_TRIAL_EDS, \
        ALLOW_ZERO_WIDTH_EDS
    G_INPUT_TRIAL_EDS = ['F', 'L', 'A']
    if sys.version_info[0] >= 3:
        PROC_MAXINT = sys.maxsize
    else:
        PROC_MAXINT = sys.maxint
    RET_WRITTEN_VARS_ONLY = False
    RET_UNWRITTEN_VARS_NONE = True
    PROC_INCL_PLUS = False
    PROC_ALLOW_NEG_BOZ = False
    PROC_PAD_CHAR = ' '
    PROC_NEG_AS_ZERO = True
    PROC_SIGN_ZERO = False
    PROC_MIN_FIELD_WIDTH = 46
    PROC_DECIMAL_CHAR = '.'
    G0_NO_BLANKS = False
    PROC_NO_LEADING_BLANK = False
    PROC_BLANKS_AS_ZEROS = False
    ALLOW_ZERO_WIDTH_EDS = True

    
WIDTH_OPTIONAL_EDS = [A]
NON_WIDTH_EDS = [BN, BZ, P, SP, SS, S, X, T, TR, TL, Colon, Slash]
FORBIDDEN_EDS = [QuotedString, H]

# Some problems without pre written input vars:
#   Cannot say when reversion conditions are met
#   Cannot determine width of A edit descriptor
#   Cannot determine complex input
#   Cannot determine proper input for G edit descriptors


def _input(eds, reversion_eds, records, num_vals=None):

    state = { \
        'position' : 0,
        'scale' : 0,
        'incl_plus' : False,
        'blanks_as_zeros' : PROC_BLANKS_AS_ZEROS,
        # TODO: Implement halt if no more record input
        'halt_if_no_vals' : False,
        'exception_on_fail' : True,
    }

    # pdb.set_trace()

    for ed in eds + reversion_eds:
        if isinstance(ed, tuple(FORBIDDEN_EDS)):
            raise InvalidFormat("%d edit descriptr not permitted on input")

    # Expand repeated edit decriptors
    eds = expand_edit_descriptors(eds)
    reversion_eds = expand_edit_descriptors(reversion_eds)
    # Assume one-to-one correspondance between edit descriptors and output
    # values if number of output values is not defined 
    num_out_eds = 0
    for ed in eds:
        if isinstance(ed, OUTPUT_EDS):
            num_out_eds += 1
    num_rev_out_eds = 0
    if num_vals is None:
        num_vals = num_out_eds
    for ed in reversion_eds:
        if isinstance(ed, OUTPUT_EDS):
            num_rev_out_eds += 1

    
    # Will loop forever is no output edit descriptors
    if (num_out_eds == 0):
        return []
    # Will loop forever if no output eds in reversion format and is more values
    # requested than in the format
    if (num_vals > num_out_eds) and (num_rev_out_eds == 0):
        raise ValueError('Not enough output edit descriptors in reversion format to output %d values' % num_vals)

    # May need to process multiple records, down to a higher function to supply
    # appropriate string for format
    if not hasattr(records, 'next'):
        records = iter(re.split('\r\n|\r|\n', records))
    record = _next(records, None)
    if record is None:
        return [] 
    
    # if a_widths is not None:
    #     a_widths = itertools.cycle(a_widths)

    vals = []
    finish_up = False
    ed_ind = -1
    while True:
        ed_ind += 1
        # Signal to stop when Colon edit descriptor or end of format or end of
        # reversion format reached. Also not to output any more data
        if len(vals) >= num_vals:
            finish_up = True
        # Select the appropriate edit descriptor
        if ed_ind < len(eds):
            ed = eds[ed_ind]
        else:
            rev_ed_ind = (ed_ind - len(eds)) % len(reversion_eds)
            # Reversion begun and has been instructed to halt
            if finish_up and (rev_ed_ind == 0):
                break
            ed = reversion_eds[rev_ed_ind]

        if isinstance(ed, QuotedString):
            raise InvalidFormat('Cannot have string literal in an input format')
        elif isinstance(ed, BN):
            state['blanks_as_zeros'] = False
        elif isinstance(ed, BZ):
            state['blanks_as_zeros'] = True
        elif isinstance(ed, P):
            state['scale'] = ed.scale
        elif isinstance(ed, SP):
            state['incl_plus'] = True
        elif isinstance(ed, SS):
            state['incl_plus'] = False
        elif isinstance(ed, S):
            state['incl_plus'] = PROC_INCL_PLUS
        elif isinstance(ed, (X, TR)):
            state['position'] = min(state['position'] + ed.num_chars, len(record))
        elif isinstance(ed, TL):
            state['position'] = max(state['position'] - ed.num_chars, 0)
        elif isinstance(ed, T):
            if (ed.num_chars - 1) < 0:
                state['position'] = 0
            elif ed.num_chars > len(record):
                state['position'] = len(record)
            else:
                state['position'] = ed.num_chars - 1
        elif isinstance(ed, Slash):
            # End of record
            record = _next(records, None)
            state['position'] = 0
            if record is None:
                break
        elif isinstance(ed, Colon):
            # Break if input value satisfied
            if finish_up:
                break
        elif isinstance(ed, (Z, O, B, I)):
            val, state = read_integer(ed, state, record)
            vals.append(val)
        elif isinstance(ed, A):
            val, state = read_string(ed, state, record)
            vals.append(val)
        elif isinstance(ed, L):
            val, state = read_logical(ed, state, record)
            vals.append(val)
        elif isinstance(ed, (F, E, D, EN, ES)):
            val, state = read_float(ed, state, record)
            vals.append(val)
        elif isinstance(ed, G):
            # Difficult to know what wanted since do not know type of input variable
            # Use the G_INPUT_TRIAL_EDS variable to try the variables
            # until one sticks
            # n.b. vals and state do not get written to if
            # exception id raised
            resolved = False
            g_trial_eds = iter(G_INPUT_TRIAL_EDS)
            while not resolved:
                ed_name = _next(g_trial_eds, '')
                if ed_name.upper() in ('F', 'E', 'D', 'EN', 'ES'):
                    trial_ed = F()
                    trial_ed.width = ed.width
                    trial_ed.decimal_places = ed.decimal_places
                    # pdb.set_trace()
                    try:
                        val, state = read_float(trial_ed, state.copy(), record)
                        vals.append(val)
                        resolved = True
                    except ValueError:
                        continue
                elif ed_name.upper() in ('Z', 'O', 'B', 'I'):
                    trial_ed = globals()[ed_name]()
                    trial_ed.width = ed.width
                    trial_ed.min_digits = ed.decimal_places
                    try:
                        val, state = read_integer(trial_ed, state.copy(), record)
                        vals.append(val)
                        resolved = True
                    except ValueError:
                        continue
                elif ed_name.upper() in ('L'):
                    trial_ed = L()
                    trial_ed.width = ed.width
                    try:
                        val, state = read_logical(trial_ed, state.copy(), record)
                        vals.append(val)
                        resolved = True
                    except ValueError:
                        continue
                elif ed_name.upper() in ('A'):
                    trial_ed = A()
                    trial_ed.width = ed.width
                    try:
                        val, state = read_string(trial_ed, state.copy(), record)
                        vals.append(val)
                        resolved = True
                    except ValueError:
                        continue
                elif ed_name in ('G'):
                    raise ValueError('G edit descriptor not permitted in G_INPUT_TRIAL_EDS')
                else:
                    raise ValueError('Unrecognised trial edit descriptor string in G_INPUT_TRIAL_EDS')
                    
    if RET_WRITTEN_VARS_ONLY:
        vals = [val for val in vals if val is not None]
    return vals[:num_vals]

def _interpret_blanks(substr, state):
    # Save leading blanks
    len_str = len(substr)
    if state['blanks_as_zeros']:
        # TODO: Are tabs blank characters?
        substr = substr.replace(' ', '0')
    else:
        substr = substr.replace(' ', '')
    # If were blanks but have been stripped away, replace with a zero
    if len(substr) == 0 and (len_str > 0):
        substr = '0'
    return substr

def _get_substr(w, record, state):
    start = max(state['position'], 0)
    end = start + w
    # if end > len(record):
    #     substr = ''
    #     # TODO: test if no chars transmitted, then poition does not change
    #     w = 0
    # else:
    substr = record[start:end]
    state['position'] = min(state['position'] + w, len(record))
    return substr, state


def _next(it, default=None):
    try:
        if IS_PYTHON3:
            val = next(it)
        else:
            val = it.next()
    except StopIteration:
        val = default
    return val


def read_string(ed, state, record):
    if ed.width is None:
        # Will assume rest of record is fair game for the
        # unsized A edit descriptor
        ed.width = len(record) - state['position']
    substr, state = _get_substr(ed.width, record, state)
    val = substr.ljust(ed.width, PROC_PAD_CHAR)
    return (val, state)


def read_integer(ed, state, record):
    substr, state = _get_substr(ed.width, record, state)
    if ('-' in substr) and (not PROC_ALLOW_NEG_BOZ) and isinstance(ed, (Z, O, B)):
        if state['exception_on_fail']:
            raise ValueError('Negative numbers not permitted for binary, octal or hex')
        else:
            return (None, state)
    if isinstance(ed, Z):
        base = 16
    elif isinstance(ed, I):
        base = 10
    elif isinstance(ed, O):
        base = 8
    elif isinstance(ed, B):
        base = 2
    # If a negative is followed by blanks, Gfortran and ifort
    # interpret as a zero
    if re.match(r'^ *- +$', substr):
        substr = '0'
    # If a negative or negative and blanks, ifort interprets as
    # zero for an I edit descriptor
    if PROC_NEG_AS_ZERO and isinstance(ed, I) and re.match(r'^( *- *| +)$', substr):
        substr = '0'
    # If string is zero length (reading off end of record?),
    # interpret as zero so as to match what would be found in an
    # unwritten FORTRAN variable
    if substr == '':
        if RET_UNWRITTEN_VARS_NONE or RET_WRITTEN_VARS_ONLY:
            return (None, state)
        else:
            substr = '0'
    teststr = _interpret_blanks(substr, state)
    try:
        val = int(teststr, base)
    except ValueError:
        if state['exception_on_fail']:
            raise ValueError('%s is not a valid input for one of integer, octal, hex or binary' % substr)
        else:
            return (None, state)
    return (val, state)


def read_logical(ed, state, record):
    substr, state = _get_substr(ed.width, record, state)
    # Deal with case where there is no more input to read from
    if (substr == '') and (RET_UNWRITTEN_VARS_NONE or RET_WRITTEN_VARS_ONLY):
        return (None, state)
    # Remove preceding whitespace and take the first two letters as
    # uppercase for testing
    teststr = substr.upper().lstrip().lstrip('.')
    if len(teststr):
        teststr = teststr[0]
    else:
        # This is case where just a preceding period is read in
        raise ValueError('%s is not a valid boolean input' % substr)
    if teststr == 'T':
        val = True
    elif teststr == 'F':
        val = False
    else:
        if state['exception_on_fail']:
            raise ValueError('%s is not a valid boolean input' % substr)
        else:
            val = None
    return (val, state)


def read_float(ed, state, record):
    substr, state = _get_substr(ed.width, record, state)
    teststr = _interpret_blanks(substr, state)
    # When reading off end of record, get empty string,
    # interpret as 0
    if teststr == '':
        if RET_UNWRITTEN_VARS_NONE or RET_WRITTEN_VARS_ONLY:
            return (None, state)
        else:
            teststr = '0'
    # Python only understands 'E' as an exponential letter
    teststr = teststr.upper().replace('D', 'E')
    # Prepend an exponential letter if only a '-' or '+' denotes an exponent
    if 'E' not in teststr:
        teststr = teststr[0] + teststr[1:].replace('+', 'E+').replace('-', 'E-')
    # ifort allows '.' to be interpreted as 0
    if re.match(r'^ *\. *$', teststr):
        teststr = '0'
    # ifort allows '-' to be interpreted as 0
    if re.match(r'^ *- *$', teststr):
        teststr = '0'
    # ifort allows numbers to end with 'E', 'E+', 'E-' and 'D'
    # equivalents
    res = re.match(r'(.*)(E|E\+|E\-)$', teststr)
    if res:
        teststr = res.group(1)
    try:
        val = float(teststr)
    except ValueError:
        if state['exception_on_fail']:
            raise ValueError('%s is not a valid input as for an E, ES, EN or D edit descriptor' % substr)
        else:
            return (None, state)
    # Special cases: insert a decimal if none specified
    if ('.' not in teststr) and (ed.decimal_places is not None):
        val = val / 10 ** ed.decimal_places
    # Apply scale factor if exponent not supplied
    if 'E' not in teststr:
        val = val / 10 ** state['scale'] 
    return (val, state)

class FortranRecordReader(object):
    '''
    Generate a reader object for FORTRAN format strings

    Typical use case ...

    >>> header_line = FortranRecordReader('(A15, A15, A15)')
    >>> header_line.read('              x              y              z')
    ['              x', '              y', '              z']
    >>> line = FortranRecordReader('(3F15.3)')
    >>> line.read('          1.000          0.000          0.500')
    [1.0, 0.0, 0.5]
    >>> line.read('          1.100          0.100          0.600')
    [1.1, 0.1, 0.6]

    Note: it is best to create a new object for each format, changing the format
    causes the parser to reevalute the format string which is costly in terms of
    performance
    '''
    
    def __init__(self, format):
        self.format = format
        self._eds = []
        self._rev_eds = []
        self._parse_format()

    def __eq__(self, other):
        if isinstance(other, FortranRecordReader):
            return self.format == other.format
        else:
            return object.__eq__(self, other)

    def match(self, record):
        try:
            self.read(record)
        except RecordError:
            return False
        else:
            return True

    def read(self, record):
        '''
        Pass a string representing a FORTRAN record to obtain the relevent
        values
        '''
        return _input(self._eds, self._rev_eds, record)

    def get_format(self):
        return self._format
    def set_format(self, format):
        self._format = format
        self._parse_format()
    format = property(get_format, set_format)

    def _parse_format(self):
        self._eds, self._rev_eds = _parser(_lexer(self.format))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
