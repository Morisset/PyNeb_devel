# This module defines a class that handles I/O using
# Fortran-compatible format specifications.
#
#
# Warning: Fortran formatting is a complex business and I don't
# claim that this module works for anything complicated. It knows
# only the most frequent formatting options. Known limitations:
#
# 1) Only A, D, E, F, G, I, and X formats are supported (plus string constants
#    for output).
# 2) No direct support for complex numbers. You have to split them into
#    real and imaginary parts before output, and for input you get
#    two float numbers anyway.
# 3) No overflow check. If an output field gets too large, it will
#    take more space, instead of being replaced by stars.
#
#
# Written by Konrad Hinsen <hinsen@cnrs-orleans.fr>
# With contributions from Andreas Prlic <andreas@came.sbg.ac.at>
# last revision: 2006-6-23
#

"""
Fortran-style formatted input/output

This module provides two classes that aid in reading and writing
Fortran-formatted text files.

Examples::

  Input::

    >>>s = '   59999'
    >>>format = FortranFormat('2I4')
    >>>line = FortranLine(s, format)
    >>>print line[0]
    >>>print line[1]

  prints::

    >>>5
    >>>9999


  Output::

    >>>format_ = FortranFormat('2D15.5')
    >>>line = FortranLine([3.1415926, 2.71828], format_)
    >>>print str(line)

  prints::

    '3.14159D+00    2.71828D+00'
"""

#
# The class FortranLine represents a single line of input/output,
# which can be accessed as text or as a list of items.
#
class FortranLine:

    """Fortran-style record in formatted files

    FortranLine objects represent the content of one record of a
    Fortran-style formatted file. Indexing yields the contents as
    Python objects, whereas transformation to a string (using the
    built-in function 'str') yields the text representation.

    Restrictions:

      1. Only A, D, E, F, G, I, and X formats are supported (plus string
         constants for output).

      2. No direct support for complex numbers; they must be split into
         real and imaginary parts before output.

      3. No overflow check. If an output field gets too large, it will
         take more space, instead of being replaced by stars according
         to Fortran conventions.
    """

    def __init__(self, line, format_, length = 80):
        """
        @param data: either a sequence of Python objects, or a string
                     formatted according to Fortran rules

        @param format_: either a Fortran-style format string, or a
                       L{FortranFormat} object. A FortranFormat should
                       be used when the same format string is used repeatedly,
                       because then the rather slow parsing of the string
                       is performed only once.

        @param length: the length of the Fortran record. This is relevant
                       only when data is a string; this string is then
                       extended by spaces to have the indicated length.
                       The default value of 80 is almost always correct.
        """
        if type(line) == type(''):
            self.text = line
            self.data = None
        else:
            self.text = None
            self.data = line
        if type(format_) == type(''):
            self.format_ = FortranFormat(format_)
        else:
            self.format_ = format_
        self.length = length
        if self.text is None:
            self._output()
        if self.data is None:
            self._input()

    def __len__(self):
        """
        @returns: the number of data elements in the record
        @rtype: C{int}
        """
        return len(self.data)

    def __getitem__(self, i):
        """
        @param i: index
        @type i: C{int}
        @returns: the ith data element
        """
        return self.data[i]

    def __getslice__(self, i, j):
        """
        @param i: start index
        @type i: C{int}
        @param j: end index
        @type j: C{int}
        @returns: a list containing the ith to jth data elements
        """
        return self.data[i:j]

    def __str__(self):
        """
        @returns: a Fortran-formatted text representation of the data record
        @rtype: C{str}
        """
        return self.text

    def isBlank(self):
        """
        @returns: C{True} if the line contains only whitespace
        @rtype: C{bool}
        """
        return len(self.text.strip()) == 0

    def _input(self):
        text = self.text
        if len(text) < self.length: text = text + (self.length-len(text))*' '
        self.data = []
        for field in self.format_:
            l = field[1]
            s = text[:l]
            text = text[l:]
            type_ = field[0]
            value = None
            if type_ == 'A':
                value = s
            elif type_ == 'I':
                s = s.strip()
                if len(s) == 0:
                    value = 0
                else:
                    # by AP
                    # sometimes a line does not match to expected format,
                    # e.g.: pdb2myd.ent.Z chain: - model: 0 : CONECT*****
                    # catch this and set value to None
                    try:
                        value = int(s)
                    except:
                        value = None
            elif type_ == 'D' or type_ == 'E' or type_ == 'F' or type_ == 'G':
                s = s.strip().lower()
                n = s.find('d')
                if n >= 0:
                    s = s[:n] + 'e' + s[n+1:]
                if len(s) == 0:
                    value = 0.
                else:
                    try:
                        value = float(s)
                    except:
                        value = None
            if value is not None:
                self.data.append(value)

    def _output(self):
        data = self.data
        self.text = ''
        for field in self.format_:
            type_ = field[0]
            if type_ == "'":
                self.text = self.text + field[1]
            elif type_ == 'X':
                self.text = self.text + field[1]*' '
            else: # fields that use input data
                length = field[1]
                if len(field) > 2: fraction = field[2]
                value = data[0]
                data = data[1:]
                if type_ == 'A':
                    self.text = self.text + (value+length*' ')[:length]
                else: # numeric fields
                    if value is None:
                        s = ''
                    elif type_ == 'I':
                        s = repr(value)
                    elif type_ == 'D':
                        s = ('%'+repr(length)+'.'+repr(fraction)+'e') % value
                        n = s.find('e')
                        s = s[:n] + 'D' + s[n+1:]
                    elif type_ == 'E':
                        s = ('%'+repr(length)+'.'+ repr(fraction)+'e') % value
                    elif type_ == 'F':
                        s = ('%'+repr(length)+'.'+repr(fraction)+'f') % value
                    elif type_ == 'G':
                        s = ('%'+repr(length)+'.'+repr(fraction)+'g') % value
                    else:
                        raise ValueError('Not yet implemented')
                    s = s.upper()
                    self.text = self.text + ((length*' ')+s)[-length:]
        self.text = self.text.strtrip()

#
# The class FortranFormat represents a format specification.
# It ought to work for correct specifications, but there is
# little error checking.
#
class FortranFormat:

    """
    Parsed Fortran-style format string

    FortranFormat objects can be used as arguments when constructing
    FortranLine objects instead of the plain format string. If a
    format string is used more than once, embedding it into a FortranFormat
    object has the advantage that the format string is parsed only once.
    """

    def __init__(self, format_, nested = False):
        """
        @param format_: a Fortran format specification
        @type format_: C{str}
        @param nested: I{for internal use}
        """
        fields = []
        format_ = format_.strip()
        while format_ and format_[0] != ')':
            n = 0
            while format_[0] in '0123456789':
                n = 10*n + int(format_[0])
                format_ = format_[1:]
            if n == 0: n = 1
            type_ = format_[0].upper()
            if type_ == "'":
                eof = format_.find("'", 1)
                text = format_[1:eof]
                format_ = format_[eof+1:]
            else:
                format_ = format_[1:].strip()
            if type_ == '(':
                subformat_ = FortranFormat(format_, 1)
                fields = fields + n*subformat_.fields
                format_ = subformat_.rest
                eof = format_.find(',')
                if eof >= 0:
                    format_ = format_[eof+1:]
            else:
                eof = format_.find(',')
                if eof >= 0:
                    field = format_[:eof]
                    format_ = format_[eof+1:]
                else:
                    eof = format_.find(')')
                    if eof >= 0:
                        field = format_[:eof]
                        format_ = format_[eof+1:]
                    else:
                        field = format_
                        format_ = ''
                if type_ == "'":
                    field = (type_, text)
                else:
                    dot = field.find('.')
                    if dot > 0:
                        length = int(float(field[:dot]))
                        fraction = int(float(field[dot+1:]))
                        field = (type_, length, fraction)
                    else:
                        if field:
                            length = int(field)
                        else:
                            length = 1
                        field = (type_, length)
                fields = fields + n*[field]
        self.fields = fields
        if nested:
            self.rest = format_

    def __len__(self):
        return len(self.fields)

    def __getitem__(self, i):
        return self.fields[i]


# Test code

if __name__ == '__main__':
    f = FortranFormat("'!!',D10.3,F10.3,G10.3,'!!'")
    l = FortranLine([1.5707963, 3.14159265358, 2.71828], f)
    print((str(l)))
    f = FortranFormat("F12.0")
    l = FortranLine('2.1D2', f)
    print((l[0]))
