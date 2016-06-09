from datetime import datetime
from time import time

class PyNebError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)    
    
class my_logging(object):
    """
    Logging object that can log events with 3 different types (messages, warnings and error).
    Date and time are save as well as the calling procedure.
    events are printed in the standard output depending on the verbosity level.
    level: -1: Quiet and no log into log_ variable
             0: Quiet
             1: Only Errors
             2: Errors and Warnings
             3: very verbose, print Errors, Warnings and Messages 
             4: debug
    
    """
    
    def __init__(self, level=None, calling=None, file_=None, print_time=False, no_exit=True):
        """
        Parameters:
            - level [int] 0, 1, 2, 3, or 4
            - calling [str] description of the calling module
            - file_ [str] file to store the events (default is None).
            - print_time [boolean] if True, the event is printed with date-time (Default is False).

        """
        self.log_ = {}
        self.log_['messages'] = []
        self.log_['warnings'] = []
        self.log_['errors'] = []
        self.log_['timer'] = []
        self.log_['debug'] = []
        
        self.caller_max_size = 0
        self.level = level if level is not None else 2
        self.__start = time()
        if calling is not None:
            self.calling = calling
        else:
            self.calling = 'None'
        self.print_time = print_time
        self.open_file(file_)
        self.no_exit = no_exit
    
    def _pprint(self, mess_type, message, calling=None):
        if calling is None:
            calling = self.calling
        if self.print_time:
            print('{0} {1} @ {2}: {3}'.format(mess_type, calling, str(datetime.now()), message))
        else:
            print('{0} {1}: {2}'.format(mess_type, calling, message))
        
    def _add2log(self, mess_type, message, calling=None):
        if calling is None:
            calling = self.calling
        if self.level > -1:
            self.log_[mess_type].append((calling, str(datetime.now()), message))
            if len(calling) > self.caller_max_size:
                self.caller_max_size = len(calling)
            if self._to_file:
                self._file_out.write("{0} - {1} @ {2} : {3} \n".format(mess_type, calling, str(datetime.now()), message))

    def debug(self, message, calling=None):
        """
        method to add a debug message to the log                
        param:
            message [str] string to log
        """
        if self.level > 3:
            self._pprint('debug', message, calling=calling)
        self._add2log('debug', message, calling=calling)

    def message(self, message, calling=None):
        """
        method to add a message to the log
        param:
            message [str] string to log
        """
        if self.level > 2:
            self._pprint('    ', message, calling=calling)
        self._add2log('messages', message, calling=calling)

    def warn(self, message, calling=None):
        """
        method to add a warning to the log                
        param:
            message [str] string to log
        """
        if self.level > 1:
            self._pprint('warng', message, calling=calling)
        self._add2log('warnings', message, calling=calling)
        
    def error(self, message, calling=None, exception=None):
        """
        method to add an error to the log
        param:
            message [str] string to log
        """
        if self.level > 0:
            self._pprint('ERROR', message, calling=calling)
        self._add2log('errors', message, calling=calling)
        if exception is None:
            if self.no_exit:
                exception = PyNebError
            else:
                exception = SystemExit
        raise exception(message)
              
    def timer(self, message, quiet=False, calling=None):
        """
        method to add a timer to the log. Print the time elapsed since last call of timer (or instanciation
        of the object if first call).
        param:
            message [str] string to log
            quiet [boolean] if True, no message printed (but timer reset)
        """
        if calling is None:
            calling = self.calling
        delta_t = str((time() - self.__start))
        if not quiet:
            print('   {0}: {1} in {2}'.format(calling, message, delta_t))
        self.log_['timer'].append((calling, delta_t, message))
        if self._to_file:
            self._file_out.write("timer - {0} @ {1}: {2} \n".format(calling, delta_t, message))
        self.__start = time()

    def print_messages(self):
        """
        print all the messages.
        """
        to_print = '{0[0]:' + str(self.caller_max_size) + 's} at {0[1]}: {0[2]}' 
        for message in self.log_['messages']:
            print(to_print.format(message))
    
    def print_warnings(self):
        """
        print all the warnings.
        """
        to_print = '{0[0]:' + str(self.caller_max_size) + 's} at {0[1]}: {0[2]}' 
        for message in self.log_['warnings']:
            print(to_print.format(message))
    
    def print_errors(self):
        """
        print all the errors.
        """
        to_print = '{0[0]:' + str(self.caller_max_size) + 's} at {0[1]}: {0[2]}' 
        for message in self.log_['errors']:
            print(to_print.format(message))
    
    def print_debug(self):
        """
        print all the debug information.
        """
        to_print = '{0[0]:' + str(self.caller_max_size) + 's} at {0[1]}: {0[2]}' 
        for message in self.log_['debug']:
            print(to_print.format(message))
    
    def print_timer(self):
        """
        print the timer messages.
        """
        print(self.log_['timer'])
        
    def open_file(self, file_):
        """
        Open the file for the output and allow printing to it.
        param:
            file_ [str] file to save the log events.
        """
        if file_ is not None:
            self._file_out = open(file_, 'w')
            self._to_file = True
        else:
            self._to_file = False

    def close_file(self):
        """
        Close the output file and avoid printing to it.
        """
        if self._to_file:
            self._file_out.close()
            self._to_file = False

