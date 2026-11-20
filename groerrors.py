import logging
logger_GROerror = logging.getLogger('pymemdyn.groerrors')

class GromacsError(BaseException):
    pass


class IOGromacsError(GromacsError):
    """
    Exception raised with "File input/output error" message
    """
    def __init__(self, command, explain):
        self.command = command
        self.explain = explain

    def __str__(self):
        return "%s failed with '%s'" % (self.command, self.explain)


class GromacsMessages(object):
    """
    Load an error message and split it along as many properties as
    possible
    """

    # Map the messages with the errors
    e = {"File input/output error": IOGromacsError,
         "Can not open file": IOGromacsError,
         # FIXME: Clearly the following is not IOError
         "srun: error: Unable to create job step": IOGromacsError,
         "Fatal error:": IOGromacsError,         
         }

    def __init__(self, gro_err="", command="", *args, **kwargs):
        """
        Pass the command and the output of that command in "command" and
        "gro_err" kwargs. The check for "error" property
        """
        self.logger_groerrors = logging.getLogger('pymemdyn.groerrors.GromacsMessages')
        self.error = False
        self.command = command
        self.gro_err = gro_err.decode('utf-8').split("\n")

        self.check()

    def check(self):
        """
        Check if the GROMACS error message has any of the known error
        messages. Set the self.error to the value of the error
        """
        for line in self.gro_err:
            for error in self.e.keys():
                if line.startswith(error):
                    if error == "Fatal error:":
                        self.logger_groerrors.error("\n\nGROMACS FATAL ERROR:\n")
                        try:
                            j = self.gro_err.index('of the License, or (at your option) any later version.')
                            self.logger_groerrors.error('\n'.join(self.gro_err[j+1:]))
                        except ValueError:
                            # Print the last 40 lines if the license string is missing
                            self.logger_groerrors.error('\n'.join(self.gro_err[-40:]))
                        self.logger_groerrors.error("\n\n GROMACS FATAL ERROR\n")

                    raise self.e[error](self.command, error)