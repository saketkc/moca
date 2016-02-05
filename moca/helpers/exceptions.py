class MocaException(Exception):
    """Base class for MoCA Exceptions"""

    def __init__(self, msg):
        self.value = msg

    def __str__(self):
        """string representation of MoCA Exception

        Returns
        -------
        mocastr: string representation

        """
        mocastr = repr(self.value)
        return mocastr
