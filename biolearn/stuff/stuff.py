from .._utils import fill_doc

##############################################################################
# Class for stuff ############################################################
##############################################################################
@fill_doc
class Stuff():
    """Implement stuff

    Parameters
    ----------
    number : int
        Integer because why not

    Notes
    -----
    So much extra information on stuff
    """

    def __init__(
        self,
        number,
    ):
        self.number = number

    def fit(self, letter):
        """do that stuff

        Parameters
        ----------
        letter : str
            An ascii letter because why not

        """
        return self
