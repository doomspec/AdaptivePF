class Task:
    """
    The base class of tasks. Task should include information for a TaskRunner to process
    By run(), a result should be returned
    """

    def __init__(self):
        self.index_of_in = -1
        self.series_id = -1
        return

    def run(self):
        return

    def __hash__(self):
        return self.__str__().__hash__()
