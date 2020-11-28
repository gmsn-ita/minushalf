"""
Function to find elements
in a array of strings
"""


def find_element(array: list, element: str, initial_index: str = None, final_index: str = None) -> int:
    """
        Find string elements in list
        comparing them without whitespace
        and new lines
    """

    if not initial_index or initial_index < 0:
        initial_index = 0

    if not final_index or final_index > len(array):
        final_index = len(array)

    for index in range(initial_index, final_index):

        if array[index].replace(" ", "").rstrip('\n') == element.replace(" ", "").rstrip('\n'):
            return index

    raise Exception("Value not found")
