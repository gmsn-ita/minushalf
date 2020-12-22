"""
Function to check if the line
a commentarie or not
"""
import re


def drop_comments(lines: list) -> list:
    """
    Function to remove comments from lines in a file

        Args:
            lines(list): list of file lines
        Returns:
            lines_without_comments (list): lines of the file without comments
    """

    comment_pattern = r"#.*"
    lines_without_comments = [
        re.sub(comment_pattern, "", line) for line in lines
    ]

    return list(
        filter(lambda x: x != "\n" and x.strip() != "",
               lines_without_comments))
