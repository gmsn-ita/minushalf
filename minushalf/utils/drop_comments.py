"""
Function to check if the line
a commentarie or not
"""
import re


def drop_comments(lines: str) -> list:
    """function to check if a line
        starts with some character.
        Here # for comment
    """

    comment_pattern = r"#.*"
    lines_without_comments = [
        re.sub(comment_pattern, "", line) for line in lines
    ]

    return list(filter(lambda x: x != "\n" and x != "",
                       lines_without_comments))
