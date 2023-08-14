"""
Test drop comments function
"""
from minushalf.utils.drop_comments import drop_comments


def test_drop_comments_first():
    """
    Test drop comments
    """
    text = ["### ssss", "#ffff", "aleatory text"]
    parsed_text = drop_comments(text)

    assert parsed_text[0] == "aleatory text"


def test_drop_comments_second():
    """
    Test drop comments
    """
    text = [
        "###comentario 1", "#comentario 2", "aleatory text## sjsjsjs",
        "#sdjkdij"
    ]
    parsed_text = drop_comments(text)

    assert parsed_text[0] == "aleatory text"


def test_drop_comments_third():
    """
    Test drop comments
    """
    text = [
        "###comentario 1", "#comentario 2", "aleatory text## sjsjsjs",
        "other text", "#sdjkdij"
    ]
    parsed_text = drop_comments(text)

    assert parsed_text[0] == "aleatory text"
    assert parsed_text[1] == "other text"


def test_drop_comments_fourth():
    """
    Test drop comments
    """
    text = [
        "###comentario 1", "#comentario 2", "aleatory text## sjsjsjs",
        "other text", " ##physics", "#sdjkdij"
    ]
    parsed_text = drop_comments(text)

    assert parsed_text[0] == "aleatory text"
    assert parsed_text[1] == "other text"
    assert len(parsed_text) == 2


def test_drop_comments_fiveth():
    """
    Test drop comments
    """
    text = [
        "###comentario 1", "#comentario 2", "aleatory text## sjsjsjs",
        "other text", " ##physics", "final test ##ddd", "#sdjkdij"
    ]
    parsed_text = drop_comments(text)

    assert parsed_text[0] == "aleatory text"
    assert parsed_text[1] == "other text"
    assert parsed_text[2] == "final test "
    assert len(parsed_text) == 3
