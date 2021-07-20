"""
Test base handler class
"""
from minushalf.handlers import BaseHandler
from pytest_mock import MockerFixture


def text_set_next():
    """
    Test function that sets the next handler
    """
    base_handler = BaseHandler()
    next_handler = BaseHandler()

    res = base_handler.set_next(next_handler)

    assert base_handler._next_handler == next_handler
    assert res == next_handler


def test_handle_when_next_handler_not_defined():
    """
    Test function that activates the next handler when
    next handler isn't defined
    """
    base_handler = BaseHandler()

    assert base_handler.handle({}) == ''


def test_handle_when_next_handler_is_defined(mocker: MockerFixture):
    """
    Test function that activates the next handler when
    next handler is defined
    """
    base_handler = BaseHandler()
    next_handler = BaseHandler()

    ## mock and spy next handler
    mocker.patch.object(next_handler,
                        'handle',
                        return_value='Test',
                        autospec=True)
    spy = mocker.spy(next_handler, 'handle')

    base_handler.set_next(next_handler)

    assert base_handler.handle({}) == 'Test'
    assert spy.assert_called_once_with({})
