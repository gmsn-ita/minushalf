"""
Base class to all Handlers
"""
from minushalf.interfaces import Handler


class BaseHandler(Handler):
    """
    The default chaining behavior is defined inside a base handler
    class.
    """

    _next_handler: Handler = None

    def set_next(self, handler: Handler) -> Handler:
        self._next_handler = handler
        # Returning a handler from here will let us link handlers in a useful way
        return handler

    def handle(self, request: any) -> str:
        request = self.action(request)
        if self._next_handler:
            return self._next_handler.handle(request)

        return request
    
    def action(self, request: any) -> any:
        return request
    
