"""
Interface for the handler class
"""
from abc import ABC, abstractmethod
from typing import Optional


class Handler(ABC):
    """
    The Handler interface declares a method for building the chain of handlers.
    It also declares a method for executing a request.
    """
    @abstractmethod
    def set_next(self, handler: any) -> any:
        pass

    @abstractmethod
    def handle(self, request) -> Optional[str]:
        pass
