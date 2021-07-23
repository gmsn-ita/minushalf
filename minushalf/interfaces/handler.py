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
        """
        Set next handler
        """
        pass

    @abstractmethod
    def handle(self, request) -> Optional[any]:
        """
        Call action and pass to the next handler
        """
        pass
    
    @abstractmethod
    def action(self,request:any) -> any:
        """
        Performs the check or the process
        """
        pass
