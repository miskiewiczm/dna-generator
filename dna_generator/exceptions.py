"""
Custom exceptions used by the DNA sequence generator.

Defines rich exceptions used across the package to signal specific error types.
"""


class GeneratorError(Exception):
    """Base exception for all DNA generator errors."""
    
    def __init__(self, message: str, **kwargs):
        """
        Initialize exception.
        
        Args:
            message: Error message
            **kwargs: Context parameters
        """
        super().__init__(message)
        self.message = message
        self.context = kwargs
    
    def __str__(self) -> str:
        """String representation with context."""
        if self.context:
            context_str = ", ".join(f"{k}={v}" for k, v in self.context.items())
            return f"{self.message} [{context_str}]"
        return self.message


class ConfigurationError(GeneratorError):
    """Configuration error."""
    pass


class ValidationError(GeneratorError):
    """Sequence validation error."""
    
    def __init__(self, message: str, sequence: str = None, **kwargs):
        """
        Initialize.
        
        Args:
            message: Error message
            sequence: Sequence that failed validation
            **kwargs: Extra parameters
        """
        super().__init__(message, sequence=sequence, **kwargs)
        self.sequence = sequence


class BacktrackingError(GeneratorError):
    """Backtracking algorithm error."""
    
    def __init__(self, message: str, attempts: int = None, position: int = None, **kwargs):
        """
        Initialize.
        
        Args:
            message: Error message
            attempts: Number of attempts
            position: Position in sequence
            **kwargs: Extra parameters
        """
        super().__init__(message, attempts=attempts, position=position, **kwargs)
        self.attempts = attempts
        self.position = position


class Primer3Error(GeneratorError):
    """primer3-related error."""
    
    def __init__(self, message: str, operation: str = None, **kwargs):
        """
        Initialize.
        
        Args:
            message: Error message
            operation: primer3 operation name
            **kwargs: Extra parameters
        """
        super().__init__(message, operation=operation, **kwargs)
        self.operation = operation


class InputError(GeneratorError):
    """Invalid input error."""
    
    def __init__(self, message: str, parameter: str = None, value: any = None, **kwargs):
        """
        Initialize.
        
        Args:
            message: Error message
            parameter: Parameter name
            value: Invalid value
            **kwargs: Extra parameters
        """
        super().__init__(message, parameter=parameter, value=value, **kwargs)
        self.parameter = parameter
        self.value = value
