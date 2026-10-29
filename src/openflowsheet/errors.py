"""Errors that explain themselves.

Every error raised while setting up or solving a flowsheet says *where* it happened
(unit or stream id), *what* went wrong, and *what to try* (``hint``).
"""

from __future__ import annotations


class SimulationError(Exception):
    def __init__(self, message: str, *, where: str | None = None, hint: str | None = None):
        super().__init__(message)
        self.message = message
        self.where = where
        self.hint = hint

    def located(self, where: str) -> SimulationError:
        """Return the same error attributed to ``where`` unless already attributed."""
        if self.where is None:
            self.where = where
        return self

    def to_dict(self) -> dict:
        return {"message": self.message, "where": self.where, "hint": self.hint}

    def __str__(self) -> str:
        text = f"{self.where}: {self.message}" if self.where else self.message
        return f"{text}\n  hint: {self.hint}" if self.hint else text
