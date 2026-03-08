"""Modular validation checks for the ITm seed leaderboard."""

from dataclasses import dataclass, field


@dataclass
class CheckResult:
    """Result of a single validation check category."""
    status: str = "pass"  # "pass", "warn", "fail"
    messages: list[str] = field(default_factory=list)

    def fail(self, msg: str):
        self.status = "fail"
        self.messages.append(msg)

    def warn(self, msg: str):
        if self.status != "fail":
            self.status = "warn"
        self.messages.append(msg)

    def info(self, msg: str):
        self.messages.append(msg)

    @property
    def passed(self) -> bool:
        return self.status != "fail"

    def to_dict(self) -> dict:
        return {"status": self.status, "messages": self.messages}
