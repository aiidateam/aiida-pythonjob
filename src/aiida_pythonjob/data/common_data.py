import datetime

from aiida import orm
from aiida.orm import Data


class NoneData(orm.Data):
    """
    A Data node that explicitly represents a Python `None`.

    - Has no repository content and (by default) no attributes.
    - All instances have identical content hash, which is desirable here:
      "None" is a single value.
    """

    # Convenience aliases to mirror simple-* nodes (Int, Bool, etc.)
    @property
    def value(self):
        return None

    @property
    def obj(self):
        return None

    def __repr__(self) -> str:
        return "NoneData()"

    def __str__(self) -> str:
        return "NoneData()"


class DateTimeData(Data):
    """AiiDA node to store a datetime.datetime object."""

    def __init__(self, value: datetime.datetime, **kwargs):
        if not isinstance(value, datetime.datetime):
            raise TypeError(f"Expected datetime.datetime, got {type(value)}")
        super().__init__(**kwargs)
        # Store as ISO string for portability
        self.base.attributes.set("datetime", value.isoformat())

    @property
    def value(self) -> datetime.datetime:
        """Return the stored datetime as a datetime object."""
        return datetime.datetime.fromisoformat(self.base.attributes.get("datetime"))

    def __str__(self):
        return str(self.value)
