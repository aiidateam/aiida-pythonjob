from aiida import orm


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
