from aiida.orm import CalcJobNode


class PythonJobNode(CalcJobNode):
    """ORM class for the PythonJob plugin."""

    FUNCTION_DATA_KEY = "function_data"

    def _setup_db_record(self):
        super()._setup_db_record()
        self.set_function_data(self.inputs.function_data)

    def set_function_data(self, function_data: dict) -> None:
        """Set the function data dictionary."""
        self.set_attribute(self.FUNCTION_DATA_KEY, function_data)
