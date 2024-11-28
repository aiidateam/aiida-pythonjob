from typing import Optional

from aiida.common.exceptions import NotExistent
from aiida.orm import Computer, InstalledCode, load_code, load_computer


def get_or_create_code(
    label: str = "python3",
    computer: Optional[str | Computer] = "localhost",
    filepath_executable: Optional[str] = None,
    prepend_text: str = "",
):
    """Try to load code, create if not exit."""

    try:
        return load_code(f"{label}@{computer}")
    except NotExistent:
        description = f"Code on computer: {computer}"
        computer = load_computer(computer)
        filepath_executable = filepath_executable or label
        code = InstalledCode(
            computer=computer,
            label=label,
            description=description,
            filepath_executable=filepath_executable,
            default_calc_job_plugin="pythonjob.pythonjob",
            prepend_text=prepend_text,
        )

        code.store()
        return code
