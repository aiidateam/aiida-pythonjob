from __future__ import annotations


def generate_script_py(
    pickled_function: bytes | None, source_code: str | None, function_name: str = "user_function"
) -> str:
    """
    Generate the script.py content as a single string with robust exception handling.

    :param pickled_function: Serialized function bytes if running in pickled mode, else None.
    :param source_code: Raw Python source code if running in source-code mode, else None.
    :param function_name: The name of the function to call when running source code mode.
    :return: A string representing the entire content of script.py.
    """
    # We build a list of lines, then join them with '\n' at the end
    script_lines = [
        "import sys",
        "import json",
        "import traceback",
        "",
        "def write_error_file(error_type, exc, traceback_str):",
        "    # Write an error file to disk so the parser can detect the error",
        "    error_data = {",
        "        'error_type': error_type,",
        "        'exception_message': str(exc),",
        "        'traceback': traceback_str,",
        "    }",
        "    with open('_error.json', 'w') as f:",
        "        json.dump(error_data, f, indent=2)",
        "",
        "def main():",
        "    # 1) Attempt to import cloudpickle",
        "    try:",
        "        import cloudpickle as pickle",
        "    except ImportError as e:",
        "        write_error_file('IMPORT_CLOUDPICKLE_FAILED', e, traceback.format_exc())",
        "        sys.exit(1)",
        "",
        "    # 2) Attempt to unpickle the inputs",
        "    try:",
        "        with open('inputs.pickle', 'rb') as handle:",
        "            inputs = pickle.load(handle)",
        "    except Exception as e:",
        "        write_error_file('UNPICKLE_INPUTS_FAILED', e, traceback.format_exc())",
        "        sys.exit(1)",
        "",
    ]

    if pickled_function:
        # Mode 1: pickled function
        script_lines += [
            "    # 3) Attempt to unpickle the function",
            "    try:",
            "        with open('function.pkl', 'rb') as f:",
            "            user_function = pickle.load(f)",
            "    except Exception as e:",
            "        write_error_file('UNPICKLE_FUNCTION_FAILED', e, traceback.format_exc())",
            "        sys.exit(1)",
            "",
            "    # 4) Attempt to run the function",
            "    try:",
            "        result = user_function(**inputs)",
            "    except Exception as e:",
            "        write_error_file('FUNCTION_EXECUTION_FAILED', e, traceback.format_exc())",
            "        sys.exit(1)",
        ]
    elif source_code:
        # Mode 2: raw source code
        # Indent each line of source_code by 4 spaces to keep correct indentation
        source_lines = [f"    {line}" for line in source_code.split("\n")]
        script_lines += [
            "    # 3) Define the function from raw source code",
            *source_lines,
            "",
            "    # 4) Attempt to run the function",
            "    try:",
            f"        result = {function_name}(**inputs)",
            "    except Exception as e:",
            "        write_error_file('FUNCTION_EXECUTION_FAILED', e, traceback.format_exc())",
            "        sys.exit(1)",
        ]
    else:
        raise ValueError("You must provide exactly one of 'source_code' or 'pickled_function'.")

    # 5) Attempt to pickle (save) the result
    script_lines += [
        "",
        "    # 5) Attempt to pickle the result",
        "    try:",
        "        with open('results.pickle', 'wb') as handle:",
        "            pickle.dump(result, handle)",
        "    except Exception as e:",
        "        write_error_file('PICKLE_RESULTS_FAILED', e, traceback.format_exc())",
        "        sys.exit(1)",
        "",
        "    # If we've made it this far, everything succeeded. Write an empty _error.json",
        "    # so the parser can always read _error.json (if it's empty, no error).",
        "    with open('_error.json', 'w') as f:",
        "        json.dump({}, f, indent=2)",
        "",
        "if __name__ == '__main__':",
        "    main()",
        "",
    ]

    # Join lines with newline
    script_content = "\n".join(script_lines)
    return script_content
