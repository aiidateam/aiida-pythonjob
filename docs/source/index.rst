AiiDA PythonJob & PyFunction
============================

`PythonJob` and `pyfunction` enable users to run Python functions—either locally or remotely—with automatic serialization, provenance tracking, and workflow integration.

These tools are designed to simplify the experience for users not familiar with AiiDA's internal types, allowing them to write pure Python functions while still benefiting from AiiDA's powerful infrastructure.

**Key Features**

1. **Remote Execution with PythonJob**  
   Seamlessly run Python functions on a remote computer via `PythonJob`. A working directory is automatically created for the job, allowing full support for:
   
   - File upload and retrieval
   - Parent/remote folders
   - Integration with external codes (e.g. ASE + DFT engines)

2. **Local Execution with PyFunction**  
   Use `@pyfunction` to run pure Python functions locally. All inputs and outputs are automatically serialized/deserialized, enabling users to work with native Python types.  
   This decorator is ideal for functions that are:
   
   .. note::
      `pyfunction` does **not** create a dedicated working directory. It executes in the current local Python environment.  
      File I/O or relative path access should be avoided unless explicitly handled.

3. **User-Friendly Interface**  
   Designed for users unfamiliar with AiiDA’s internal `Data` classes. The decorators handle all conversion and provenance tracking behind the scenes.

4. **Workflow Management with Checkpoints**  
   Combine `pyfunction` and `pythonjob` in `WorkGraph` workflows to build robust, restartable workflows with full checkpoint support.

5. **Full Data Provenance**  
   All inputs, outputs, and intermediate results are stored in the AiiDA provenance graph, ensuring reproducibility and traceability.

**When to Use What**

+----------------+--------------------+---------------------+------------------+---------------------------------------------+
| Decorator      | Execution          | Input/Output Types  | AiiDA Provenance | Recommended Use Case                        |
+================+====================+=====================+==================+=============================================+
| ``@pyfunction``| Local              | Python native       | ✅ Yes           | Pure functions                              |
+----------------+--------------------+---------------------+------------------+---------------------------------------------+
| ``@pythonjob`` | Remote             | Python native       | ✅ Yes           | Remote jobs with file handling              |
+----------------+--------------------+---------------------+------------------+---------------------------------------------+

.. toctree::
   :maxdepth: 1
   :caption: Contents:
   :hidden:

   installation
   autogen/pythonjob
   autogen/pyfunction
   tutorial/index
