# AiiDA-PythonJob
[![PyPI version](https://badge.fury.io/py/aiida-pythonjob.svg)](https://badge.fury.io/py/aiida-pythonjob)
[![Unit test](https://github.com/aiidateam/aiida-pythonjob/actions/workflows/ci.yml/badge.svg)](https://github.com/aiidateam/aiida-pythonjob/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/aiidateam/aiida-pythonjob/branch/main/graph/badge.svg)](https://codecov.io/gh/aiidateam/aiida-pythonjob)
[![Docs status](https://readthedocs.org/projects/aiida-pythonjob/badge)](http://aiida-pythonjob.readthedocs.io/)

`PythonJob` allows users to run Python functions on a remote computer. It is designed to enable users from non-AiiDA communities to run their Python functions remotely and construct workflows with checkpoints, maintaining all data provenance. For instance, users can use ASE's calculator to run a DFT calculation on a remote computer directly.

## Key Features

1. **Remote Execution**: Seamlessly run Python functions on a remote computer.
2. **User-Friendly**: Designed for users who are not familiar with AiiDA, simplifying the process of remote execution.
3. **Workflow Management**: Construct workflows using WorkGraph with checkpoints, ensuring that intermediate states and results are preserved.
4. **Data Provenance**: Maintain comprehensive data provenance, tracking the full history and transformations of data.



## Installation

```console
    pip install aiida-pythonjob
```

To install the latest version from source, first clone the repository and then install using `pip`:

```console
git clone https://github.com/aiidateam/aiida-pythonjob
cd aiida-pythonjob
pip install -e .
```


## Documentation
Explore the comprehensive [documentation](https://aiida-pythonjob.readthedocs.io/en/latest/) to discover all the features and capabilities of AiiDA PythonJob.


## License
[MIT](http://opensource.org/licenses/MIT)
