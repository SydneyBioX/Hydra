# 📖 Contributing 

[![PyPI version](https://img.shields.io/pypi/v/Hydra-tools?color=orange)](https://pypi.org/project/Hydra-tools/)
[![Docs](https://img.shields.io/badge/docs-passing-brightgreen)](https://sydneybiox.github.io/Hydra/)
![Python](https://img.shields.io/badge/python-%3E%3D3.8-blue)
![R](https://img.shields.io/badge/R-%3E%3D4.0-blueviolet)
[![License](https://img.shields.io/badge/license-MIT-green)](https://github.com/SydneyBioX/Hydra?tab=MIT-1-ov-file#readme)


Thank you for using <code><span style="color: red;">Hydra</span></code> - We welcome bug fixes, documentation improvements, and other new features.

## Ways to contribute
- Report bugs in the [issue tracker](https://github.com/SydneyBioX/Hydra/issues)
- Improve the documentation by fixing typos, clarifying steps, or adding examples

## Reporting an issue
Before opening a new issue, search existing issues to avoid duplicates. If you are considering opening an issue, please include:

- Hydra version, OS, Python and R versions, GPU and other relevant info
- Exact code snippet to reproduce the problem
- Full error message or traceback

## Development setup
1. Fork the repo on GitHub and clone your fork:
   ```bash
   git clone https://github.com/SydneyBioX/Hydra.git
   cd Hydra
   ```
2. Create a isolated environment (Ex: Mamba)
3. Install Hydra in editable mode with development extras:
   ```bash
   pip install -e ".[dev]"
   ```
4. Verify your setup:
   ```bash
   pytest
   ```

## Making changes
- Work on a feature branch:
   ```bash
   git checkout -b feature/short-description
   ```
- Please write clear messages before committing
- Update documentation if the user-facing interface changes

## License
By contributing, you agree your contributions are licensed under the project’s <a href="https://github.com/SydneyBioX/Hydra?tab=MIT-1-ov-file#readme">MIT License</a>.

<br>

---
<p style="text-align: left; font-size: 15px">
  Documentation by <a href="http://manojmw.github.io" target="_blank">Manoj M Wagle</a>
</p>
