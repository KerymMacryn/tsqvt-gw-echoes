# Contributing to TSQVT GW Echoes

Thank you for your interest in contributing to this project! This document provides guidelines for contributing.

## Code of Conduct

This project adheres to a code of conduct of respect, inclusivity, and scientific rigor. By participating, you are expected to uphold these values.

## How to Contribute

### Reporting Bugs

If you find a bug:
1. Check if the issue already exists in the [Issues](https://github.com/yourusername/tsqvt-gw-echoes/issues)
2. If not, create a new issue with:
   - Clear title and description
   - Steps to reproduce
   - Expected vs. actual behavior
   - System information (OS, Python version, package versions)
   - Minimal code example demonstrating the bug

### Suggesting Enhancements

Enhancement suggestions are welcome! Please:
1. Check existing issues/pull requests first
2. Provide clear use case and rationale
3. If possible, outline implementation approach

### Pull Requests

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/amazing-feature`
3. Make your changes following the style guide below
4. Add tests for new functionality
5. Ensure all tests pass: `pytest tests/`
6. Update documentation as needed
7. Commit with clear messages: `git commit -m 'Add amazing feature'`
8. Push to your fork: `git push origin feature/amazing-feature`
9. Open a Pull Request with:
   - Clear description of changes
   - Reference to related issues
   - Summary of testing performed

## Development Setup

```bash
# Clone your fork
git clone https://github.com/yourusername/tsqvt-gw-echoes.git
cd tsqvt-gw-echoes

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install in development mode
pip install -e ".[dev]"

# Run tests
pytest tests/ -v
```

## Style Guide

### Python Code

- Follow [PEP 8](https://pep8.org/)
- Use [Black](https://black.readthedocs.io/) for formatting: `black src/ tests/`
- Maximum line length: 100 characters
- Use type hints where appropriate
- Write docstrings for all public functions (Google style)

Example:
```python
def compute_echo_delay(
    M_total: float,
    M_star: float,
    xi: float,
    n: int
) -> float:
    """
    Compute the time delay for the n-th gravitational wave echo.
    
    Args:
        M_total: Total black hole mass in solar masses
        M_star: Spectral gap scale in GeV
        xi: Spectral coupling parameter
        n: Echo number (n=1 for first echo)
        
    Returns:
        Echo delay in seconds
        
    Raises:
        ValueError: If n < 1 or M_total <= 0
    """
    ...
```

### Documentation

- Update README.md for user-facing changes
- Update docs/ for theoretical/technical details
- Add docstrings to all new functions/classes
- Include usage examples where appropriate

### Testing

- Write tests for all new features
- Maintain >80% code coverage
- Use pytest fixtures for common setups
- Test edge cases and error conditions

Example:
```python
def test_echo_delay_positive():
    """Test that echo delays are always positive."""
    model = TSQVTEchoModel(M_total=65.0)
    for n in range(1, 10):
        assert model.echo_delay(n) > 0
```

## Commit Messages

- Use present tense: "Add feature" not "Added feature"
- Use imperative mood: "Move cursor to..." not "Moves cursor to..."
- First line: brief summary (50 chars or less)
- Blank line, then detailed explanation if needed
- Reference issues: "Fixes #123"

Example:
```
Add frequency modulation to echo model

- Implement beta parameter for frequency shifts
- Add tests for frequency_shift() method
- Update documentation with new parameter

Fixes #42
```

## Areas for Contribution

We welcome contributions in:

### Code
- Additional echo models for comparison
- Improved PSD estimation methods
- Bayesian parameter estimation (MCMC/nested sampling)
- GPU acceleration for matched filtering
- Advanced glitch removal techniques

### Documentation
- Tutorial notebooks for specific use cases
- Theoretical background explanations
- API reference improvements
- Video tutorials

### Testing
- Additional unit tests
- Integration tests
- Performance benchmarks
- Validation against known results

### Analysis
- Application to new GW events
- Statistical methods improvements
- Multi-detector coherent analysis
- Null hypothesis testing

## Questions?

Feel free to:
- Open an issue for questions
- Email the maintainer: mhamed34@alumno.uned.es
- Join discussions in existing issues/PRs

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

Thank you for helping improve this project!
