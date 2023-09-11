# Developer Guide to Minushalf

**Welcome to the Developer Guide for MinusHalf!**

In this markdown guide, we'll delve into everything you need to know about Minushalf. Whether you're an experienced developer or just starting out, this guide will walk you through the features, capabilities, and best practices for working with MinusHalf. 

## Installing pyenv and Switching to Python 3.6.15

[`pyenv`](https://github.com/pyenv/pyenv) is a popular tool for managing multiple versions of Python on your system. Here's how you can install `pyenv` and switch to Python 3.6.15 for your MinusHalf development:

1. **Install pyenv**: Open a terminal window and execute the following commands:

   ```bash
   # Install pyenv using the official installer script
   curl -L https://github.com/pyenv/pyenv-installer/raw/master/bin/pyenv-installer | bash

   # Add pyenv initialization to your shell profile (e.g., .bashrc or .zshrc)
   echo 'export PATH="$HOME/.pyenv/bin:$PATH"' >> ~/.bashrc
   echo 'eval "$(pyenv init --path)"' >> ~/.bashrc
   echo 'eval "$(pyenv virtualenv-init -)"' >> ~/.bashrc
   source ~/.bashrc
   ```

   This will install `pyenv` and set up the necessary environment variables in your shell profile.

2. **Install Python 3.6.15**: Use `pyenv` to install Python 3.6.15 by executing the following command:

   ```bash
   pyenv install 3.6.15
   ```

   This will download and install Python 3.6.15.

3. **Set Global Python Version**: Set the global Python version to 3.6.15 with the following command:

   ```bash
   pyenv global 3.6.15
   ```

   This ensures that Python 3.6.15 is used as the default version across your system.

4. **Verify Python Version**: Confirm that Python 3.6.15 is now the active version by running:

   ```bash
   python --version
   ```

   You should see the output indicating Python 3.6.15.

Now you have `pyenv` installed and Python 3.6.15 set as your active version. You can proceed with your MinusHalf development using this specific Python version. Remember to refer to the MinusHalf documentation and the earlier sections of this guide for a complete development setup.


## Installation of Poetry

Poetry is a powerful dependency management and packaging tool for Python. Follow these steps to install Poetry and set up your MinusHalf development environment:

1. **Install Poetry**: Open a terminal window and execute the following command:

   ```bash
   curl -sSL https://install.python-poetry.org | python3 - --version 1.3.2
   ```

   This will download and install Poetry on your system.

2. **Verify Installation**: After the installation is complete, verify that Poetry was installed successfully by running:

   ```bash
   poetry --version
   ```

   You should see the version number of Poetry displayed.

3. **Clone MinusHalf Repository**: If you haven't already, clone the MinusHalf repository from your preferred source code repository using Git:

   ```bash
   git clone https://github.com/gmsn-ita/minushalf
   ```

4. **Navigate to Project Directory**: Change into the project directory:

   ```bash
   cd minushalf
   ```

5. **Activate Virtual Environments**: While Poetry creates a virtual environment automatically, you can manually activate it for your current session:

   ```bash
   make virtual_env
   ```

   This ensures you're working within the isolated environment created by Poetry.

6. **Install Dependencies**: Use Poetry to install the project dependencies:

   ```bash
   poetry install
   ```

   Poetry will read the `pyproject.toml` file in your project directory and set up the necessary dependencies and virtual environment.

6. **Run tests**: To make sure you setup the env correctly, run the command:

   ```bash
   pytest
   ```

You're now all set to start developing with Minushalf using the power of Poetry! Remember to consult the MinusHalf documentation and this guide for further details on using MinusHalf effectively.

Happy coding!
