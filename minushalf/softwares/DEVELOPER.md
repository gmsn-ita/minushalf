
# Software Interface Implementation Guide

This guide outlines the steps for developers to implement a new software interface using the `SoftwaresAbstractFactory` and how to organize the implementation details within specific folders.

## SoftwaresAbstractFactory Interface

The `SoftwaresAbstractFactory` is a factory interface that provides methods for creating various software components. Developers should follow this interface when implementing a new software interface. The factory methods should return instances of specific software components.


## Implementation Steps

To implement a new software interface, follow these steps:

1. **Create a New Folder**: Create a new folder with the name of the software. For example, if you are implementing an interface for "VASP," create a folder named `vasp`.

2. **Implement Software Factory**: Inside the software folder, implement a class that adheres to the `SoftwaresAbstractFactory` interface. This class should provide concrete implementations for the factory methods.

3. **Organize Implementation Details**: Within the software folder (`vasp` in this example), organize the implementation details, code, and resources related to the software interface.

   ```
   vasp/
   ├── eigenval.py
   ├── proar.py
   ├── potcar.py
   ├── README.md
   └── ... other implementation files and resources
   ```

4. **Integration**: Integrate the implemented software interface and its factory class into your application's main codebase. To do that you should add a new entry on the `Softwares` Enum and add it to the function `get_software_factory`.