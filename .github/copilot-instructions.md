# C++ Coding Conventions

Please follow these rules when generating or modifying code:

  **No `auto` keyword**: Always use explicit types. Do not use type deduction.
  **No Inheritance**: Avoid class inheritance. Use composition or plain structs/functions.
  **Style**: Follow **Google C++ Style** conventions.
  **“Plain code” preference** (avoid unnecessary built-in abstractions)
- Avoid built in abstractions unless it is a routine one like <stdio> and <vector> or it is absolutely necessary(like random number generator and timer...) or the same abstraction appears many times accross the project(like smart pointer).
 
