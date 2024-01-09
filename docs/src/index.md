````@eval
using Markdown
Markdown.parse("""
$(read("../../README.md",String))
""")
````


```@docs
AMGSolver
RLXSolver
AMGPrecon
RLXPrecon
error_state
blocksize_instantiated
```


