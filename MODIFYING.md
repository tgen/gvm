## gvm source code layout

Important files:

```
project-root/
    src/
        gvm.c               <- main source code
        variant_table.h     <- data structure definitions/explanation
        config_keys.h       <- configuration options definitions
        cigar.c             <- alignment processor
        gengetopt/
            gvm.ggo         <- commandline arguments
```

To learn how to accomplish certain tasks, read the comment at the top
of `gvm.c` and follow the instructions.
