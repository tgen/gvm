## Preprocessing program for LumosVar

`gvm` is a tool that replaces the old LumosVar preprocess program.

### Dependencies

* `htslib` (preferably a recent version)
* `libyaml`
* `uthash` (bundled as submodule)

### Cloning

Make sure to use `--recursive` when cloning.

### Building

I use gcc 5.1. Older versions create issues, newer ones probably won''t.

To build, run:

```
$ ./setup.sh
$ make
```

and hope nothing goes wrong. This will place an executable in the `prod` directory.

### Running

To run:

```
$ ./prod/gvm conf.yaml <chromosome>
```

This will generate the pos and exon files and name them as per the config file.

Use the same config you use for LumosVar.

### Warning

This has only been tested on PNAP.

