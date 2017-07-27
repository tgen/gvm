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

```
$ ./setup.sh
$ make
```

and hope nothing goes wrong.

### Running

Use your LumosVar `conf.yaml` here too. The executable will be placed in the `prod`
directory.

```
$ prod/gvm conf.yaml <chromosome>
```

This will generate the pos and exon files and name them as per the config file.

### Warning

This has only been tested on PNAP.

