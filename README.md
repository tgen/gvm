## Preprocessing program for LumosVar

`gvm` is a tool that replaces the old LumosVar preprocess program.

### Dependencies

* `htslib` (preferably a recent version)
* `libyaml`
* `uthash` (for building, bundled as submodule)
* `gengetopts` (for building)

### Cloning

Make sure to use `--recursive` when cloning.

### Building

If you're on the dback cluster, all you should need to do is run:

```
$ ./build.sh --init
```

And everything should get built. The `gvm` executable should be placed
in the `src` directory.

### Running

Run

```
$ src/gvm --help
```

to see the options available. Generally you're going to want to do something
like this:

```
$ src/gvm -c <path to your config.yaml> -C <chromosome to run on>
```

#### The slurm script

There is also a slurm script to automate running gvm on every chromosome
(1,2,...,22,X,Y) but you should change some of the variables in there to avoid
messing with my scratch directory. They are declared near the top.

Namely:
* `$SC`
* `$CONF`

The script will substitute any mention of `$OUTFILE` in your config with a
newly created directory in your `$SC` (scratch folder if you want). This
directory will share the name of the job that is created, so it will look like
`110205` or something. It will then copy the substituted config to that
directory and run gvm there.

If you don't want all of this then just run it like above.
