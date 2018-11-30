## Preprocessing program for LumosVar

`gvm` is a tool that replaces the old LumosVar preprocess program.

### Dependencies

Make sure the following are installed on your system:

* `autotools` (for building)
* `htslib` (>=1.4.1)
* `libyaml`
* `gsl` (>=2.4)
* `gengetopt` (for building)

The following dependencies are bundled as submodules:

* `uthash`

If you're on the dback cluster, running 

```
$ source ./setup.sh
```

will load the modules that satisfy these dependencies. If not, you have to
install them yourself.

### Cloning

Make sure to use `--recursive` when cloning.

### Building

#### Scripted

Run the following:

```
$ ./build.sh --init
```

And everything should get built. The `gvm` executable should be placed
in the `src` directory. If dependencies were not met, the configure step will
fail and tell you to install something. You can omit the `--init` on subsequent
runs of `build.sh`

You can supply `--debug` to build a debug version and `--clean` to remove
all generated files.

#### Manual

Run

```
$ autoreconf --install
$ ./configure
$ make
```

### Running

#### Standard

Run

```bash
$ src/gvm --help
```

to see the options available. Generally you're going to want to do something
like this:

```bash
$ src/gvm -c <path to your config.yaml> -C <chromosome to run on>
```

#### Normal Metrics

To calculate normal metrics, you have to set the flag for it
(`-N`). However, you might also want to turn pos and exon file
generation off and that is done with their respective flags as well
(`-P` and `-E`). Note that all of these are toggles but `-P` and `-E`
are on by default and `-N` is off.

If the YAML configuration is malformed for the task at hand, you might
get a segmentation fault because I haven't had the chance to verify
that I wrote the checks for every condition necessary to calculate
normal metrics.

There is also a helper script that provides arguments with which `-N`
must co-occur. The additional arguments are printed to standard
output.

```bash
$ python scripts/nmconf.py conf.yaml 6
--ploidystr=222
```

So, a complete invocation of gvm for normal metrics calculation would
be, given a `conf.yaml` that is made for normal metrics:

```bash
CONF=conf.yaml
CHR=6
NMOPTS=$(python scripts/nmconf.py $CONF $CHR)
gvm -c $CONF -C $CHR -P -E -N $NMOPTS
```

Note that there is nothing stopping you from running normal metrics as
well as pos/exon file generation (supply `-N` but not `-P` nor `-E`)
if you can provide the configuration that satisfies both tasks. This,
however, has not been tested.
