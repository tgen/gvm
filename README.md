## Preprocessing program for LumosVar

`gvm` is a tool that replaces the old LumosVar preprocess program.

### Dependencies

Make sure the following are installed on your system:

* `autotools` (for building)
* `htslib` (>=1.4.1)
* `libyaml`
* `gengetopts` (for building)

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

Run

```
$ src/gvm --help
```

to see the options available. Generally you're going to want to do something
like this:

```
$ src/gvm -c <path to your config.yaml> -C <chromosome to run on>
```
