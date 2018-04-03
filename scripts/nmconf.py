import sys, yaml

def die(msg):
  sys.stderr.write("%s: " % sys.argv[0])
  sys.stderr.write(msg)
  sys.stderr.write("\n")
  exit(1)

def chr_lists(autosomes, sexs):
  a_start, a_end = map(int, autosomes.split(":"))

  sex_list = sexs.split(",")

  autosomes_list = list(range(a_start, a_end + 1))

  return list(map(str, autosomes_list)), sex_list
  

if __name__ == "__main__":
  if len(sys.argv) < 2:
    die("missing positional argument: config")

  if len(sys.argv) < 3:
    die("missing positional argument: chromosome")

  with open(sys.argv[1]) as f:
    config = yaml.load(f)

  chro = sys.argv[2]

  autosomes, sexs = chr_lists(config["autosomes"], config["sexChr"])
  chrs = autosomes + sexs

  sexlist = config["sexList"].split(",")
  nsample = len(sexlist)

  mlist = list(map(int, config["M"].split(",")))
  flist = list(map(int, config["F"].split(",")))

  try:
    index = sexs.index(chro)
    ploidy = [None] * nsample
    for i, g in enumerate(sexlist):
      if g == "F":
        ploidy[i] = str(flist[index])
      else:
        ploidy[i] = str(mlist[index])

    print("--ploidystr=%s" % ''.join(ploidy))
  except ValueError:
    print("--ploidystr=%s" % ('2' * nsample))
