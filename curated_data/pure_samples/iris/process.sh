grep 'U133B' GSE22886_descriptions.txt  | awk '{print "rm " $1 ".CEL.gz"}' | sh