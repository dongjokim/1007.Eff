GREPARG=${1:-"effTrigg-"}
FIGDIR=${2:-"figs"}
mkdir -p $FIGDIR
ls eff/* | grep -i "$GRAPDIR" | while read x; do echo root -b -q -l DrawEffTrigg.C\'\(\"$x\",1\)\';done
