#!/bin/sh -e

if [ "$FASTERCAP" != "" ]; then
  fastercap=$(realpath $FASTERCAP)
else
  fastercap=FasterCap
  if ! which $fastercap >/dev/null 2>&1; then
    echo "*** ERROR: could not find $fastercap in PATH"
    exit 1
  else
    fastercap=$(which $fastercap)
  fi
fi

if [ "$KLAYOUT" != "" ]; then
  klayout=$(realpath $KLAYOUT)
else
  klayout=klayout
  if ! which $klayout >/dev/null 2>&1; then
    echo "*** ERROR: could not find $klayout in PATH"
    exit 1
  else
    klayout=$(which $klayout)
  fi
fi

cd $(dirname $(which $0))

$klayout -b \
    -r ../lvs_sky130_fastercap.lylvs \
    vppcap.gds \
    -rd fastercap="$fastercap" \
    -rd no_lvs=1 \
    -rd schematic=vppcap.cir \
    -rd no_simplify=1 

