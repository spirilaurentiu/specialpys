#!/bin/bash
python ../lin_bead.py --dir ./ --inFNRoots  lin_bead.in --skiprows 300 --nbins 50 --xmin -180 --xmax 180 --makehist --makeplot
if [ "$?" -eq 0 ] 
then
  echo "Test 0 for lin_bead.py passed."
fi
rm lin_bead.pdf
