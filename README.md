# soci_gnc_groundpass

Example output for an ISS orbit:
![Example output for an ISS orbit](figs/sample_groundtrack.png)

## Known issues
- the ground station circle that is plotted is an approximation computed using the s/c orbit's semimajor axis. The actual altitude of the s/c when above the ground station will affect this curve. Hence a ground pass will not start or end at the exact instant the s/c crosses the cyan line, but it will be close.
