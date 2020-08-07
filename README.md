# Droplet_FRAP
Droplet FRAP macros

David Brown 
2020-August-07

These [ImageJ macros](https://github.com/davidbrown2324/Droplet_FRAP/tree/master/ImageJ) and [R scripts](https://github.com/davidbrown2324/Droplet_FRAP/tree/master/R) were written to analyse Fluorescence Recovery After Photobleaching (FRAP) experiments.

There is a short [FRAPline_macro.ijm](https://github.com/davidbrown2324/Droplet_FRAP/blob/master/ImageJ/FRAPline_macro.ijm) to generate the original bleach ROI, and an alternate FRAPline_3pix_off_center.roi which was used when the objective was lightly misaligned.

The simplest script to understand the FRAP normalisation process is [FRAP_Normalisation1.ijm](https://github.com/davidbrown2324/Droplet_FRAP/blob/master/ImageJ/FRAP_Normalisation1.ijm), which performs the same calculation used in my [CFP1 paper](https://www.cell.com/cell-reports/pdf/S2211-1247(17)31137-3.pdf), and reported in [Methods in Molecular Biology](https://link.springer.com/protocol/10.1007%2F978-1-61779-477-3_11).

The file imports and file paths are ugly, and need tidying.

I hope this is useful to someone.


