The fitler sets

Decam
=====

Decam filter transmission tables::

      DECam_filters.xlsx: DECam_filters.xlsxTable with total
      throughput for u,g,r,i,z,Y filters (thanks to William Wester
      (DES Collaboration)). If you want the system throughput without
      the atmosphere, divide by the "atm" column

      DECam_filters_transmission.txt: Filter transmission (area
      weighted response) from Asahi for u,g,r,i,z,Y and VR (plot)

These data come from the following web page::

      http://www.ctio.noao.edu/noao/content/Dark-Energy-Camera-DECam

HSC
===

Filter transmission come from the following web page::

       [https://www.naoj.org/Projects/HSC/forobservers.html](https://www.naoj.org/Projects/HSC/forobservers.html)

The following plot comes from the same webpage

[](https://www.naoj.org/Projects/HSC/filterData/fig.png)

Megacam
=======

Megacam 1
---------

Original filter set for MegaCam. First column of the files used:

> through airmass 1.3, point course

Megacam 2
---------

New `i` filter (`i2`). See here for details::

    [http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/ifilt.html](http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/ifilt.html)
    
Megacam 3
---------

From http://www.cfht.hawaii.edu/Instruments/Filters/megaprime.html

All the (ugriz)'

u' 9302 8/18/2014
g' 9402 10/21/2014
r' 9602 9/8/2014
i' 9703 9/8/2014
z' 9901 10/23/2014

http://www.cfht.hawaii.edu/Instruments/Imaging/Megacam/specsinformation.html#P2

> Each filter was scanned 17 times by the manufacturer at different places

These 17 measruements correspond to the 17 columns of the filter data
files. The average transmission will thus be used. The `(ugriz).dat`
files correspond to the average over the different measurements for
each filter..

Pan-Starrs
==========

From here::

     [https://confluence.stsci.edu/display/PANSTARRS/PS1+Filter+properties#PS1Filterproperties-Filterdescriptions](https://confluence.stsci.edu/display/PANSTARRS/PS1+Filter+properties#PS1Filterproperties-Filterdescriptions)
