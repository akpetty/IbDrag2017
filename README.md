## IceBridge Arctic sea ice drag scripts

Data processing and plotting used to estimate the neutral drag coefficient over Arctic sea ice using high-resolution IceBridge laser (ATM) data, from the following publication (in review currently):

Petty, A. A., M. C. Tsamados, N. T. Kurtz, Atmospheric form drag over Arctic sea ice using remotely sensed ice topography observations, JGR Earth's Surface (in review)

Note that individual descriptions should be included at the top of each script.

The processed data can be downloaded from X.
Place this in the DataOutput folder. Note that this data is in a binary Python format so data readers (included in IB_functions.py) are needed.

The 2D ice topography dataset is from an earlier publication: Petty, A. A., M. C. Tsamados, N. T. Kurtz, S. L. Farrell, T. Newman, J. Harbeck, D. L. Feltham, J. A. Richter-Menge (2016), Characterizing Arctic sea ice topography using high-resolution IceBridge data, The Cryosphere, doi:10.5194/tcd-9-6495-2015.

The scripts used in that publication are hosted on another GitHub repo: https://github.com/akpetty/ibtopo2016

The following raw IceBridge datasets are used to create the 1D and 2D topography data:
The IceBridge ATM data: https://nsidc.org/data/docs/daac/icebridge/ilatm1b/. 
The IceBridge DMS imagery: http://nsidc.org/data/iodms1b}. 
The IceBridge IDCSI4 and quick-look sea ice data: http://nsidcorg/data/docs/daac/icebridge/evaluation_products/sea ice-freeboard-snowdepth-thickness-quicklook-index.html and http://nsidc.org/data/idcsi4.html. 

For some of the extra processing/plotting, the following datasets are also required:
The daily OSI-SAF ice type data: http://saf.met.no/p/ice/
The nearest coastline proximity data: http://oceancolor.gsfc.nasa.gov/DOCS/DistFromCoast/.

which have been included in the data link given above, along with the DMS image used in figure 3.


Note also that Python 2.7 was used for all processing. I have not tested these scripts in Python 3.



