q2 Python package
=================

q2 allows you to use MOOG (in its SILENT mode) to calculate elemental abundances of stars and/or determine their atmospheric parameters using the standard techniques of iron line excitation/ionization equlibrium. It also allows you to calculate other fundamental stellar parameters such as mass and age using isochrones. A tutorial is available <a href="https://github.com/astroChasqui/q2_tutorial">here</a>.

This version is mostly identical to <a href="https://github.com/astroChasqui/q2">the main repo</a> written by Ivan Ramirez, with a few under-the-hood tweaks. See that repo for more information about installing and running q2. <b>Multiprocessing must be installed for this fork to work.</b>

Functionality of this fork
------------

Here are the things I (M.B.) have changed:
- naming of the temporary files that MOOG runs with has been changed, so it is now possible to run multiple sessions of q2 on parallel processors within the same directory.
- looking at abundance trends with condensation temperature is made easier.
- galactic chemical evolution corrections can be made to the abundances once the stellar age is determined.


Contact info
------------

Megan Bedell (mbedell@oddjob.uchicago.edu)
