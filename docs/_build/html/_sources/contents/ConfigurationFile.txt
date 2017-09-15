Configuration
=============

Three configuration files
-------------------------

Three configuration files control the code.

The file ``setup.ini``
----------------------

This file includes several sections that are read by the code.

The ``Modelling`` section
^^^^^^^^^^^^^^^^^^^^^^^^^

As many line ...Models and ...DateRages as locations observatory (Earth, Spiter, K2...). Two line mandatory. If do not know, psbl_par. tmin-tmax with tmin<tmax.

The ``Archive`` section
^^^^^^^^^^^^^^^^^^^^^^^

The following options are currently supported. ::

   [Archive]
   
   Name : trial4
   Path : Archive/

Every time the code is run, a ``zip`` archive is created in the working directory. The field ``Name`` corresponds to the name given to the archive. It is allowed to choose a personal path by filling the field ``path`` [#f1]_.

In the example above, a file ``trial4.zip`` is created in the directory ``path-to-working-directory/Archive/``.


Summary
^^^^^^^

Here is an example of an actual configuration file ``setup.ini``.

::

   [EventDescription]

   Name : OB150966
   RA : 17h55m01.02s
   DEC : -29d02m49.6s


   [Observatories]

   Reference : OGLE-I

   OGLE-I : 7100-7157, 7159-7230
   Spitzer-I : 7100-7230


   [Plotting]

   Data : True
   Type : plottype1
   TimeRange : 7100-7230, 5000


   [Modelling]

   Fit : False
   Verbose : 2
   Method : grid_dmcmc

   EarthModels : psbl_par, esblparall
   EarthDateRanges : 7100-7150, 7200-7300 

   SpitzerModels : psbl_par
   SpitzerDateRanges : 7100-7150

   t0    : fit, 7205.1937
   u0    : fix, -0.01152
   tE    : fit, 57.7
   rho   : grid, 0.01, 0.1, 20
   gamma : fix, 0.0
   piEE  : fix, -0.237
   piEN  : fit, -0.0412
   s     : fix, 1.1147
   q     : fix, 0.000168
   alpha : fit, -2.2619426535897933
   dalpha : fix, 0.0
   ds : fix, 0.0

   tp : 7208.88

   thetaEN : 0.5
   thetaEE : 0.2


   [FitSetupDMCMC]

   Threads : 7
   Chains : 16
   ChainLength : 100000
   Resume : False
   Path : chains/


   [Archive]

   Name : AZERTY
   Path : Archives/





The file ``observatories.ini``
------------------------------

The ``PlotExcludedData`` section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The only option in this section is ::

   PlotExcludedData : True

that makes the excluded data semi-transparents (``True``) or remove them from the plots (``False``).

The ``ObservatoriesDetails`` section
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section includes all the observational information about the data files. As many options as data files are defined. The values of each option include the name of the observatory, the color associated to it, its location, and an additional list of data points that have to be removed from the data file.

Be carefull, the data file should be dates, magnitude, e_magn, seeing, bkground.

**Example** ::

   OGLE-I : OGLE, 000000, I, Earth, 341-372, 233
   Spitzer-I : Spitzer 0.85m, FF0000, I, Spitzer

In this example, the code will look for the files ``OGLE-I.dat`` and ``Spitzer-I.dat`` [#f2]_. It is highly recommended that the name of the file includes the filter. The corresponding observatories are ``OGLE`` and ``Spitzer 0.85m``. The colours used will be ``#000000`` and ``#FF0000``. The filters are both ``I``. The first observatory is on Earth whereas the second one is space-based. It is mandatory that the two files named ``Hori-Earth_whatever-you-want.dat`` and ``Hori-Spitzer_whatever-you-want.dat`` including ephemerids exist. Finally, in the case of OGLE data, the observations lines 341 to 372 and 233 will be removed. You can include as many data to remove as it is necessary, or nothing.


Summary
^^^^^^^

Here is an example of an actual configuration file ``setup.ini``.

::

   [PlotOptions]

   PlotExcludedData : False


   [ObservatoriesDetails]

   OGLE-I : OGLE, 000000, I, Earth, 341-372, 233
   MOA-I : MOA, 7F0000, I, Earth
   DanishLuckyCam-I : Danish 1.54m LuckyCam, 0000FF, I, Earth
   DanishDFOSC-I : Danish 1.54m DFOSC, FF7F00, I, Earth
   FaulkesNorth-I : Faulkes North 2.0m, 00FFFF, I, Earth
   FaulkesSouth-I : Faulkes South 2.0m, 007F00, I, Earth
   Liverpool-I : Liverpool 2.0m, 00A0A0, I, Earth
   MonetNorth-I : MONET North 1.2m, C0C0C0, I, Earth
   MonetSouth-I : MONET South 1.2m, BF0F00, I, Earth
   LcogtCTIOa-I : LCOGT CTIO 1m A, FF00FF, I, Earth
   LcogtCTIOb-I : LCOGT CTIO 1m B, FF00FF, I, Earth
   LcogtCTIOc-I : LCOGT CTIO 1m C, FF00FF, I, Earth
   LcogtSAAOa-I : LCOGT SAAO 1m A, FFAF00, I, Earth
   LcogtSAAOb-I : LCOGT SAAO 1m B, FFAF00, I, Earth
   LcogtSAAOc-I : LCOGT SAAO 1m C, FFAF00, I, Earth
   LcogtSSOa-I : LCOGT SSO 1m A, 7F007F, I, Earth
   LcogtSSOb-I : LCOGT SSO 1m B, 7F007F, I, Earth
   Utas-I : UTas 1.0m, C07F7F, I, Earth
   Perth-I : Perth 0.6m, 00007F, I, Earth
   SAAO-I : SAAO 1.0m, 00FF00, I, Earth
   CTIO13-I : CTIO 1.3m, 7F7F00, I, Earth
   CTIO10-I : CTIO 1.0m, 7F7F00, I, Earth
   Hereford-I : Hereford Arizona 0.35m, 007070, I, Earth
   Lemmon-I : Mt Lemmon 1.0m, B0FFB0, I, Earth
   MDM-I : MDM 2.4m, 00FFFF, I, Earth
   Palomar-I : Palomar 60'', FF7F00, I, Earth
   Regent-I : Regent Lane, 7F7FC0, I, Earth
   Possum-I : Possum 11'', FFAF00, I, Earth
   Auckland-I : Auckland 0.4m, 007F00, I, Earth
   Hunters-I : Hunters Hill 0.35m, C07F7F, I, Earth
   SouthernStars-I : Southern Stars 11'', C0C0C0, I, Earth
   FarmCove-I : Farm Cove 0.35m, 0000FF, I, Earth
   Kumeu-I : Kumeu Obs 0.35m, 00007F, I, Earth
   VintageLane-I : Vintage Lane 0.4m, 7F007F, I, Earth
   CBAPerth-I : CBA Perth 0.25m, FF00FF, I, Earth
   WiseE2V-I : Wise 1.0m E2V, BF0F00, I, Earth
   WiseSITe-I : Wise 1.0m SITe, BF0F00, I, Earth
   Bronberg-I : Bronberg 0.35m, 00A0A0, I, Earth
   Salerno-I : Salerno 0.35m, 00FF00, I, Earth
   Spitzer-I : Spitzer 0.85m, FF0000, I, Spitzer,





.. rubric:: Footnotes

.. [#f1] Do not forget the character ``/`` at the end of any path.
.. [#f2] The extention can be what ever you like.