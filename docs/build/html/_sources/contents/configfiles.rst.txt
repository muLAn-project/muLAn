Configuration files
===================

**muLAn** is designed to use two parameters files: one file
``setup.ini`` that includes all the basic parameters a modeler wishes
to change often, and one file ``advanced_setup.ini`` that includes all
the information that should be chosen only one, before starting to
work on an event.

Basic parameters 
----------------

List of the sections
^^^^^^^^^^^^^^^^^^^^

In the file ``setup.ini``, the parameters are divided in 
the following sections:

- ``EventDescription``: information about the event;
- ``Controls``: choice of the general mode of **muLAn**;
- ``Observatories``: data to include in the modeling, and rescaling
  factors;
- ``Plotting``: plotting mode and corresponding options; 
- ``Modelling``: physical parameters related to the fit;
- ``FitSetupDMCMC``: computational parameters related to the fit. 

Full Explained Example
^^^^^^^^^^^^^^^^^^^^^^

In the file, a section starts with its name following the syntax:
``[NameOfSection]``. Here below is an example of a full ``setup.ini``
file. ::

   [EventDescription]
   
   Name : OB151737      # name of the event
   RA   : 18h08m05.83s  # right ascention
   DEC  : -25d36m41.0s  # declination
    
   [Controls]
   
   Modes   : Fit, Plot  # fit and plot the results 
   Archive : run1       # name of the fit
   
   [Observatories]
   # All dataset to be used must be included in this section.
   
   # The line below includes all data in the time interval [6900;7640] of the file OGLE-I.dat
   OGLE-I : (1.0, 0.0), 6900-7640
   
   [Plotting]
   
   Data   : True  # True | False, wether or not plot data
   Models :       # list of models to plot, e.g. "1, 5, 7" to plot 
                  # models 1, 5 and 7, or "1-5" to plot all the models 
                  # from 1 to 5. 
   
   Type    :          fitmag,         fitflux  # name of the plot routines
   Options : 6900-7640/5000/, 6900-7640/5000/  # options for the plot routines
                                               # one or several plot routines can be listed
   
   [Modelling]
   
   Verbose : 5           # from 0 (print nothing) to 5 (print everything)
   Method  : grid_dmcmc  # exploration method, only one option for now
   
   # Line below shows how several algorithms to compute the magnification can be used:
   # * Models_Earth means a telescope in the center of the Earth 
   # * 0-10000/PSPL means that the method PSPL is used in [0;10000] time interval
   # Important: for overlaping intervals, priority is given to the right-most algorithm.
   Models_Earth : 0-10000/PSPL, 7000-7400/fhexaBL, 7220-7300/BLcontU
   
   # Lines below shows how the initial conditions are chosen:
   #    param : fit, a, b, m
   # means using a differential MCMC with a random initial condition
   # in the interval [m-a;m+a]
   t0     : fit,  1e-3, 1e-3, 7.2657375353e+03
   u0     : fit,  1e-4, 1e-4, 4.0952516233e-02
   tE     : fit,  3e-2, 3e-2, 5.4396069563e+01
   rho    : fit,  1e-5, 1e-5, 7.3298960206e-04
   gamma  : fix,      ,     , 0.0
   piEN   : fit,  1e-3, 1e-3, 0.0
   piEE   : fit,  1e-3, 1e-3, 0.0
   s      : fit,  3e-2, 3e-2, 2.2277240797e+00
   q      : fit,  1e-3, 1e-3, 7.7902505266e-02
   alpha  : fit,  1e-5, 1e-5, 1.3903531793e+00
   dalpha : fit,  1e-3, 1e-3, 0.0
   ds     : fit,  1e-3, 1e-3, 0.0
   
   tp   : 7265.0
   tb : 7266.0
   
   thetaEN : fix, -2, 2, 0.5  # currently not used
   thetaEE : fix, -2, 2, 0.2  # currently not used
   
   [FitSetupDMCMC]
   
   Threads : 24         # any number of threads
   Chains : 180         # number of parallel chains in MCMC
   ChainLength : 10000  # maximum length of each chain
   Resume : False       # wether or not resume a run

Advanced parameters
-------------------

List of the sections
^^^^^^^^^^^^^^^^^^^^

In the file ``advancedsetup.ini``, the parameters are divided in 
the following sections:

- ``FullPaths``: information about the paths;
- ``RelativePaths``: paths of the inputs/outputs.

Full Explained Example
^^^^^^^^^^^^^^^^^^^^^^

In the file, a section starts with its name following the syntax:
``[NameOfSection]``. Here below is an example of a full ``advancedsetup.ini``
file. ::

   [FullPaths]
   
   Code : path/to/muLAn/ 
   PathLocalPythonPackages :  # not used any longer
   
   [RelativePaths]
   
   Data : Data/
   Plots : Plots/
   Chains : Chains/
   Outputs : Outputs/
   Archives : Archives/
   ModelsHistory : ./
   
   # The next sections correspond to the advanced parameters used by
   # the modules computing the magnification.
   
   [ESBLauto]
   
   TriggerNextMethod : 1e-3
   PrecisionGoalRayshooting : 1e-3
   LocalMarginRayshooting : 3
   
   [IRS]
   
   PrecisionGoalRayshooting : 1e-3
   LocalMarginRayshooting : 2

Observatory properties and data definition
------------------------------------------

List of the sections
^^^^^^^^^^^^^^^^^^^^

In the file ``observatories.ini``, the parameters are divided in 
the following sections:

- ``PlotOptions``: options related to the plots;
- ``ObservatoriesDetails``: information about each *dataset*.

Full Explained Example
^^^^^^^^^^^^^^^^^^^^^^

In the file, a section starts with its name following the syntax:
``[NameOfSection]``. Here below is an example of a full
``observatories.ini`` file. ::

   [PlotOptions]
   
   PlotExcludedData : False  # wether or not plot the excuded data
   
   [ObservatoriesDetails]
   
   # The format for each dataset is:
   #    Name: displayed name, HTML color, filter, "Magnitude" | "Flux", reference frame, list
   # The list is a sequence of rejected data (the number refers to the data ID in the data file.
   OGLE-I : OGLE, 000000, I, Magnitude, Earth
   OGLE-V : OGLE, 000000, I, Magnitude, Earth
   MOA-I : MOA, 7F0000, I, Flux, Earth, 1572, 1575, 1581, 1584
   MOA-V : MOA, 7F0000, V, Magnitude, Earth
   KMTCTIO-I : KMT CTIO 1.6m, 09364E, I, Magnitude, Earth
   KMTCTIO-V : KMT CTIO 1.6m, 3DAE11, V, Magnitude, Earth
   KMTSAAO-I : KMT SAAO 1.6m, 8A084B, I, Magnitude, Earth
   KMTSSO-I : KMT SSO, 174405, I, Magnitude, Earth
   DanishLuckyCam-I : Danish 1.54m LuckyCam A, 0000FF, I, Magnitude, Earth
   DanishLuckyCam-Z : Danish 1.54m LuckyCam B, 0000FF, Z, Magnitude, Earth
   DanishDFOSC-I : Danish 1.54m DFOSC, FF7F00, I, Magnitude, Earth
   FaulkesNorth-I : Faulkes North 2.0m, 00FFFF, I, Magnitude, Earth
   FaulkesSouth-I : Faulkes South 2.0m, 007F00, I, Magnitude, Earth
   Liverpool-I : Liverpool 2.0m, 00A0A0, I, Magnitude, Earth
   MonetNorth-I : MONET North 1.2m, C0C0C0, I, Magnitude, Earth
   MonetSouth-I : MONET South 1.2m, BF0F00, I, Magnitude, Earth
   LcogtCTIOa-I : LCOGT CTIO 1m A, FF00FF, I, Magnitude, Earth
   LcogtCTIOb-I : LCOGT CTIO 1m B, FF00FF, I, Magnitude, Earth
   LcogtCTIOc-I : LCOGT CTIO 1m C, FF00FF, I, Magnitude, Earth
   LcogtSAAOa-I : LCOGT SAAO 1m A, FFAF00, I, Magnitude, Earth
   LcogtSAAOb-I : LCOGT SAAO 1m B, FFAF00, I, Magnitude, Earth
   LcogtSAAOc-I : LCOGT SAAO 1m C, FFAF00, I, Magnitude, Earth
   LcogtSSOa-I : LCOGT SSO 1m A, 7F007F, I, Magnitude, Earth
   LcogtSSOb-I : LCOGT SSO 1m B, 7F007F, I, Magnitude, Earth
   Utas-I : UTas 1.0m, C07F7F, I, Magnitude, Earth
   Perth-I : Perth 0.6m, 00007F, I, Magnitude, Earth
   SAAO-I : SAAO 1.0m, 00FF00, I, Magnitude, Earth
   CTIO13-I : CTIO 1.3m, 7F7F00, I, Magnitude, Earth
   CTIO10-I : CTIO 1.0m, 7F7F00, I, Magnitude, Earth
   Hereford-I : Hereford Arizona 0.35m, 007070, I, Magnitude, Earth
   Lemmon-I : Mt Lemmon 1.0m, B0FFB0, I, Magnitude, Earth
   MDM-I : MDM 2.4m, 00FFFF, I, Magnitude, Earth
   Palomar-I : Palomar 60'', FF7F00, I, Magnitude, Earth
   Regent-I : Regent Lane, 7F7FC0, I, Magnitude, Earth
   Possum-I : Possum 11'', FFAF00, I, Magnitude, Earth
   Auckland-I : Auckland 0.4m, 007F00, I, Magnitude, Earth
   Hunters-I : Hunters Hill 0.35m, C07F7F, I, Magnitude, Earth
   SouthernStars-I : Southern Stars 11'', C0C0C0, I, Magnitude, Earth
   FarmCove-I : Farm Cove 0.35m, 0000FF, I, Magnitude, Earth
   Kumeu-I : Kumeu Obs 0.35m, 00007F, I, Magnitude, Earth
   VintageLane-I : Vintage Lane 0.4m, 7F007F, I, Magnitude, Earth
   CBAPerth-I : CBA Perth 0.25m, FF00FF, I, Magnitude, Earth
   WiseE2V-I : Wise 1.0m E2V, BF0F00, I, Magnitude, Earth
   WiseSITe-I : Wise 1.0m SITe, BF0F00, I, Magnitude, Earth
   Bronberg-I : Bronberg 0.35m, 00A0A0, I, Magnitude, Earth
   Salerno-I : Salerno 0.35m, 00FF00, I, Magnitude, Earth
   Spitzer-I : Spitzer 0.85m, FF0000, I, Magnitude, Spitzer
   LcogtSBIG1m-I : LCOGT SBIG cameras 1m, FFAF00, I, Magnitude, Earth
   LcogtSinistro1m-I : LCOGT Sinistro cameras 1m, FF00FF, I, Magnitude, Earth
   Lcogt2m-I : LCOGT 2m, 00FFFF, I, Magnitude, Earth
   


