# Landslide_Tsunami_Benchmark

The objective of this repo is to provide programs and experimental data for the testing
and verification of dispersive wave models for submarine landslide generated tsunami.
The experiment is based on the research and thesis of Langford Sue at the University of
Canterbury, Christchurch, New Zealand. The basic experiment is an elliptical slider
that travels down a sloping surface in a water-filled flume. Each experiment is a
different combination of slider weight and initial location. 

While there are many experiments that use sliders, there are two important attributes 
that make this experiment unique. First, the slope is continuously variable allowing
the slider to runout without an abrupt termination at the base of the slope. Hence, the 
dispersive wavefield can be followed along the entire flume and resolve the propagating
wave packet. Second, a continuous surface elevation is measured over time rather than
discrete points at water level gauges. This feature allows additional analysis of phase
and group velocity, as well as evolution of the wave packet.

In summary, the experiments provide a high resolution framework in which to test
numerical models.

This repo is organised into the following folders:

code: Contains programs to reproduce the slider motion.

data: Contains the experimental data for comparison with model results.

docs: Contains Langford Sue's thesis with complete details of the experments and a 
   shorter summary paper.

grid: Contains a spreadsheet with the equations and output for calculationg the flume 
   bottom elevation and a sample grid.
