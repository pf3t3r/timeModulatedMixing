# timeModulatedMixing
Time modulation of mixing at moorings IC3 and M1 over Reykjanes Ridge.

Running Instructions

1. Ensure that the following data files are placed in the Data folder:
merged_hourly_unfiltered_data_20142020.mat
Mixing_parameterization_fields.mat
WOCE_climatology.mat

2. Find the tide
extractVelocity -> extractTide

3. Find the buoyancy frequency
extractSalinity + extractTemperature -> extractBuoyancy

4. Find the near-inertial wave energy:
extractNiw
(extractTide should be run before this)

5. Apply the time-modulation to the parametrisation of HIL and NIW
paramEhil -> param_Fzt -> param_EpsHil -> param_Diffusivity
param_Ewwi -> param_Fzt_wwi -> param_EpsWwi -> param_Diffusivity_wwi
