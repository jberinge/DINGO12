version_name = "OzFluxQC"
version_number = "V2.7.0"
# V2.7.0   - major bug fixes as for V2.6.3 above
#          - minor bug fixes to clean up use of switches in ['Options'] section
#          - minor fixes to check that requested files exist
#          - implemented compare_ep.py to automate comparison of OzFluxQC and
#            EddyPro results
# V2.6.3 - clean up of code after comparing EddyPro and OzFluxQC
#          - fix bugs in mf.vapourpressure and mf.RHfromabsolutehumidity
#          - deprecated WPLcov after finding an error in Fc_WPLcov (units of wA term)
#          - rationalised calculation of dry and moist air densities and parial
#            density of water vapour
#          - implemented EddyPro method of calculating rho*Cp
#          - tidyied up use of densities in WPL correction and used flux
#            form (WPL80 Eqn 42a and 44) for Fe and Fc
#          - implemented EddyPro method of calculating Fh from Fhv
#          - rationalised use of densities and rho*Cp when calculating fluxes
# V2.6.2 - implement Lloyd-Taylor Reco
# V2.6.1 - fix 2D coordinate rotation for momentum covariances
# V2.6.0 - fix ConvertCO2 units bug
# V2.5.2 - implement sofm/solo/seqsolo
# V2.5.1 - post-Cairns 2013 version
