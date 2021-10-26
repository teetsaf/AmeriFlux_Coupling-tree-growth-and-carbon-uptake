## AmeriFlux_Coupling-tree-growth-and-carbon-uptake

### Description of code

##### The code ("AmeriFlux_woodyNPP_corMat.R") will import and plot data used for the 2021 submission to JGR Biogeosciences titled "Coupling of tree growth and carbon uptake in six North American forests". The majority of the code is a large for loop that generates correlation matrices. Correlation matrices (as presented in Fig. 1) can be used to identify the temporal periods of carbon uptake with the strongest correlations to woody NPP. This has been done previously in studies linking carbon uptake with tree woody NPP to find sub-annual integrals of carbon uptake most strongly correlated to woody NPP (Babst et al., 2014b; Lagergren et al., 2019). Moving windows of carbon uptake integrations can be used to identify relationships that might not be apparent based on calendar-year integrals. Using this approach allows us to connect our results to our competing hypotheses (see Fig. 1) and visualize the relationship between numerous periods of carbon uptake woody NPP simultaneously.



### Description of datasets

### "woodyNPP.csv" is the annual woody NPP estimated from diameter reconstructions of trees in randomly allocated forest plots 
######       "site" : US.Bar = Bartlett Experimental Forest, NH; US.Ho1 = Howland Forest, ME; US.Ha1 = Harvard Forest; US.MMS = Morgan Monroe State Forest, IN; US.NR1 = Niwot Ridge, CO; and US.UMB = University of Michigan Biological Station, MI
######       "year" : year of biomass growth
######       "inc" : woody biomass increment (reported in grams of carbon m-2)
######       "sd" : standard deviation of 1000 Monte Carlo simulations of tree plots

### "daily_carbon_uptake.csv" is the gap-filled and partioned eddy covariance data daily 
######       "site" : US.Bar = Bartlett Experimental Forest, NH; US.Ho1 = Howland Forest, ME; US.Ha1 = Harvard Forest; US.MMS = Morgan Monroe State Forest, IN; US.NR1 = Niwot Ridge, CO; and US.UMB = University of Michigan Biological Station, MI
######       "year" : year of measurement
######       "doy" : day of measurement (reported in day of year, 1-365) 
######       "NEE_U05" : daily net ecosystem exchange, raw data filtered based on the 5% quantile of the estimated u* distributions generated in REddyProc package (Wutzler et al. 2016) (reported in grams of carbon m-2)
######       "NEE_U50" : daily net ecosystem exchange, raw data filtered based on the 50% quantile of the estimated u* distributions generated in REddyProc package (Wutzler et al. 2016) (reported in grams of carbon m-2)
######       "NEE_U95" : daily net ecosystem exchange, raw data filtered based on the 95% quantile of the estimated u* distributions generated in REddyProc package (Wutzler et al. 2016) (reported in grams of carbon m-2)     
######       "GPP_U05" : daily gross primary productiviy, raw data filtered based on the 5% quantile of the estimated u* distributions generated in REddyProc package (Wutzler et al. 2016) (reported in grams of carbon m-2)
######       "GPP_U50" : daily gross primary productiviy, raw data filtered based on the 50% quantile of the estimated u* distributions generated in REddyProc package (Wutzler et al. 2016) (reported in grams of carbon m-2)
######       "GPP_U95" : daily gross primary productiviy, raw data filtered based on the 95% quantile of the estimated u* distributions generated in REddyProc package (Wutzler et al. 2016) (reported in grams of carbon m-2)     
######       "Reco_U05" : daily ecosystem respiration, raw data filtered based on the 5% quantile of the estimated u* distributions generated in REddyProc package (Wutzler et al. 2016) (reported in grams of carbon m-2)
######       "Reco_U50" : daily ecosystem respiration, raw data filtered based on the 50% quantile of the estimated u* distributions generated in REddyProc package (Wutzler et al. 2016) (reported in grams of carbon m-2)
######       "Reco_U95" : daily ecosystem respiration, raw data filtered based on the 95% quantile of the estimated u* distributions generated in REddyProc package (Wutzler et al. 2016) (reported in grams of carbon m-2)     



### "monthly_climate.csv" is the monthly climate data for each site 
######       "site" : US.Bar = Bartlett Experimental Forest, NH; US.Ho1 = Howland Forest, ME; US.Ha1 = Harvard Forest; US.MMS = Morgan Monroe State Forest, IN; US.NR1 = Niwot Ridge, CO; and US.UMB = University of Michigan Biological Station, MI
######       "year" : year of measurement
######       "month" : month of measurement (1-12 corresponding to Jan-Dec)
######       "date" : YYYY-MM
######       "temp" : gap filled mean monthly temperature (degrees C) measured from the eddy covariance tower
######       "Rg" : gapf-filled mean monthly incoming short-wave radiation (W m-2) measured from the eddy covariance tower
######       "precip" : monthly precipitation totals from the PRISM interpolated dataset (Daly et al., 2008)



### Citations

##### Babst F, Bouriaud O, Papale D, Gielen B, Janssens I a., Nikinmaa E, Ibrom A, Wu J, Bernhofer C, Köstner B, et al. 2014b. Above-ground woody carbon sequestration measured from tree rings is coherent with net ecosystem productivity at five eddy-covariance sites. New Phytologist 201: 1289–1303.

##### Daly C, Halbleib M, Smith JI, Gibson WP, Doggett MK, Taylor GH, Pasteris PP. 2008. Physiographically sensitive mapping of climatological temperature and precipitation across the conterminous United States. International Journal of Climatology.

##### Lagergren F, Jönsson AM, Linderson H, Lindroth A. 2019. Time shift between net and gross CO2 uptake and growth derived from tree rings in pine and spruce. Trees - Structure and Function 33: 765–776.

##### Wutzler T, Lucas-Moffat A, Migliavacca M, Knauer J, Sickel K, Šigut L, Menzer O, Reichstein M. 2018. Basic and extensible post-processing of eddy covariance flux data with REddyProc. Biogeosciences 15: 5015–5030.

