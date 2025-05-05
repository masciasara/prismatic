#################################################
# README file for CAPERS master yield files v2. #
#################################################

Complete yield files compiling all the targets included in the MSA in any of 
the observed pointings, separated by survey field (UDS, EGS, COSMOS).

When MSA overlap between pointings is present, the weight of any previously 
selected source in the overlapping region was downgraded to avoid duplication,
except for top-priority targets.

Note that the format of the yield file is simpler for the UDS observations, 
where there is no MSA overlap at all. For the EGS and COSMOS observations, the
yield file is a semicolon-separated table in an attempt to minimize the number
of required columns to contain all the possible pointing and/or configuration
combinations.

The UDS yield file v2 includes only the sources in the observations that were
successfully executed in December 2024 – January 2025. About 60% of the 
scheduled UDS observations were missed for several reasons, including a JWST
star tracker failure, MSA target acquisition failure, and limited JWST 
communications due to the evaculation of JPL during the Januray 2025 California 
wildfires. The skipped observations will be rescheduled in the future.

Observation no. 60 for EGS pointing 6 also failed due to a problem with one of
JWST's fine guidance sensors. The EGS yield file here excludes pointing 6, 
which will also be replanned in the future.

Unexpected stuck closed shutters occurring during the UDS observations have
been identified and properly accounted for in the relevant columns 
(n_nods_visX, eff_exp_time).
The yield files for EGS and COSMOS will be updated after those observations are 
completed and these stuck closed shutters are identified.


- Column descriptions:

                 ID: MSA ID in the observations.
                 ra: Right Ascension.
                dec: Declination.
              z_med: Only for very high redshift candidates. Median photometric
                     redshift among the contributed high-z samples including 
                     the source.
          z_UNICORN: Photometric redshift in the UNICORN catalogs.
           flux_XXX: Flux in AB mag.
         MSA_weight: Assigned linear weight during the MSA planning.
           pointing: Integer or array indicating the pointing/s in which the 
                     source is observed.
        n_nods_visX: Integer or array indicating the number of nods (or just 
                     exposure units) in which the source is observed. If an 
                     array is given, the values correspond to the different 
                     pointings in the pointing array. The maximum number of 
                     exposure units a source can get is 6, accounting for 
                     three nods in two different (dithered or not) exposures.
                     Note that dithering was employed in some of the pointings
                     but not all of them. Moreover, there are rare cases where
                     a source gets shifted out of its slitlet when dithering,
                     which can result in odd n_nods_visX values. If the 
                     unexpected stuck closed shutters during the observations
                     have been identified (see above), the effect of those 
                     are also taken into account for the n_nods_visX values 
                     here reported.
       eff_exp_time: Total effective exposure time (s) the source receives 
                     combining all the pointings and configurations in which
                     it is observed.
  shutter_centering: String or array of strings indicating the centering 
                     constraint of the source within its shutter, following the
                     MPT nomenclature. If an array is given, the values 
                     correspond to the different pointings in the pointing 
                     array. The MSA configurations were first optimized 
                     using a 'Midpoint' constraint and some extra sources were
                     afterwards included in "Entire Open Shutter Area" ('EOSA')
                     constraint (less centered within their shutters).