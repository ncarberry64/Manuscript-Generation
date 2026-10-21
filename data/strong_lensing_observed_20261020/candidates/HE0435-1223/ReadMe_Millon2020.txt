J/A+A/640/A105      R-band light curves of 23 lensed QSOs        (Millon+, 2020)
================================================================================
COSMOGRAIL XIX: Time delays in 18 strongly lensed quasars from 15 years of
optical monitoring.
    Millon M., Courbin F., Bonvin V., Paic E., Meylan G., Tewes M., Sluse D.,
    Magain P., Chan J.H.H., Galan A., Joseph R., Lemon C., Tihhonova O.,
    Anderson R.I., Marmier M., Chazelas B., Lendl M., Triaud A.H.M.J.,
    Wyttenbach A.
    <Astron. Astrophys. 640, A105 (2020)>
    =2020A&A...640A.105M        (SIMBAD/NED BibCode)
================================================================================
ADC_Keywords: QSOs ; Gravitational lensing ; Photometry
Keywords: methods: data analysis - gravitational lensing: strong -
          cosmological parameters

Abstract:
    We present the results of 15 years of monitoring lensed quasars, which
    was conducted by the COSMOGRAIL programme at the Leonhard Euler 1.2m
    Swiss Telescope. The decade-long light curves of 23 lensed systems are
    presented for the first time. We complement our data set with other
    monitoring data available in the literature to measure the time delays
    in 18 systems, among which nine reach a relative precision better than
    15% for at least one time delay. To achieve this, we developed an
    automated version of the curve-shifting toolbox PyCS to ensure robust
    estimation of the time delay in the presence of microlensing, while
    accounting for the errors due to the imperfect representation of
    microlensing. We also re-analysed the previously published time delays
    of RXj1131-1231 and HE0435-1223, by adding six and two new seasons of
    monitoring, respectively, and confirming the previous time-delay
    measurements. When the time delay measurement is possible, we
    corrected the light curves of the lensed images from their time delay
    and present the difference curves to highlight the microlensing signal
    contained in the data. To date, this is the largest sample of
    decade-long lens monitoring data, which is useful to measure $H_0$ and
    the size of quasar accretion discs with microlensing as well as to
    study quasar variability.

Description:
    R-band cosmograil light curves of lensed quasars taken at the Euler
    1.2m Swiss telescope with the EulerCAM and EulerC2 instrument.
    Additional data taken at the SMARTS 1.3m telescope with the ANDICAM
    optical/ infrared camera are also included for SDSS J0924+0219
    (MacLeod et al., 2015ApJ...806..258M, Cat. J/ApJ/806/258).

    For HE 0435-1223 and RX J1131-1231, data taken at the 1.5m telescope
    at the Maidanak Observatory, the Mercator Belgian telescope and at the
    SMARTS 1.3 m telescope are included (Courbin et al.
    2011A&A...536A..53C, Cat, J/A+A/536/A53, Tewes et al.
    2013A&A...553A.120T, 2013A&A...556A..22T, Cat. J/A+A/556/A22). For
    SDSS J1515+1511, we had the data from Shalyapin & Goicoechea
    (2017ApJ...836...14S) from the Liverpool Telescope.

File Summary:
--------------------------------------------------------------------------------
 FileName      Lrecl  Records   Explanations
--------------------------------------------------------------------------------
ReadMe            80        .   This file
stars.dat        100       23   List of studied stars
lcab/*             .       17   R-band light curves with A and B images
lcabcd/*           .        4   R-band light curves with images ABCD or ABC
lcapd/*            .        1   R-band light curves of SDSS J0924+0219
lcapb/*            .        1   R-band light curves of HE 0230-2130
--------------------------------------------------------------------------------

See also:
 J/A+A/488/481    : Time delays and H0 from WFI J2033-4723 (Vuissoz+, 2008)
 J/A+A/536/A53    : R-band light curves of HE 0435-1223 (Courbin+, 2011)
 J/A+A/553/A121   : 2 lensed quasars light curves (Eulaers+, 2013)
 J/A+A/556/A22    : R-band light curves of RX J1131-1231 (Tewes+, 2013)
 J/A+A/557/A44    : Doubly lensed quasar SDSS J1001+5027 (Rathna Kumar+, 2013)
 J/MNRAS/465/4914 : R-band light curves of HE 0435-1223 (Bonvin+, 2017)
 J/A+A/609/A71    : R-band light curves of DES J0408-5359 (Courbin+ 2018)
 J/A+A/616/A183   : R-band light curves of PG 1115+080 (Bonvin+, 2018)
 J/A+A/629/A97    : R-band light curves of WFI2033-4723 (Bonvin+, 2019)

Byte-by-byte Description of file: stars.dat
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1- 18  A18   ---     Name      Star name
  20- 21  I2    h       RAh       Simbad right ascension (J2000)
  23- 24  I2    min     RAm       Simbad right ascension (J2000)
  26- 30  F5.2  s       RAs       Simbad right ascension (J2000)
      32  A1    ---     DE-       Simbad declination sign (J2000)
  33- 34  I2    deg     DEd       Simbad declination (J2000)
  36- 37  I2    arcmin  DEm       Simbad declination (J2000)
  39- 42  F4.1  arcsec  DEs       Simbad declination (J2000)
  44- 67  A24   ----    FileName  Name of the R light curve file
  69-100  A32   ---     SName     Simbad name
--------------------------------------------------------------------------------

Byte-by-byte Description of file: lcab/*
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1- 11  F11.5 d       HJD       Median Heliocentric Julian Date (HJD-2400000)
  13- 21  F9.5  mag     magA      Relative magnitude of image A
  23- 29  F7.5  mag   e_magA      Relative magnitude error of image A
  31- 39  F9.5  mag     magB      Relative magnitude of image B
  41- 47  F7.5  mag   e_magB      Relative magnitude error of image B
  49- 59  A11   ---     Tel       Telescope used for the observation
--------------------------------------------------------------------------------

Byte-by-byte Description of file: lcabcd/*
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1- 11  F11.5 d       HJD       Median Heliocentric Julian Date  (HJD-2400000)
  13- 21  F9.5  mag     magA      Relative magnitude of image A
  23- 29  F7.5  mag   e_magA      Relative magnitude error of image A
  31- 39  F9.5  mag     magB      Relative magnitude of image B
  41- 47  F7.5  mag   e_magB      Relative magnitude error of image B
  49- 57  F9.5  mag     magC      Relative magnitude of image C
  59- 65  F7.5  mag   e_magC      Relative magnitude error of image C
  67- 75  F9.5  mag     magD      ? Relative magnitude of image D (1)
  77- 83  F7.5  mag   e_magD      ? Relative magnitude error of image D  (1)
  85- 97  A13   ---     Tel       Telescope used for the observation
--------------------------------------------------------------------------------
Note (1): Not for WFI2026_Euler.dat
--------------------------------------------------------------------------------

Byte-by-byte Description of file: lcapd/*
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1- 11  F11.5 d       HJD       Median Heliocentric Julian Date  (HJD-2400000)
  13- 21  F9.5  mag     magA+D    Relative magnitude of image A+D
  23- 29  F7.5  mag   e_magA+D    Relative magnitude error of image A+D
  31- 39  F9.5  mag     magB      Relative magnitude of image B
  41- 47  F7.5  mag   e_magB      Relative magnitude error of image B
  49- 57  F9.5  mag     magC      Relative magnitude of image C
  59- 65  F7.5  mag   e_magC      Relative magnitude error of image C
  67- 79  A13   ---     Tel       Telescope used for the observation
--------------------------------------------------------------------------------

Byte-by-byte Description of file: lcapb/*
--------------------------------------------------------------------------------
   Bytes Format Units   Label     Explanations
--------------------------------------------------------------------------------
   1- 11  F11.5 d       HJD       Median Heliocentric Julian Date  (HJD-2400000)
  13- 21  F9.5  mag     magA+B    Relative magnitude of image A+B (1)
  23- 29  F7.5  mag   e_magA+B    Relative magnitude error of image A+B (1)
  31- 38  F8.5  mag     magC      Relative magnitude of image C
  40- 46  F7.5  mag   e_magC      Relative magnitude error of image C
  48- 55  F8.5  mag     magD      Relative magnitude of image D
  57- 63  F7.5  mag   e_magD      Relative magnitude error of image D
  65- 77  A13   ---     Tel       Telescope used for the observation
--------------------------------------------------------------------------------
Note (1): Flux of image A and B have been summed due to low image separation.
--------------------------------------------------------------------------------

Acknowledgements:
    Martin Millon, martin.millon(at)epfl.ch

References:
 Eigenbrod et al.,    Paper I     2005A&A...436...25E
 Eigenbrod et al.,    Paper II    2006A&A...451..747E
 Eigenbrod et al.,    Paper III   2006A&A...451..759E
 Sama et al.,         Paper IV    2006A&A...450..461S
 Vuissoz et al.,      Paper V     2007A&A...464..845V
 Eigenbrod et al.,    Paper VI    2007A&A...465...51E
 Vuissoz et al.,      Paper VII   2008A&A...488..481V, Cat. J/A+A/488/481
 Chantry et al.,      Paper VIII  2010A&A...522A..95C
 Courbin et al.,      Paper IX    2011A&A...536A..53C, Cat. J/A+A/536/A53
 Sluse et al.,        Paper X     2012A&A...538A..99S
 Tewes et al.,        Paper XI    2013A&A...553A.120T
 Eulaers et al.,      Paper XII   2013A&A...553A.121E, Cat. J/A+A/553/A121
 Tewes et al.,        Paper XIII  2013A&A...556A..22T, Cat. J/A+A/556/A22
 Rathna Kumar et al., Paper XIV   2013A&A...557A..44R, Cat. J/A+A/557/A44
 Bonvin et al.,       Paper XV    2016A&A...585A..88B
 Courbin et al.,      Paper XVI   2018A&A...609A..71C, Cat. J/A+A/609/A71
 Bonvin et al.,       Paper XVII  2018A&A...616A.183B, Cat. J/A+A/616/A183
 Bonvin et al.,       Paper XVIII 2019A&A...629A..97B, Cat. J/A+A/629/A97

================================================================================
(End)  Martin Millon [EPFL, Switzerland], Patricia Vannier [CDS]    16-Jun-2020
