

**A (huge) collection of interchangeable empirical models to compute the thermosphere density, temperature and composition**

### Description

Thermosphere models are usually employed in space activities such as mission planning and analysis, orbit determination, attitude control design, maneuvering planning, fuel estimation and orbit decay studies, among others.

The Orbital Dynamics group of INPE coded during the 80’s several thermospheric computer programs in FORTRAN. Those models were based on Jacchia works (70 and 77) and their analytical versions, since Jacchia assumed numeric integration of the diffusion equation in order to compute density and composition as function of height. Numeric integration normally demands computing time, and therefore analytical means “fast” in this sense.

Presently there is a set of 5 thermospheric models:

* Jacchia 70 (in fact it is almost the same as the 64 model)
* Jacchia-Roberts
* Jacchia 77
* LaFontaine-Hughes
* Jacchia-Lineberry (77 model)

 All these programs are interchangeable, which means that any model, can be easily changed to any other. Of course there are several models not included in this list, but, with few exceptions, they do not satisfy the requirement of "interchangeability" since they need different input data, like, for instance, the Jacchia-Bowman JB2006 and JB2008 models ([Bowman et al. 2008a](#Bowman_2006), [2008b](#Bowman_2008)), which considers the effect of solar UV radiation.

The interchangeable models present 3 entry points:

* Static profile (function of the exospheric temperature and height)
* Static density (same as static profile, but function of solar flux F10.7 and geomagnetic activity besides height).
* Dynamic model.

 In order to make it easy to remember and to exchange between models, they were respectively named as:

* **x**mowei
* **x**smade
* **x**sdamo

where **x** can be “j” for Jacchia 70, “r” for  Robert’s analytical version of Jacchia 70, “i” for Jacchia 77 model (numerical integration), “a” for the  analytical version from LaFontaine and Hughes of Jacchia’s 77, and “p” for the Lineberry  polynomial analytical version of Jacchia’s 77\. “mowei” stands for MOecular WEIght, since the static profile produces mainly the molecular weight as function of height. “smade” means Static Model for Atmospheric DEnsity, and finally “sdamo” is an acronym for Static and Dynamic Atmospheric MOdel, as seem in Table 1.

### Table 1 - Profile, Static and Dynamic functions

<table style="text-align: left; width: 800px;" border="1" cellpadding="2" cellspacing="2">
<tbody>
<tr>
<td style="vertical-align: top;"><br>
</td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;"><span style="font-weight: bold;">Static
profile</span></span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;"><span style="font-weight: bold;">Static
density</span></span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;"><span style="font-weight: bold;">Dynamic
model</span></span></td>
</tr>
<tr>
<td style="vertical-align: top;"><span style="font-family: ms sans serif;">Jacchia 70</span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;">jmowei</span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;"><span style="font-weight: bold;"></span>jsmade</span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;"><span style="font-weight: bold;"></span>jsdamo</span></td>
</tr>
<tr>
<td style="vertical-align: top;"><span style="font-family: ms sans serif;">Jacchia-Roberts</span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;">rmowei</span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;"><span style="font-weight: bold;"></span>rsmade</span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;">rsdamo</span></td>
</tr>
<tr>
<td style="vertical-align: top;"><span style="font-family: ms sans serif;">Jacchia 77</span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;">imowei</span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;">ismade</span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;">isdamo</span></td>
</tr>
<tr>
<td style="vertical-align: top;"><span style="font-family: ms sans serif;">LaFontaine-Hughes</span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;">amowei</span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;">asmade</span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;">asdamo</span></td>
</tr>
<tr>
<td style="vertical-align: top;"><span style="font-family: ms sans serif;">Jacchia-Lineberry (77 model)</span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;">pmowei</span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;">psmade</span></td>
<td style="vertical-align: top; text-align: center;"><span style="font-family: ms sans serif;">psdamo</span></td>
</tr>
</tbody>
</table>


### References

 * The orignal source of this information: http://www.dem.inpe.br/~val/atmod/default.html

 * Bowman, B. R.; Tobiska, W. K.; Marcos, F. A.; Valladares, C. The JB2006 empirical thermospheric density model. In:  Journal of Atmospheric and Solar-Terrestrial Physics, 70, pp. 774-793, 2008a. Available in: [http://sol.spacenvironment.net/~JB2006/pubs/JASTP_Bowman_JB2006_2008.pdf](http://sol.spacenvironment.net/%7EJB2006/pubs/JASTP_Bowman_JB2006_2008.pdf)

 * Bowman, B. R.; Tobiska, W. K. Marcos, F. A.; Huang, C. Y.; Lin, C. S.; Burke, W. J.,  A New Empirical Thermospheric Density Model JB2008 Using New Solar and Geomagnetic Indices.  AIAA/AAS Astrodynamics Specialist Conference.  Proceedings. Honolulu, Aug. 2008b. Available in: [http://sol.spacenvironment.net/~JB2008/pubs/AIAA_2008-6438_JB2008_Model.pdf](http://sol.spacenvironment.net/%7EJB2008/pubs/AIAA_2008-6438_JB2008_Model.pdf)

 * Carrara, V. [Implementações de modelos atmosféricos para uso em propagadores de órbita e atitude](http://www2.dem.inpe.br/val/publicacoes/carrara_inpe_5094_rpi_231.pdf). S. J. Campos, INPE, maio 1990 (INPE-5094-RPI/231).

 * CIRA 1986,  Part I: Thermosphere Model, D. Rees (ed.), Adv. Space Res. Vol. 8, n. 5, n. 6, 1988.

 * DRAO – Dominion Radio Astrophysical Observatory – Herzberg Institute of Astrophysics.  Solar Radio Monitoring Programme. [http://www.drao-ofr.hia-iha.nrc-cnrc.gc.ca](http://www.drao-ofr.hia-iha.nrc-cnrc.gc.ca/)  2011.

 * ISGI – International Service of Geomagnetic Indices. Indices: Data [http://isgi.cetp.ipsl.fr/](http://isgi.cetp.ipsl.fr/) . Access in 2011.

 * Jacchia, L. G. [ Static Diffusion Models of the Upper Atmosphere With Empirical Temperature Profiles](http://adsabs.harvard.edu/full/1965SCoA....8..215J). Cambridge, MA, SAO, 1964 (SAO Special Report 170).

 * Jacchia, L. G. Atmospheric models in the region from 110 to 2000 km. In:  COMMITTEE ON SPACE RESEARCH (COSPAR). CIRA 1972. Berlim, AkademicVerlag, 1972\. Part 3, p. 227-338.

 * Jacchia, L. G. [Thermospheric Temperature, Density and Compostion: New Models](http://adsabs.harvard.edu/full/1977SAOSR.375.....J). Cambridge,  Ma,  SAO, 1977\. (SAO Special Report No 375).

 * Lafontaine, J.; Hughes, P. [An Analytic Version of Jacchia's 1977 model atmosphere](http://www.springerlink.com/content/vp15605350n07x56/).  Celestial Mechanics, 29(3-26) 1983.

 * Mattos, B. S. Modelamento da Atmosfera Superior da Terra. In: II COLÓQUIO DE DINÂMICA ORBITAL E RESSONÂNCIA GRAVITACIONAL, 1984. ATIBAIA-SP.

 * Mueller, A. C. [ Jacchia-Lineberry upper atmosphere density model](http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19830012203_1983012203.pdf). Huston NASA, 1982\. (NASA-CR-167824).

 * Roberts Jr., C. E. [An analytical model for upper atmospheric densities based upon Jacchia's 1970 models](http://www.springerlink.com/content/l431jwg273w08u8g/).  Celestial Mechanics, 4(3/4):368-377, Dec. 1971.

 * SPDIR - Space Physics Interactive Data Resources - National Geophysical Data Center.  Geomagnetic Indices Data Set. [http://spidr.ngdc.noaa.gov/spidr/](http://spidr.ngdc.noaa.gov/spidr/)  2011.

