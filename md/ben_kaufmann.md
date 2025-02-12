# GNSS Reception at different latitudes


In your mail dating February 4, 2025, you wrote:

>Looks actually pretty good.
>
>The French have done now a test campaign in Canada and Groenland. Regarding the coverage of use of the rover, are there geographical limitations in the use of this equipment ?
>
>We on our side identified (LL from the French) a limitation for the coverage of some GPS based equipments of the aircraft.

The coverage and quality of a GNSS system depends on latitude, higher (absolute) latitudes influence the local view of the satellites. A satellite of a GNSS reaches as maximum latitude the inclination angle of its orbit. For GPS, Galileo and Beidou, the inclination angle[^1] is 55 degrees, while Glonass uses 64 degrees since the Russian territory covers higher latitudes.

This effects the Dilution of Precision (DOP) parameters[^2] of the constellation (or their combination). The DOP value is to be interpreted as the factor by which the GNSS system measurement error is multiplied to obtain the position precision.
$$\sigma_{\text{pos}} = \text{PDOP} \cdot \sigma_{\text{GNSS}}$$
So the lower the DOP value, the better position estimates. For calculating a position a minimal of 4 satellites are needed, the best observed geometry comprises 1 satellite near the zenithal position of the observer, and 3 satellites at elevation angles of about 15 degrees above horizon separated by about 120 degrees in azimuth. Usually DOP values greater than 6 are considered to be too high for reliable position estimation. 

For the rover we use standard only GPS and Galileo. For a location with (absolute) latitude below the inclination angle, satellites will be visible near the zenithal direction allowing to near the optimal geometry. When the (absolute) latitude exceeds the inclination angle, the satellites will be visible at lower elevation angles, which reduces the DOP value.


[^1]: The inclination angle is the angle between the orbital plane of the satellite and the equatorial plane of the Earth.
[^2]: The DOP is a measure for the geometry of the constellation. It is a function of the relative positions of the satellites and the receiver. There are 5 DOP parameters, which are the GDOP, PDOP, HDOP, VDOP and TDOP. PDOP is used for evaluating the 3D position, HDOP for the horizontal position, VDOP for the vertical position and TDOP for the time. The effect of the (absolute) latitude is, due to the high number of satellites in each individual constellations moderate. To demonstrate this, we have conducted a simulation with the following constellation:

- GPS + Galileo
- GPS + Galileo + Glonass

for different latitudes:

- RMA $\varphi=50^d50^m, \lambda=4^d23^m$ 
- Fairbanks $\varphi=64^d50^m, \lambda=-147^d42^m$ 

The  figures 1 and 2 show the DOP values and the total number of visible satellites for respectively RMA and Fairbanks. The skyview shows the satellites visible for February 7, 2025. The central point corresponds to the location of the observer, the outer circle is your horizon while the inner circles are the elevation angles of the satellites.

![Satellite view at RMA](./RMA-gps-ops-galileo-20250207-skyview.png){width=95%}

![Satellite view at fairbanks](./Fairbanks-gps-ops-galileo-20250207-skyview.png){width=95%}

As you can see, at RMA we have Galileo and GPS satellites overhead (at our zenith) while at Fairbanks we have none satellites overhead. Since the inclination angle of Glonass is higher, the following figure 3 shows you the Fairbanks view for GPA, Galileo and Glonass.

![Satellite view at fairbanks including Glonass](./Fairbanks-gps-ops-galileo-glo-ops-20250207-skyview.png){width=95%}

We see that Glonass satellites do reach the zenith of Fairbanks.

\clearpage

The following figures (4-6) show the number of visible satellites (on right y-axis) and the DOP values (on left y-axis) for respectively RMA and Fairbanks,  plot 6 adds the Glonass constellation.

![DOP values at RMA](./RMA-gps-ops-galileo-20250207-DOP.png){width=75%}

![DOP values at fairbanks](./Fairbanks-gps-ops-galileo-20250207-DOP.png){width=75%}

![DOP values at fairbanks including Glonass](./Fairbanks-gps-ops-galileo-glo-ops-20250207-DOP.png){width=75%}

Between RMA and Fairbanks, the DOP values are slightly higher at Fairbanks, but still very good (below 2 overall). This is due to not having overhead satellites. The last figure shows at the DOP values at Fairbanks including Glonass with DOP values still remaining very similar.

\clearpage

### Conclusion

Though the GNSS reception at different latitudes is influenced by the inclination angle of the satellites, with the current number of satellites available in each constellation, the DOP values remain very good. So I do not see the problem the French had while operating at higher latitudes.

If needed I can make similar plots for whatever location on Earth and for any GNSS constellation alone or combined.

