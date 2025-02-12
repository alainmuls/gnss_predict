#!/usr/bin/env python
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from rich import print

from gnss_tracking.GNSS import gpstime
from gnss_tracking.output_files import (
    create_dop_file,
    create_geodetic_file,
    create_tle_arc_file,
    create_visible_file,
)
from gnss_tracking.plotting import (
    plot_dop_visible,
    plot_groundtracks,
    plot_sky_view,
    plot_visible_sats,
)
from gnss_tracking.station import Station
from gnss_tracking.tle.tle_class import TLEManager
from gnss_tracking.utils.cli import parse_arguments

matplotlib.rcParams["backend"] = "TkAgg"

__author__ = "amuls"

# exit codes
E_SUCCESS = 0
E_FILE_NOT_EXIST = 1
E_NOT_IN_PATH = 2
E_UNKNOWN_OPTION = 3
E_TIME_PASSED = 4
E_WRONG_OPTION = 5
E_SIGNALTYPE_MISMATCH = 6
E_DIR_NOT_EXIST = 7
E_TIMING_ERROR = 8
E_REQUEST_ERROR = 9
E_NO_TLE_FILE = 10
E_FAILURE = 99


# main starts here
def main():
    # matplotlib.use("TkAgg")
    # matplotlib.use("Agg")

    args = parse_arguments()

    # if ',' in excludeSats:
    #     sats_to_exclude = excludeSats.split(',')
    if len(args.exclude) > 0:
        sats_to_exclude = args.exclude.split(",")
    else:
        sats_to_exclude = None

    # import tle data from NORAD if internet_on(), save as sat=ephem.readtle(...)-----
    tle_file = TLEManager.get_tle_from_norad(args.gnss, args.verbose)

    # Initialize Station class defaults
    Station.initialize_defaults()

    # read in the observer info (name, latitude, longitude, date
    marker_info = Station.set_observer_data(args.observer, args.date, args.verbose)
    if args.verbose:
        print(marker_info)

    # read in the list of satellites from the TLE
    tle_satellite_list = TLEManager.load_tle(tle_file, args.verbose)

    # calculate the interval settings for a full day prediction starting at 00:00:00 hr of predDate
    prediction_dt, prediction_count = marker_info.set_observation_times(
        args.start, args.end, args.interval, args.verbose
    )

    # calculate the information for each SVs in tle_satellite_list
    sub_latitude = np.empty([prediction_count, np.size(tle_satellite_list)])
    sub_longitude = np.empty([prediction_count, np.size(tle_satellite_list)])
    azim = np.empty([prediction_count, np.size(tle_satellite_list)])
    elev = np.empty([prediction_count, np.size(tle_satellite_list)])
    dist = np.empty([prediction_count, np.size(tle_satellite_list)])
    dist_velocity = np.empty([prediction_count, np.size(tle_satellite_list)])
    eclipsed = np.empty([prediction_count, np.size(tle_satellite_list)])
    xDOP = np.empty([prediction_count, 4])  # order is HDOP, VDOP, TDOP and NrOfSVsUsed
    nrExcluded = np.empty(prediction_count)

    for i, dt in enumerate(prediction_dt):

        marker_info.date = dt
        elevTxt = ""
        for j, sat in enumerate(tle_satellite_list):
            sat.compute(marker_info)
            sub_latitude[i, j] = np.rad2deg(sat.sublat)
            sub_longitude[i, j] = np.rad2deg(sat.sublong)
            azim[i, j] = np.rad2deg(sat.az)
            elev[i, j] = np.rad2deg(sat.alt)
            dist[i, j] = sat.range
            dist_velocity[i, j] = sat.range_velocity
            eclipsed[i, j] = sat.eclipsed

            # elevTxt += "%6.1f  " % elev[i][j]

        # determine the visible satellites at this instance
        indexVisSats = []
        indexVisSats = np.where(elev[i, :] >= args.cutoff)
        indexVisSatsUsed = indexVisSats[0]
        index2Delete = []

        # exclude the SVs which have no valid signal (cfr sats_to_exclude list)
        if sats_to_exclude is not None:
            for k, prn in enumerate(tle_satellite_list):
                if k in indexVisSats[0]:
                    for jj, prnX in enumerate(sats_to_exclude):
                        if prnX in prn.name:
                            index2Delete.append(np.where(indexVisSats[0] == k)[0])

            index2Delete.sort(reverse=True)
            nrExcluded[i] = np.size(index2Delete)

            for k in index2Delete:
                indexVisSatsUsed = np.delete(indexVisSatsUsed, k)

        # calculate xDOP values when at least 4 sats are visible above cutoff angle
        if np.size(indexVisSatsUsed) >= 4:
            A = np.matrix(np.empty([np.size(indexVisSatsUsed), 4], dtype=float))
            elevation_visible_sats_rad = np.radians(elev[i, indexVisSatsUsed])
            azimuth_visible_sat_rad = np.radians(azim[i, indexVisSatsUsed])

            for j in range(np.size(indexVisSatsUsed)):
                A[j, 0] = np.cos(azimuth_visible_sat_rad[j]) * np.cos(
                    elevation_visible_sats_rad[j]
                )
                A[j, 1] = np.sin(azimuth_visible_sat_rad[j]) * np.cos(
                    elevation_visible_sats_rad[j]
                )
                A[j, 2] = np.sin(elevation_visible_sats_rad[j])
                A[j, 3] = 1.0

            # calculate ATAInv en get the respective xDOP parameters (HDOP, VDOP and TDOP)
            AT = A.getT()
            ATA = AT * A
            ATAInv = ATA.getI()

            xDOP[i, 0] = np.sqrt(ATAInv[0, 0] + ATAInv[1, 1])  # HDOP
            xDOP[i, 1] = np.sqrt(ATAInv[2, 2])  # VDOP
            xDOP[i, 2] = np.sqrt(ATAInv[3, 3])  # TDOP
            xDOP[i, 3] = np.size(indexVisSatsUsed)
        else:  # not enough visible satellites
            xDOP[i] = [np.nan, np.nan, np.nan, np.nan]

    # set all elev < cutoff to NAN
    elev[elev < args.cutoff] = np.nan

    # write to results file
    create_visible_file(
        observer=marker_info,
        list_satellites=tle_satellite_list,
        prediction_dt=prediction_dt,
        elevation=elev,
        azimuth=azim,
        excluded_satellites=sats_to_exclude,
        cli_args=args,
    )

    create_dop_file(
        observer=marker_info,
        prediction_dt=prediction_dt,
        elevation=elev,
        xDOPs=xDOP,
        cli_args=args,
    )

    create_geodetic_file(
        observer=marker_info,
        list_satellites=tle_satellite_list,
        prediction_dt=prediction_dt,
        latitudes=sub_latitude,
        longitudes=sub_longitude,
        cli_args=args,
    )

    create_tle_arc_file(
        observer=marker_info,
        list_satellites=tle_satellite_list,
        prediction_dt=prediction_dt,
        elevation=elev,
        cli_args=args,
    )

    # create plots
    plot_visible_sats(
        observer=marker_info,
        list_satellites=tle_satellite_list,
        prediction_dt=prediction_dt,
        elevations=elev,
        excluded_satellites=sats_to_exclude,
        cli_args=args,
    )

    plot_sky_view(
        observer=marker_info,
        list_satellites=tle_satellite_list,
        prediction_dt=prediction_dt,
        elevations=elev,
        azimuths=azim,
        excluded_satellites=sats_to_exclude,
        cli_args=args,
    )

    plot_groundtracks(
        observer=marker_info,
        list_satellites=tle_satellite_list,
        prediction_dt=prediction_dt,
        sat_latitudes=sub_latitude,
        sat_longitudes=sub_longitude,
        excluded_satellites=sats_to_exclude,
        cli_args=args,
    )

    plot_dop_visible(
        observer=marker_info,
        list_satellites=tle_satellite_list,
        prediction_dt=prediction_dt,
        elevations=elev,
        xDOPs=xDOP,
        excluded_satellites=sats_to_exclude,
        cli_args=args,
    )

    # show all plots
    if args.verbose:
        plt.show()

    # end program
    return E_SUCCESS


if __name__ == "__main__":
    main()
