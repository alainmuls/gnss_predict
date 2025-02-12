import argparse
import csv
import math
from datetime import datetime
from typing import List

import ephem
import numpy as np
from ephem import EarthSatellite

from gnss_tracking.GNSS import gpstime

from .station import Station


def format_satellite_header(
    page_satellites: List[EarthSatellite], include_vis_column: bool = False
) -> str:
    """
    Format satellite names into header lines with optional visibility column.

    Args:
        page_satellites: List of satellites for current page
        include_vis_column: Whether to include the #Vis column in header

    Returns:
        Formatted header string with satellite names
    """
    satellite_line1 = ""
    satellite_line2 = ""

    for i, sat in enumerate(page_satellites):
        if len(sat.name) < 12:
            if include_vis_column:
                satellite_line1 += f"  {sat.name:11s}"
            else:
                satellite_line1 += f"   {sat.name:12s}"
        else:
            if include_vis_column:
                satellite_line1 += f"  {sat.name[:10]:11s}"
                endChar = min(20, len(sat.name))
                satellite_line2 += f"  {sat.name[10:endChar]:11s}"
            else:
                if i == 0:
                    satellite_line1 = "   "
                    satellite_line2 = "  "
                satellite_line1 += f"  {sat.name[:10]:12s}"
                endChar = min(20, len(sat.name))
                satellite_line2 += f"   {sat.name[10:endChar]:11s}"

    # Build complete header
    header = ""
    if include_vis_column:
        header += f"      |#Vis|{satellite_line1}\n"
        if satellite_line2:
            header += f"            {satellite_line2}\n"
    else:
        header += f"      {satellite_line1}\n"
        if satellite_line2:
            header += f"      {satellite_line2}\n"

    header += "\n"

    return header


def create_dop_file(
    observer: Station,
    prediction_dt: List[datetime],
    elevation: np.ndarray,
    xDOPs: np.ndarray,
    cli_args: argparse.Namespace,
):
    """
    Creates a DOP (Dilution of Precision) file for a given observer and prediction data.

    Parameters:
    observer (Station): The observer station containing information such as name, latitude, and longitude.
    prediction_dt (List[datetime]): List of datetime objects representing the prediction times.
    elev (np.ndarray): 2D numpy array containing elevation data of satellites.
    xDOPs (np.ndarray): 2D numpy array containing DOP values (HDOP, VDOP, PDOP, TDOP, GDOP).
    cli_args (argparse.Namespace, optional): Command line arguments, if any. Defaults to None.

    Returns:
    None

    Raises:
    IOError: If there is an issue accessing or writing to the file.

    """

    xdop_ofn = f"/tmp/{observer.name}_{cli_args.gnss.replace(',', '_')}_{observer.get_ymd_str()}_DOP.txt"
    if cli_args.verbose:
        print(f"\tCreating DOP file: {xdop_ofn}")
    try:
        ofd = open(xdop_ofn, "w")
        # write the observer info out
        ofd.write(f"Observer: {observer.name}\n")
        ofd.write(f"     lat: {ephem.degrees(observer.lat)}\n")
        ofd.write(f"     lon: {ephem.degrees(observer.lon)}\n")
        ofd.write(
            f"    date: {observer.get_ymd()[0]:04d}/{observer.get_ymd()[1]:02d}/{observer.get_ymd()[2]:02d}\n"
        )
        ofd.write(f"  cutoff: {cli_args.cutoff:2d}\n\n")

        ofd.write("      |#Used/#Vis|   HDOP   VDOP   PDOP   TDOP   GDOP\n\n")

        for i, predDate in enumerate(prediction_dt):
            ofd.write(f"{predDate.hour:02d}:{predDate.minute:02d}")

            # number of visible satellites
            if ~np.isnan(xDOPs[i, 3]):
                ofd.write(
                    f" | {xDOPs[i, 3]:3.0f} / {np.count_nonzero(~np.isnan(elevation[i, :])):2d} |"
                )
            else:
                ofd.write(
                    f" |  -- / {np.count_nonzero(~np.isnan(elevation[i, :])):2d} |"
                )

            # write the DOP values in order
            if ~np.isnan(xDOPs[i, 0]):
                PDOP2 = xDOPs[i, 0] * xDOPs[i, 0] + xDOPs[i, 1] * xDOPs[i, 1]
                ofd.write(
                    f" {xDOPs[i, 0]:6.1f} {xDOPs[i, 1]:6.1f} {np.sqrt(PDOP2):6.1f} "
                    f"{xDOPs[i, 2]:6.1f} {np.sqrt(PDOP2 + xDOPs[i, 2] * xDOPs[i, 2]):6.1f}"
                )
            else:
                ofd.write(" ------ ------ ------ ------ ------")
            ofd.write("\n")

        # close the file
        ofd.close()
    except IOError:
        print("  Access to file %s failed" % xdop_ofn)


def create_geodetic_file(
    observer: Station,
    list_satellites: List[EarthSatellite],
    prediction_dt: List[datetime],
    latitudes: np.ndarray,
    longitudes: np.ndarray,
    cli_args: argparse.Namespace,
):
    """
    Create a geodetic file with satellite tracking information, split into pages of 12 satellites each.
    """
    geodetic_ofn = f"/tmp/{observer.name}_{cli_args.gnss.replace(',', '_')}_{observer.get_ymd_str()}_GEOD.txt"
    if cli_args.verbose:
        print(f"\tCreating sub-stellar file: {geodetic_ofn}")

    SATS_PER_PAGE = 12
    total_pages = (len(list_satellites) + SATS_PER_PAGE - 1) // SATS_PER_PAGE

    try:
        with open(geodetic_ofn, "w") as ofd:
            # Write header information once
            ofd.write(f"Observer: {observer.name}\n")
            ofd.write(f"     lat: {ephem.degrees(observer.lat)}\n")
            ofd.write(f"     lon: {ephem.degrees(observer.lon)}\n")
            ofd.write(f"    date: {observer.get_ymd_str()}\n\n")

            # Process each page of satellites
            for page in range(total_pages):
                start_idx = page * SATS_PER_PAGE
                end_idx = min(start_idx + SATS_PER_PAGE, len(list_satellites))
                page_satellites = list_satellites[start_idx:end_idx]

                ofd.write(
                    format_satellite_header(page_satellites, include_vis_column=False)
                )

                # Write position data for this page's satellites
                for i, predDate in enumerate(prediction_dt):
                    ofd.write(f"{predDate.hour:02d}:{predDate.minute:02d} |")

                    # Write lat/lon for satellites on this page
                    for j in range(start_idx, end_idx):
                        ofd.write(f"  {latitudes[i, j]:5.1f} {longitudes[i, j]:6.1f}")
                    ofd.write("\n")

                # Add page separator if not the last page
                if page < total_pages - 1:
                    ofd.write("\n\n\n")

    except IOError:
        print(f"\tAccess to file {geodetic_ofn} failed")


def create_visible_file(
    observer: Station,
    list_satellites: List[EarthSatellite],
    prediction_dt: List[datetime],
    elevation: np.ndarray,
    azimuth: np.ndarray,
    cli_args: argparse.Namespace,
    excluded_satellites: List[str] = None,
):
    """
    Writes visibility information for a list of satellites to a file, organized in pages of 12 satellites each.
    """
    visible_ofn = f"/tmp/{observer.name}_{cli_args.gnss.replace(',', '_')}_{observer.get_ymd_str()}_VIS.txt"
    if cli_args.verbose:
        print(f"\tCreating visibility file: {visible_ofn}")

    SATS_PER_PAGE = 12
    total_pages = (len(list_satellites) + SATS_PER_PAGE - 1) // SATS_PER_PAGE

    # Mark excluded satellites
    indexIncluded = np.ones(np.size(list_satellites), dtype=bool)
    if excluded_satellites is not None:
        for i, sat in enumerate(list_satellites):
            for PRN in excluded_satellites:
                if PRN in sat.name:
                    indexIncluded[i] = False

    try:
        with open(visible_ofn, "w") as ofd:
            # Write header information once
            ofd.write(f"Observer: {observer.name}\n")
            ofd.write(f"     lat: {ephem.degrees(observer.lat)}\n")
            ofd.write(f"     lon: {ephem.degrees(observer.lon)}\n")
            ofd.write(f"    date: {observer.get_ymd_str()}\n")
            ofd.write(f"  cutoff: {cli_args.cutoff:2d}\n\n")

            # Process each page of satellites
            for page in range(total_pages):
                start_idx = page * SATS_PER_PAGE
                end_idx = min(start_idx + SATS_PER_PAGE, len(list_satellites))
                page_satellites = list_satellites[start_idx:end_idx]

                ofd.write(
                    format_satellite_header(page_satellites, include_vis_column=True)
                )

                # Write visibility data for this page's satellites
                for i, predDate in enumerate(prediction_dt):
                    ofd.write(f"{predDate.hour:02d}:{predDate.minute:02d}")
                    ofd.write(f" | {np.count_nonzero(~np.isnan(elevation[i, :])):2d} |")

                    for j in range(start_idx, end_idx):
                        if indexIncluded[j]:
                            if math.isnan(elevation[i, j]):
                                ofd.write("  ---- ----- ")
                            else:
                                ofd.write(
                                    f"  {elevation[i, j]:4.1f} {azimuth[i, j]:5.1f} "
                                )
                        else:
                            if math.isnan(elevation[i, j]):
                                ofd.write(" (---- -----)")
                            else:
                                ofd.write(
                                    f" ({elevation[i, j]:4.1f} {azimuth[i, j]:5.1f})"
                                )
                    ofd.write("\n")

                # Add page separator if not the last page
                if page < total_pages - 1:
                    ofd.write("\n\n\n")

    except IOError:
        print(f"\tAccess to file {visible_ofn} failed")


def create_tle_arc_file(
    observer: Station,
    list_satellites: List[EarthSatellite],
    prediction_dt: List[datetime],
    elevation: np.ndarray,
    cli_args: argparse.Namespace,
):
    """
    Creates arc files for multiple GNSS systems based on TLE predictions.

    Args:
        gnss_systems: Comma-separated string of GNSS systems (e.g. "gps-ops,galileo")
        list_satellites: List of satellites to process
        prediction_dt: List of prediction datetimes
        elevation: Satellite elevation array
        verbose: Enable verbose output
    """
    # Get first prediction date for filename
    base_date = prediction_dt[0]

    # Convert prediction times to GPS time
    pred_time_seconds = []
    for pred_date in prediction_dt:
        pyUTC = gpstime.mkUTC(
            pred_date.year,
            pred_date.month,
            pred_date.day,
            pred_date.hour,
            pred_date.minute,
            pred_date.second,
        )
        WkNr, TOW = gpstime.wtFromUTCpy(pyUTC, leapSecs=0)
        pred_time_seconds.append(TOW)

    pred_time_seconds = np.array(pred_time_seconds)

    # Create arc file for each GNSS system
    for gnss in cli_args.gnss.split(","):
        arc_ofn = f"/tmp/{observer.name}_{gnss}_{observer.get_ymd_str()}_TLE_arc.csv"

        if cli_args.verbose:
            print(f"\tCreating arc file for {gnss}: {arc_ofn}")

        with open(arc_ofn, "w") as arc_ofd:
            arc_writer = csv.writer(arc_ofd, delimiter=",", lineterminator="\n")

            # Write header row
            header = [
                "PRN",
                "GPS_Week",
                "Start_TOW",
                "End_TOW",
                "Start_Time",
                "End_Time",
                "Duration",
            ]
            arc_writer.writerow(header)

            # Process satellites belonging to current GNSS
            for i, sat in enumerate(list_satellites):
                if gnss[:3].lower() in sat.name.lower():
                    # Transform satellite name to PRN number
                    # try:
                    #     prn_nr = int(sat.name[-3:-1]) + 70
                    # except ValueError:
                    prn_nr = sat.name

                    # Find elevation arcs looking for non NaN values
                    pos_elev_indices = np.where(np.isfinite(elevation[:, i]))[0]

                    # if satellite is visible, determine arcs based on continuous elevation
                    if np.size(pos_elev_indices) > 0:
                        # Add -10 index to catch first arc
                        pos_elev_indices_extend = np.append(
                            np.array([-10]), pos_elev_indices
                        )
                        new_arc_indices = np.where(
                            np.ediff1d(pos_elev_indices_extend) > 1
                        )[0]

                        # Calculate arc times
                        arc_start_times = []
                        arc_end_times = []

                        for j, new_arc in enumerate(new_arc_indices):
                            arc_start_times.append(
                                pred_time_seconds[pos_elev_indices[new_arc]]
                            )

                            if j < len(new_arc_indices) - 1:
                                arc_end_times.append(
                                    pred_time_seconds[
                                        pos_elev_indices[new_arc_indices[j + 1] - 1]
                                    ]
                                )
                            else:
                                arc_end_times.append(
                                    pred_time_seconds[pos_elev_indices[-1]]
                                )

                        # Write arcs to file
                        for t in range(len(arc_start_times)):
                            _, _, _, start_hh, start_mm, _ = gpstime.UTCFromGps(
                                WkNr, arc_start_times[t], leapSecs=0
                            )
                            start_hhmm = f"{start_hh:02d}:{start_mm:02d}"

                            _, _, _, end_hh, end_mm, _ = gpstime.UTCFromGps(
                                WkNr, arc_end_times[t], leapSecs=0
                            )
                            end_hhmm = f"{end_hh:02d}:{end_mm:02d}"

                            # calculate the duration and express in hours
                            arc_duration_sec = arc_end_times[t] - arc_start_times[t]
                            # Convert to hours and minutes
                            duration_hours = int(arc_duration_sec // 3600)
                            duration_minutes = int((arc_duration_sec % 3600) // 60)

                            # Format as HH:MM
                            arc_duration_hhmm = (
                                f"{duration_hours:02d}:{duration_minutes:02d}"
                            )

                            arc_row = [
                                prn_nr,
                                WkNr,
                                arc_start_times[t],
                                arc_end_times[t],
                                start_hhmm,
                                end_hhmm,
                                arc_duration_hhmm,
                            ]

                            # write to CSV file
                            arc_writer.writerow(arc_row)
