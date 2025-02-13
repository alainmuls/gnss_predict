import argparse
import copy
from datetime import datetime
from typing import List

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import ephem
import matplotlib
import matplotlib.cm as cm
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np

from ephem import EarthSatellite
from matplotlib.pyplot import figure, grid, rc, rcParams

from gnss_predict.config import rich_console

from .GNSS import gpstime
from .station import Station

matplotlib.use("TkAgg")


def plot_visible_sats(
    observer: Station,
    list_satellites: List[EarthSatellite],
    prediction_dt: List[datetime],
    elevations: np.ndarray,
    cli_args: argparse.Namespace,
    excluded_satellites: List[str] = None,
):
    """
    Plots the visibility of satellites from a given observation station.

    Parameters:
    observer (Station): The observation station.
    list_satellites (List[EarthSatellite]): List of EarthSatellite objects representing the satellites to be plotted.
    prediction_dt (List[datetime]): List of datetime objects representing the prediction times.
    elevations (np.ndarray): 2D array of elevations angles for the satellites.
    cli_args (argparse.Namespace): Command line arguments.
    excluded_satellites (List[str], optional): List of satellite names to be excluded from the plot. Defaults to None.

    Returns:
    None
    """
    with rich_console.status("Plotting satellite visibility", spinner="aesthetic"):

        plt.style.use("ggplot")

        # create a figure
        fig = plt.figure(figsize=(14.0, 9.0))
        ax1 = plt.gca()

        colors = iter(cm.jet(np.linspace(0, 1, len(list_satellites))))

        # local copy to work with
        elev2 = copy.deepcopy(elevations)

        # plot the lines for visible satellites
        for i, sat in enumerate(list_satellites):
            elev2[~np.isnan(elev2)] = i + 1  # create a horizontal line @ height i+1
            line_style = "-"
            if excluded_satellites is not None:
                for j, PRN in enumerate(excluded_satellites):
                    print(
                        f"DEBUG:   PRN[{j} of {np.size(excluded_satellites)}] = {PRN}"
                    )
                    if PRN in sat.name:
                        line_style = "--"

            if len(list_satellites) < 50:
                plt.plot(
                    prediction_dt,
                    elev2[:, i],
                    linewidth=5,
                    color=next(colors),
                    linestyle=line_style,
                    label=sat.name,
                )
            else:
                plt.plot(
                    prediction_dt,
                    elev2[:, i],
                    linewidth=5,
                    color=next(colors),
                    linestyle=line_style,
                )

        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax1.xaxis.set_major_locator(
            mdates.MinuteLocator(byminute=0, interval=1, tz=None)
        )

        # create array for setting the satellite labels on y axis
        satNames = []
        for i in range(len(list_satellites)):
            satNames.append(list_satellites[i].name)

        # set the tick marks
        plt.xticks(rotation=50, size="medium")
        plt.yticks(range(1, len(list_satellites) + 1), satNames, size="x-small")

        # color the sat labels ticks
        colors2 = iter(cm.jet(np.linspace(0, 1, len(list_satellites))))
        for i, tl in enumerate(ax1.get_yticklabels()):
            tl.set_color(next(colors2))

        plt.grid(True)
        # set the limits for the y-axis
        plt.ylim(0, len(list_satellites) + 2)
        # ax1.set_xlabel('Time of Day', fontsize='x-large')
        # plot title
        plt.title(
            f"{cli_args.gnss.replace(',', ' & ').upper()} Satellite Visibility",
            fontsize="x-large",
        )
        yyyy, mm, dd = observer.get_ymd()
        annotateTxt = (
            f"Station: {observer.name} @ (ϕ {ephem.degrees(observer.lat)}, λ {ephem.degrees(observer.lon)}) - "
            f"Date {yyyy:04d}/{mm:02d}/{dd:02d} - Cutoff {cli_args.cutoff:2d}"
        )
        plt.text(
            0.5,
            0.99,
            annotateTxt,
            horizontalalignment="center",
            verticalalignment="top",
            transform=ax1.transAxes,
            fontsize="medium",
        )

        plot_ofn = f"/tmp/{observer.name}_{cli_args.gnss.replace(',', '_')}_{observer.get_ymd_str()}_visibility.png"
        fig.savefig(plot_ofn, dpi=fig.dpi)
        if cli_args.verbose:
            print(f"\tSaved figure {plot_ofn}")

        if cli_args.verbose:
            plt.draw()


def plot_sky_view(
    observer: Station,
    list_satellites: List[EarthSatellite],
    prediction_dt: List[datetime],
    elevations: np.ndarray,
    azimuths: np.ndarray,
    cli_args: argparse.Namespace,
    excluded_satellites: List[str] = None,
):
    """
    Plots a sky view from a given observation station.

    Parameters:
        observer (Station): The observation station.
        list_satellites (List[EarthSatellite]): List of EarthSatellite objects representing the satellites to be plotted.
        prediction_dt (List[datetime]): List of datetime objects representing the prediction times.
        elevations (np.ndarray): 2D array of elevations angles for the satellites.
        azimuths (np.ndarray): 2D array of azimuths angles for the satellites.
        cli_args (argparse.Namespace): Command line arguments.
        excluded_satellites (List[str], optional): List of satellite names to be excluded from the plot. Defaults to None.

    Returns:
    None
    """
    with rich_console.status("Plotting satellite visibility", spinner="aesthetic"):

        plt.style.use("ggplot")
        # rc('grid', color='#999999', linewidth=1, linestyle='-', alpha=[0].6)
        rc("xtick", labelsize="x-small")
        rc("ytick", labelsize="x-small")

        # force square figure and square axes looks better for polar, IMO
        width, height = rcParams["figure.figsize"]
        size = min(width, height) * 2

        # make a square figure
        fig = figure(figsize=(size, size))

        # set the axis (0 azimuth is North direction, azimuth indirect angle)
        ax = fig.add_axes([0.10, 0.15, 0.8, 0.8], projection="polar")
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)

        # Define the xticks
        ax.set_xticks(np.linspace(0, 2 * np.pi, 13))
        xLabel = [
            "N",
            "30",
            "60",
            "E",
            "120",
            "150",
            "S",
            "210",
            "240",
            "W",
            "300",
            "330",
            "??",
        ]
        ax.set_xticklabels(xLabel)

        # Define the yticks
        ax.set_yticks(np.linspace(0, 90, 7))
        y_labels = ["", "75", "60", "45", "30", "15", ""]
        ax.set_yticklabels(y_labels)

        # draw a grid
        grid(True)

        # plot the sky-tracks for each PRN
        colors = iter(cm.jet(np.linspace(0, 1, len(list_satellites))))
        satellite_labels = []

        # find full hours in date to set the elevations/azimaccordingly
        prediction_time_seconds = []
        hourTxt = []

        for t, predDate in enumerate(prediction_dt):
            prediction_time_seconds.append(gpstime.hms_to_seconds(predDate.time()))
        prediction_time_seconds = np.array(prediction_time_seconds)

        indexHour = np.where(np.fmod(prediction_time_seconds, 3600.0) == 0)
        hourTxt.append(prediction_time_seconds[indexHour])

        for i, prn in enumerate(list_satellites):
            satellite_labels.append("%s" % prn.name)
            satColor = next(colors)
            azimuths_rad = [np.radians(az) for az in azimuths[:, i]]
            zenithal_distances = [(90 - el) for el in elevations[:, i]]

            line_style = "-"
            if excluded_satellites is not None:
                for j, exlPRN in enumerate(excluded_satellites):
                    if exlPRN in prn.name:
                        line_style = "--"

            if len(list_satellites) < 50:
                ax.plot(
                    azimuths_rad,
                    zenithal_distances,
                    color=satColor,
                    marker=".",
                    markersize=4,
                    linestyle=line_style,
                    linewidth=1,
                    label=satellite_labels[-1],
                )
            else:
                ax.plot(
                    azimuths_rad,
                    zenithal_distances,
                    color=satColor,
                    marker=".",
                    markersize=4,
                    linestyle=line_style,
                    linewidth=1,
                )

            # annotate with the hour labels
            prn_azimuth_hours = azimuths[:, i][indexHour]
            prn_elevation_hours = elevations[:, i][indexHour]

            prn_azimuth_rad_hours = [np.radians(az + 2) for az in prn_azimuth_hours]
            prn_zenithal_distances = [(90 - el) for el in prn_elevation_hours]

            for j, hour in enumerate(hourTxt[0]):
                prn_zenithal_distance = prn_zenithal_distances[j]
                if ~np.isnan(prn_zenithal_distance):
                    prn_azimuth_rad_hour = prn_azimuth_rad_hours[j]
                    hour = int(float(hour) / 3600.0)
                    plt.text(
                        prn_azimuth_rad_hour,
                        prn_zenithal_distance,
                        hour,
                        fontsize="x-small",
                        color=satColor,
                    )

        # adjust the legend location
        if len(list_satellites) < 50:
            my_legend = ax.legend(
                bbox_to_anchor=(0.5, -0.15),
                loc="lower center",
                ncol=min(np.size(satellite_labels), 7),
                fontsize="xx-small",
                markerscale=0.6,
            )
            for legobj in my_legend.legend_handles:
                legobj.set_linewidth(5.0)

        plt.title(
            f"{cli_args.gnss.replace(',', ' & ').upper()} Satellite Visibility",
            fontsize="x-large",
            x=0.5,
            y=0.99,
            horizontalalignment="center",
        )
        yyyy, mm, dd = observer.get_ymd()
        annotateTxt = f"Station: {observer.name} @ (ϕ {ephem.degrees(observer.lat)}, λ {ephem.degrees(observer.lon)})"
        plt.text(
            -0.075,
            0.975,
            annotateTxt,
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax.transAxes,
            fontsize="medium",
        )
        annotateTxt = f"Date {yyyy:04d}/{mm:02d}/{dd:02d} - Cutoff {cli_args.cutoff:2d}"
        plt.text(
            -0.075,
            0.950,
            annotateTxt,
            horizontalalignment="left",
            verticalalignment="top",
            transform=ax.transAxes,
            fontsize="medium",
        )

        # needed for having radial axis span from 0 => 90 degrees and y-labels along north axis
        ax.set_rmax(90)
        ax.set_rmin(0)
        ax.set_rlabel_position(0)

        # ax2 = ax1.twinx()
        plot_ofn = f"/tmp/{observer.name}_{cli_args.gnss.replace(',', '_')}_{observer.get_ymd_str()}_skyview.png"
        fig.savefig(plot_ofn, dpi=fig.dpi)
        if cli_args.verbose:
            print(f"\tSaved figure {plot_ofn}")

        if cli_args.verbose:
            plt.draw()


def plot_groundtracks(
    observer: Station,
    list_satellites: List[EarthSatellite],
    prediction_dt: List[datetime],
    sat_latitudes: np.ndarray,
    sat_longitudes: np.ndarray,
    cli_args: argparse.Namespace,
    excluded_satellites: List[str] = None,
):
    """
    Plots the ground tracks of the satellites on a map using Cartopy.

    Args:
    Parameters:
        observer (Station): The observation station.
        list_satellites (List[EarthSatellite]): List of EarthSatellite objects representing the satellites to be plotted.
        prediction_dt (List[datetime]): List of datetime objects representing the prediction times.
        sat_latitudes (np.ndarray): array of latitudes for the satellites.
        sat_longitudes (np.ndarray): array of longitudes for the satellites.
        cli_args (argparse.Namespace): Command line arguments.
        excluded_satellites (List[str], optional): List of satellite names to be excluded from the plot. Defaults to None.

    Return:
        None
    """
    with rich_console.status("Plotting groundtracks", spinner="aesthetic"):
        plt.style.use("ggplot")
        fig = plt.figure(figsize=(16.0, 10.5))

        # Create map using Miller projection
        ax = plt.axes(projection=ccrs.Miller())

        # Add map features
        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS, linestyle=":")
        ax.add_feature(cfeature.LAND, alpha=0.5)
        ax.add_feature(cfeature.OCEAN)
        ax.gridlines()

        # Plot the observer location
        ax.plot(
            observer.lon / ephem.pi * 180.0,
            observer.lat / ephem.pi * 180.0,
            color="blue",
            marker="o",
            markersize=5,
            transform=ccrs.PlateCarree(),
        )

        # Add observer name
        plt.text(
            observer.lon / ephem.pi * 180.0 + 0.5,
            observer.lat / ephem.pi * 180.0 + 0.5,
            observer.name,
            fontsize="small",
            color="blue",
            transform=ccrs.PlateCarree(),
        )

        # Set colormap
        colors = iter(cm.jet(np.linspace(0, 1, len(list_satellites))))
        satellite_labels = []

        # Calculate time indices for hour labels
        prediction_time_seconds = []
        for t, predDate in enumerate(prediction_dt):
            prediction_time_seconds.append(gpstime.hms_to_seconds(predDate.time()))
        prediction_time_seconds = np.array(prediction_time_seconds)
        indexHour = np.where(np.fmod(prediction_time_seconds, 3600.0) == 0)
        hourTxt = [prediction_time_seconds[indexHour]]

        # Plot satellite tracks
        for i, SV in enumerate(list_satellites):
            satellite_labels.append(f"{SV.name}")
            satColor = next(colors)
            line_style = (
                "--"
                if excluded_satellites
                and any(PRN in SV.name for PRN in excluded_satellites)
                else "-"
            )

            # Handle longitude wrapping
            longitude_differences = np.abs(np.diff(sat_longitudes[:, i]))
            longitude_indices = np.where(longitude_differences > 300)

            if np.size(longitude_indices) > 0:
                # Split tracks at date line crossings
                for k, longitude_index in enumerate(longitude_indices[0]):
                    startIndex = 0 if k == 0 else longitude_indices[0][k - 1] + 1
                    end_index = longitude_index + 1

                    if len(list_satellites) < 50:
                        ax.plot(
                            sat_longitudes[startIndex:end_index, i],
                            sat_latitudes[startIndex:end_index, i],
                            linewidth=2,
                            color=satColor,
                            linestyle=line_style,
                            marker=".",
                            markersize=6,
                            transform=ccrs.PlateCarree(),
                            label=satellite_labels[-1],
                        )
                    else:
                        ax.plot(
                            sat_longitudes[startIndex:end_index, i],
                            sat_latitudes[startIndex:end_index, i],
                            linewidth=2,
                            color=satColor,
                            linestyle=line_style,
                            marker=".",
                            markersize=6,
                            transform=ccrs.PlateCarree(),
                        )
            else:
                # Plot continuous tracks
                if len(list_satellites) < 50:
                    ax.plot(
                        sat_longitudes[:, i],
                        sat_latitudes[:, i],
                        linewidth=2,
                        color=satColor,
                        linestyle=line_style,
                        marker=".",
                        markersize=6,
                        transform=ccrs.PlateCarree(),
                        label=satellite_labels[-1],
                    )
                else:
                    ax.plot(
                        sat_longitudes[:, i],
                        sat_latitudes[:, i],
                        linewidth=2,
                        color=satColor,
                        linestyle=line_style,
                        marker=".",
                        markersize=6,
                        transform=ccrs.PlateCarree(),
                    )

            # Add hour labels
            prn_hour_latitudes = sat_latitudes[:, i][indexHour]
            prn_longitude_hours = sat_longitudes[:, i][indexHour]

            for j, hr in enumerate(hourTxt[0]):
                if ~np.isnan(prn_hour_latitudes[j]):
                    hr = int(float(hr) / 3600.0)
                    plt.text(
                        prn_longitude_hours[j],
                        prn_hour_latitudes[j],
                        hr,
                        fontsize="x-small",
                        color=satColor,
                        transform=ccrs.PlateCarree(),
                    )

        # Set map extent to cover all data
        ax.set_global()

        # Add legend
        if len(list_satellites) < 50:
            my_legend = plt.legend(
                bbox_to_anchor=(0.5, 0.05),
                loc="lower center",
                ncol=min(np.size(satellite_labels), 8),
                fontsize="xx-small",
                markerscale=0.6,
            )
            for legobj in my_legend.legend_handles:
                legobj.set_linewidth(5.0)

        # Add title
        yyyy, mm, dd = observer.get_ymd()
        plt.title(
            f'{cli_args.gnss.replace(",", " & ").upper()} Satellite Groundtracks - Date {yyyy:04d}/{mm:02d}/{dd:02d}',
            fontsize="x-large",
        )

        # Save figure
        plot_ofn = f'/tmp/{observer.name}_{cli_args.gnss.replace(",", "_")}-{observer.get_ymd_str()}_groundtrack.png'
        plt.savefig(plot_ofn)
        if cli_args.verbose:
            print(f"\tSaved figure {plot_ofn}")

        if cli_args.verbose:
            plt.draw()


def plot_dop_visible(
    observer: Station,
    list_satellites: List[EarthSatellite],
    prediction_dt: List[datetime],
    elevations: np.ndarray,
    xDOPs: np.ndarray,
    cli_args: argparse.Namespace,
    excluded_satellites: List[str] = None,
):
    """Plots the Dilution of Precision (DOP) and the number of visible satellites over time.

    Args:
        observer (Station): The observing station.
        list_satellites (List[EarthSatellite]): List of Earth satellites.
        prediction_dt (List[datetime]): List of prediction datetime objects.
        elevations (np.ndarray): Array of satellite elevations.
        xDOPs (np.ndarray): Array of DOP values (HDOP, VDOP, TDOP).
        cli_args (argparse.Namespace): Command line arguments.
        excluded_satellites (List[str], optional): List of excluded satellites. Defaults to None.

    Returns:
        None
    """
    with rich_console.status(
        "Plotting xDOPS and number of visible satellites", spinner="aesthetic"
    ):

        plt.style.use("ggplot")

        max_DOP = 10

        fig = plt.figure(figsize=(14.0, 9.0))
        ax1 = plt.gca()
        ax2 = ax1.twinx()  # second y-axis needed, so make the x-axis twins
        # set colormap

        # plot the number of visible satellites
        nr_visible_satellites = []
        for i, el in enumerate(elevations):
            nr_visible_satellites.append(np.count_nonzero(~np.isnan(el)))
            ax2.set_ylim(0, max(nr_visible_satellites) + 1)

        ax2.plot(
            prediction_dt,
            nr_visible_satellites,
            linewidth=3,
            color="black",
            drawstyle="steps-post",
            label="#Visible",
            alpha=0.6,
        )

        # draw the line representing the number of SVs used
        if excluded_satellites is not None:
            nr_used_satellites = nr_visible_satellites - len(excluded_satellites)
        else:
            nr_used_satellites = nr_visible_satellites
        ax2.plot(
            prediction_dt,
            nr_used_satellites,
            linewidth=3,
            color="green",
            drawstyle="steps-post",
            label="#Visible",
            alpha=0.6,
        )

        # plot the xDOPS on first axis
        ax1.set_ylim(0, max_DOP)
        colors = iter(cm.jet(np.linspace(0, 1, len(xDOPs[0, :]) + 2)))
        labels = ["HDOP", "VDOP", "TDOP"]

        for i in range(0, 3):
            xDOP = xDOPs[:, i]
            dopColor = next(colors)
            transparency = 0.5 - i * 0.1
            ax1.fill_between(prediction_dt, 0, xDOP, color=dopColor, alpha=transparency)
            ax1.plot(prediction_dt, xDOP, linewidth=2, color=dopColor, label=labels[i])

            # add PDOP
            if i == 1:
                PDOP2 = xDOPs[:, 0] * xDOPs[:, 0] + xDOPs[:, 1] * xDOPs[:, 1]
                dopColor = next(colors)
                transparency = 0.2
                ax1.fill_between(
                    prediction_dt, 0, np.sqrt(PDOP2), color=dopColor, alpha=transparency
                )
                ax1.plot(
                    prediction_dt,
                    np.sqrt(PDOP2),
                    linewidth=2,
                    color=dopColor,
                    label="PDOP",
                )

            # add GDOP
            if i == 2:
                GDOP = np.sqrt(PDOP2 + xDOPs[:, 2] * xDOPs[:, 2])
                dopColor = next(colors)
                transparency = 0.1
                ax1.fill_between(
                    prediction_dt, 0, GDOP, color=dopColor, alpha=transparency
                )
                ax1.plot(prediction_dt, GDOP, linewidth=2, color=dopColor, label="GDOP")

        ax1.legend(loc="upper left", frameon=True)
        ax1.set_ylim(0, max_DOP)

        plt.title(
            f"{cli_args.gnss.replace(',', ' & ').upper()} Satellite Visibility",
            fontsize="x-large",
        )
        annotateTxt = (
            f"Station: {observer.name} @ (ϕ {ephem.degrees(observer.lat)}, λ {ephem.degrees(observer.lon)}) - "
            f"Date {observer.get_ymd_str()} - cutoff {cli_args.cutoff:2d}"
        )
        plt.text(
            0.5,
            0.99,
            annotateTxt,
            horizontalalignment="center",
            verticalalignment="top",
            transform=ax1.transAxes,
            fontsize="medium",
        )
        ax1.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        ax1.xaxis.set_major_locator(
            mdates.MinuteLocator(byminute=0, interval=1, tz=None)
        )

        plt.xticks(rotation=50, size="medium")
        ax1.set_ylabel("DOP [-]", fontsize="x-large")
        ax2.set_ylabel("# Visible satellites [-]", fontsize="x-large")
        ax1.set_xlabel("Time of Day", fontsize="x-large")

        plot_ofn = f"/tmp/{observer.name}_{cli_args.gnss.replace(',', '_')}_{observer.get_ymd_str()}_DOP.png"
        fig.savefig(plot_ofn, dpi=fig.dpi)
        if cli_args.verbose:
            print(f"\tSaved figure {plot_ofn}")

        if cli_args.verbose:
            plt.draw()
