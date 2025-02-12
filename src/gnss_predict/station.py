import re
import sys
from datetime import datetime, timedelta
from typing import List, Optional, Tuple

import ephem
import numpy as np
from rich import print


class Station(ephem.Observer):
    """
    Station class for managing observer coordinates and data.

    Inherits from ephem.Observer to provide astronomical calculations.
    Handles station initialization, parsing and data representation.
    """

    # Create default stations
    RMA = None
    BERTRIX = None
    _initialized = False

    @classmethod
    def initialize_defaults(cls):
        if not cls._initialized:
            # Initialize RMA
            cls.RMA = Station()
            cls.RMA.init(
                "RMA", "50:50:38.4551", "4:23:34.5421", ephem.date(ephem.now())
            )

            # Initialize BERTRIX
            cls.BERTRIX = Station()
            cls.BERTRIX.init(
                "BERTRIX", "49.894275", "5.241417", ephem.date(ephem.now())
            )

            cls._initialized = True

    def __init__(self) -> None:
        """Initialize a new Station instance with empty name."""
        self.name: str = ""
        super(Station, self).__init__()

    def init(self, name: str, lat: str, lon: str, date: ephem.Date) -> None:
        """
        Initialize station with location and date information.

        Args:
            name: Station identifier
            lat: Latitude in degrees
            lon: Longitude in degrees
            date: Observation date
        """
        self.name = name
        self.lat = ephem.degrees(lat)
        self.lon = ephem.degrees(lon)
        self.date = ephem.date(date)

    def parse(self, text: str) -> None:
        """
        Parse station information from comma-separated string.

        Args:
            text: Comma-separated string containing name, latitude, longitude
            Format: "name,latitude,longitude"
        """
        stations_elements: list[str] = re.split(",", text)
        Location: Station = Station()
        Location.init(
            stations_elements[0],
            stations_elements[1],
            stations_elements[2],
            ephem.date(ephem.now()),
        )

        if np.size(stations_elements) == 3:
            self.name = stations_elements[0]
            self.lat = ephem.degrees(stations_elements[1])
            self.lon = ephem.degrees(stations_elements[2])
        else:
            sys.stderr.write("wrong number of elements to parse\n")

    def get_ymd(self) -> Tuple[int, int, int]:
        """
        Extract year, month, day from date structure.

        Returns:
            Tuple containing (year, month, day)
        """
        dateTxt: Tuple = ephem.date(self.date).triple()
        year: int = int(dateTxt[0])
        month: int = int(dateTxt[1])
        day: int = int(dateTxt[2])

        return year, month, day

    def get_ymd_str(self) -> str:
        """
        Get date as string in "YYYYMMDD" format.

        Returns:
            Date string in "YYYYMMDD" format
        """
        yr, mm, dd = self.get_ymd()
        return f"{yr:04d}{mm:02d}{dd:02d}"

    def __str__(self) -> str:
        """Return formatted string representation of station."""
        formatted_date = self.date.datetime().strftime('%Y/%m/%d')
        lat_dms = str(ephem.degrees(self.lat)).split(':')
        lon_dms = str(ephem.degrees(self.lon)).split(':')
        
        return (
            f"Station: {self.name}\n"
            f"\t Latitude: {lat_dms[0]}° {lat_dms[1]}' {float(lat_dms[2]):.1f}\"\n"
            f"\t Longitude: {lon_dms[0]}° {lon_dms[1]}' {float(lon_dms[2]):.1f}\"\n"
            f"\t Date: {formatted_date}"
        )

    def __repr__(self) -> str:
        """Return formatted string representation of station."""
        formatted_date = self.date.datetime().strftime('%Y/%m/%d')
        lat_dms = str(ephem.degrees(self.lat)).split(':')
        lon_dms = str(ephem.degrees(self.lon)).split(':')
        
        return (
            f"Station: {self.name}\n"
            f"\t Latitude: {lat_dms[0]}° {lat_dms[1]}' {float(lat_dms[2]):.1f}\"\n"
            f"\t Longitude: {lon_dms[0]}° {lon_dms[1]}' {float(lon_dms[2]):.1f}\"\n"
            f"\t Date: {formatted_date}"
        )

    def create_observer(
        self,
        station: Optional[str] = None,
        prediction_date: Optional[str] = None,
        verbose: bool = False,
    ) -> "Station":
        """
        Create and configure a Station observer instance.

        Args:
            station: Station info string "name,lat,lon" or None for RMA default
            prediction_date: Date for prediction or None for current date
            verbose: Enable detailed output printing

        Returns:
            Configured Station observer instance
        """
        observer: Station = cls()

        if station is None:
            observer = self.RMA
        else:
            observer.parse(station)

        if prediction_date is None:
            observer.date = ephem.date(ephem.now())
        else:
            observer.date = ephem.Date(prediction_date)

        return observer

    @classmethod
    def set_observer_data(
        cls,
        station: Optional[str] = None,
        prediction_date: Optional[str] = None,
        verbose: bool = False,
    ) -> "Station":
        """
        Sets the station information for calculations.

        Args:
            station: Station info string "name,lat,lon" or None for RMA default
            prediction_date: Date for prediction or None for current date
            verbose: Enable detailed output

        Returns:
            Configured Station observer
        """
        if station is None:
            observer = cls.RMA
        else:
            observer.parse(station)

        if prediction_date is None:
            observer.date = ephem.date(ephem.now())
        else:
            observer.date = ephem.Date(prediction_date)

        if verbose:
            formatted_date = observer.date.datetime().strftime("%Y/%m/%d")
            print(f"Prediction for {formatted_date}")

        return observer

    def set_observation_times(
        self, time_start: str, time_end: str, interval_min: int, verbose: bool = False
    ) -> Tuple[List[datetime], int]:
        """
        Calculates prediction times based on start/end times and interval.

        Args:
            time_start: Start time in "HH:MM" format
            time_end: End time in "HH:MM" format
            interval_min: Time interval in minutes
            verbose: Enable detailed output

        Returns:
            Tuple containing:
                - List of prediction times
                - Number of predictions

        Raises:
            SystemExit: If end time is before start time
        """
        yyyy, mm, dd = self.get_ymd()
        start_hour, start_min = map(int, time_start.split(":"))
        end_hour, end_min = map(int, time_end.split(":"))

        start_datetime = datetime(
            yyyy,
            mm,
            dd,
            hour=start_hour,
            minute=start_min,
            second=0,
            microsecond=0,
            tzinfo=None,
        )

        end_datetime = datetime(
            yyyy,
            mm,
            dd,
            hour=end_hour,
            minute=end_min,
            second=0,
            microsecond=0,
            tzinfo=None,
        )

        if end_datetime <= start_datetime:
            sys.stderr.write(
                f"end time [red]{end_datetime}[/red] is less than start time [red]{start_datetime}[/red].\nProgram exits.\n"
            )
            sys.exit(1)

        dt_minutes = (end_datetime - start_datetime).total_seconds() / 60
        nr_predictions = int(dt_minutes / float(interval_min)) + 1

        obs_dates = [
            start_datetime + timedelta(minutes=(interval_min * x))
            for x in range(nr_predictions)
        ]

        if verbose:
            print(
                f"Observation time span from {obs_dates[0]} to {obs_dates[-1]} "
                f"with interval {interval_min} min (#{len(obs_dates)})"
            )

        return obs_dates, nr_predictions
