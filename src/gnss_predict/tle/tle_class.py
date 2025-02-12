import fileinput
import sys
from pathlib import Path
from typing import List

import ephem
import requests
from rich import print


class TLEManager:
    """Downloads and manages TLE files from NORAD/CelesTrak."""

    BASE_URL = "https://www.celestrak.com/NORAD/elements"
    TMP_DIR = Path("/tmp")

    @staticmethod
    def get_tle_from_norad(gnss_systems: str, verbose: bool = False) -> str:
        """
        Downloads latest TLE for satellite system or reuses existing local file.

        # Flow:
        # 1. Split input into satellite systems if multiple provided
        # 2. Download TLE for each system
        # 3. Combine files if multiple systems requested
        # 4. Return path to final TLE file

        Args:
            gnss: comma separated basename of TLE files (from NORAD site)
            verbose: Enable verbose output

        Returns:
            Path to downloaded/reused TLE file

        Raises:
            SystemExit: On connection/download failure
        """
        if verbose:
            print(f"Downloading TLEs from NORAD for satellite systems: {gnss_systems}")

        # Split constellation list if multiple provided
        sat_systems = gnss_systems.split(",")
        tle_filenames: List[Path] = []

        # Download TLE for each satellite system
        for gnss in sat_systems:
            url = f"{TLEManager.BASE_URL}/{gnss}.txt"
            output_file = TLEManager.TMP_DIR / f"{gnss}.txt"
            print(f"\nDownloading TLEs from {url} to {output_file}")

            try:
                TLEManager._download_tle(url, output_file)
                tle_filenames.append(output_file)
            except requests.exceptions.RequestException:
                print(f"Connection to NORAD could not be established for TLE {gnss}.")

                if output_file.exists():
                    print(f"Using local file {output_file}")
                    tle_filenames.append(output_file)

        # Combine multiple TLE files if needed
        if len(tle_filenames) > 0:
            if len(tle_filenames) > 1:
                combined_name = gnss_systems.replace(",", "-") + ".txt"
                output_file = TLEManager.TMP_DIR / combined_name

                with output_file.open("w") as fout:
                    with fileinput.input(files=tle_filenames) as fin:
                        for line in fin:
                            fout.write(line)
            else:
                output_file = tle_filenames[0]

            if verbose:
                print(f"(Combined) TLE file saved in {output_file}")
        else:
            print("No TLE files found. Exiting program.")
            sys.exit(10)

        return str(output_file)

    @staticmethod
    def _download_tle(url: str, output_file: Path) -> None:
        """
        Downloads TLE file from given URL.

        Uses streaming to handle large files efficiently.
        Writes chunks directly to disk to minimize memory usage.

        Args:
            url: URL to download TLE file from
            output_file: Path to save downloaded file
        """
        response = requests.get(url, stream=True)
        response.raise_for_status()

        with output_file.open("wb") as f:
            for chunk in response.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
                    f.flush()

    @staticmethod
    def load_tle(
        tle_filename: Path | str, verbose: bool = False
    ) -> List[ephem.EarthSatellite]:
        """
        Loads a TLE file and creates a list of satellites.

        Args:
            tle_filename: Path to TLE file
            verbose: Enable detailed output

        Returns:
            List of decoded satellite objects
        """
        satellites = []

        with open(tle_filename) as f:
            line1 = f.readline()
            while line1:
                line2 = f.readline()
                line3 = f.readline()
                satellite = ephem.readtle(line1, line2, line3)
                satellites.append(satellite)
                line1 = f.readline()

        if verbose:
            print(f"{len(satellites)} satellites loaded\n")

        return satellites
