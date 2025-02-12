
# GNSS Tracking

A Python package for tracking Global Navigation Satellite Systems (GNSS) satellites and analyzing their positions.

## Features

- Real-time GNSS satellite tracking
- Satellite position visualization using Cartopy
- TLE (Two-Line Element) data processing
- Rich command-line interface
- Comprehensive satellite data analysis tools

## Installation

Install the package using pip:

```bash
pip install gnss_predict
```

## Requirements
The package requires the following Python libraries:

- ephem
- numpy
- matplotlib
- cartopy
- rich
- pandas


## Usage

After installation, you can run the GNSS tracking tool from the command line:
```python
± gnss_predict -h
usage: gnss_predict [-h] -g GNSS [-x EXCLUDE] [-i INTERVAL] [-c CUTOFF] [-o OBSERVER] [-d DATE] [-s START] [-e END] [-m MAX_DOP] [-v]

cli.py predicts GNSS orbits based on TLEs

options:
  -h, --help            show this help message and exit
  -g GNSS, --gnss GNSS  Name of GNSSs as comma separated list (cfr NORAD naming)
  -x EXCLUDE, --exclude EXCLUDE
                        Comma separated list of satellite PRNs to exclude from DOP calculation (eg. E18,E14,E20)
  -i INTERVAL, --interval INTERVAL
                        Interval in minutes
  -c CUTOFF, --cutoff CUTOFF
                        Cutoff angle in degrees
  -o OBSERVER, --observer OBSERVER
                        Station info "name,latitude,longitude" (units = degrees, defaults to RMA)
  -d DATE, --date DATE  Prediction date (YYYY/MM/DD), defaults to today
  -s START, --start START
                        Start time (hh:mm)
  -e END, --end END     End time (hh:mm)
  -m MAX_DOP, --max-dop MAX_DOP
                        Maximum xDOP value to display
  -v, --verbose         Display interactive graphs and increase output verbosity
```

## Data Sources
The package includes:

- __GNSS__ satellite time functions
- __TLE__ (Two-Line Element) files
- __Utility__ functions for data processing

## Output

- Data files
  - [Visible satellites file](./md/RMA_beidou_20250211_VIS.txt)
  - [xDOP file](./md/RMA_beidou_20250211_DOP.txt)
  - [Coordinates for groundtrack of satellites file](./md/RMA_beidou_20250211_GEOD.txt)
  - [Satellites overview arc file](./md/RMA_beidou_20250211_TLE_arc.csv)
- Plots
  - [Visible satellites plot](./md/RMA_beidou_20250211_visibility.png)
  - [xDOP plot](./md/RMA_beidou_20250211_DOP.png)
  - [Groundtrack plot](./md/RMA_beidou-20250211_groundtrack.png)
  - [Satellites overview arc plot](./md/RMA_beidou_20250211_skyview.png)



## Contributing
1. Fork the repository
1. Create your feature branch (git checkout -b feature/amazing-feature)
1. Commit your changes (git commit -m 'Add amazing feature')
1. Push to the branch (git push origin feature/amazing-feature)
1. Open a Pull Request


## License

BSD License.

## Contact

[Alain Muls](alain.muls@gmail.com)
