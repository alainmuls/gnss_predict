from typing import Optional
import argparse
import os

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments for GNSS orbit prediction."""
    parser = argparse.ArgumentParser(
        description=f'{os.path.basename(__file__)} predicts GNSS orbits based on TLEs'
    )
    
    parser.add_argument(
        '-g', '--gnss',
        help='Name of GNSSs as comma separated list (cfr NORAD naming)',
        required=True
    )
    
    parser.add_argument(
        '-x', '--exclude',
        help='Comma separated list of satellite PRNs to exclude from DOP calculation (eg. E18,E14,E20)',
        default=''
    )
    
    parser.add_argument(
        '-i', '--interval',
        help='Interval in minutes',
        type=int,
        default=20
    )
    
    parser.add_argument(
        '-c', '--cutoff',
        help='Cutoff angle in degrees',
        type=int,
        default=10
    )
    
    parser.add_argument(
        '-o', '--observer',
        help='Station info "name,latitude,longitude" (units = degrees, defaults to RMA)',
        default=None
    )
    
    parser.add_argument(
        '-d', '--date',
        help='Prediction date (YYYY/MM/DD), defaults to today',
        default=None
    )
    
    parser.add_argument(
        '-s', '--start',
        help='Start time (hh:mm)',
        default='00:00'
    )
    
    parser.add_argument(
        '-e', '--end',
        help='End time (hh:mm)',
        default='23:59'
    )
    
    parser.add_argument(
        '-m', '--max-dop',
        help='Maximum xDOP value to display',
        type=int,
        default=10
    )
    
    parser.add_argument(
        '-v', '--verbose',
        help='Display interactive graphs and increase output verbosity',
        action='store_true'
    )

    return parser.parse_args()
