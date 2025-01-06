"""
UFS-Coastal Meteorological Forcing Downloader
=============================================

This script downloads meteorological data (HRRR or GFS) to be used within the UFS-Coastal infrastructure.
It supports both GRIB2 and NetCDF output formats and allows specification of
spatial domains either through a hgrid.gr3 file or explicit bounding box coordinates.

Features:
---------
- Supports both HRRR and GFS weather models (ERA5 will be included in later versions)
- Downloads data in either GRIB2 or NetCDF format
- Automatic domain extraction from hgrid.gr3 (currently effective for NetCDF format)
- Custom bounding box specification
- Flexible variable selection
- Temporary file management


Usage Examples:
--------------
1. Basic HRRR download for CONUS domain:
   ```
   python download_met.py --start-date "2024-01-05 00:00" --rnday 2 --model hrrr
   ```

2. GFS download with custom bounding box:
   ```
   python download_met.py --start-date "2024-01-05 00:00" --rnday 5 \
                         --model gfs \
                         --bbox -80 25 -70 35
   ```

3. Download specific variables in GRIB2 format:
   ```
   python download_met.py --start-date "2024-01-05 00:00" --rnday 1 \
                         --format grib2 \
                         --var uwind vwind prmsl
   ```

4. Using hgrid.gr3 for domain specification:
   ```
   python download_met.py --start-date "2024-01-05 00:00" --rnday 3 \
                         --hgrid /path/to/hgrid.gr3 \
                         --outdir ./sflux_dir
   ```

Arguments:
---------
Required:
  --start-date    Start date in format YYYY-mm-dd HH:MM
  --rnday         Number of days to download

Optional:
  --model         Model selection [hrrr|gfs] (default: hrrr)
  --format        Output format [netcdf|grib2] (default: netcdf)
  --hgrid         Path to hgrid.gr3 file for domain
  --bbox          Bounding box [lon_min lat_min lon_max lat_max]
  --region        Region for HRRR [conus|alaska] (default: conus)
  --outdir        Output directory (default: ./sflux)
  --pscr          Scratch directory for temporary files
  --var           Variables to download (default: all)
                  Options: uwind, vwind, prmsl, spfh, stmp, prate, dlwrf, dswrf
  --level         Output level (default: 2)
  --debug         Enable debug logging

Available Variables:
------------------
- uwind:  U-component of wind
- vwind:  V-component of wind
- prmsl:  Pressure at mean sea level
- spfh:   Specific humidity
- stmp:   Surface temperature
- prate:  Precipitation rate
- dlwrf:  Downward longwave radiation flux
- dswrf:  Downward shortwave radiation flux

Output Structure:
---------------
For NetCDF format:
./sflux/
├── air/
│   └── sflux_air_YYYYMM.nc
├── rad/
│   └── sflux_rad_YYYYMM.nc
└── prc/
    └── sflux_prc_YYYYMM.nc

For GRIB2 format:
./sflux/
└── {model}/
    └── {model}_YYYYMMDD_HH.grib2

Notes:
-----
1. The script automatically handles temporary file cleanup unless --pscr is specified
2. For HRRR, data is available hourly for up to 48 hours
3. For GFS, data is available at 3-hour intervals for up to 384 hours
4. Default bounding box for CONUS: [-130, 20, -60, 55]
5. When using hgrid.gr3, a 0.1-degree buffer is added to the domain

Requirements:
------------
- Python 3.6+
- numpy
- pyschism
- logging

Contact:

Mansur Ali Jisan (mansur.jisan@noaa.gov)
NOAA/NOS/CO-OPS

Version: 1.0
Last Updated: January 2025

"""

import argparse
from datetime import datetime, timedelta
import logging
import os
import shutil
import sys
import tempfile
import numpy as np
from types import SimpleNamespace
from pyschism.forcing.nws.nws2.gfs2 import GFS, AWSGrib2Inventory as GFSGrib2Inventory
from pyschism.forcing.nws.nws2.hrrr3 import HRRR, AWSGrib2Inventory as HRRRGrib2Inventory

logger = logging.getLogger(__name__)

def setup_logging(level=logging.INFO):
    """Setup basic logging configuration"""
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def create_bbox_object(bbox_coordinates):
    """Create a bbox object with proper attributes from coordinates"""
    if bbox_coordinates is None:
        # Default bbox for CONUS if none provided
        return SimpleNamespace(
            xmin=-130.0,
            ymin=20.0,
            xmax=-60.0,
            ymax=55.0
        )
        
    return SimpleNamespace(
        xmin=float(bbox_coordinates[0]),
        ymin=float(bbox_coordinates[1]),
        xmax=float(bbox_coordinates[2]),
        ymax=float(bbox_coordinates[3])
    )

def read_hgrid(hgrid_file):
    """Read hgrid.gr3 file and extract bounding box"""
    try:
        with open(hgrid_file, 'r') as f:
            next(f)
            ne, np = map(int, f.readline().strip().split())
            
            lons = []
            lats = []
            for i in range(np):
                line = f.readline().strip().split()
                lon = float(line[1])
                lat = float(line[2])
                lons.append(lon)
                lats.append(lat)
                
        buffer = 0.1
        bbox_coords = [
            min(lons) - buffer,
            min(lats) - buffer,
            max(lons) + buffer,
            max(lats) + buffer
        ]
        
        logger.info(f'Bounding box from hgrid: {bbox_coords}')
        return create_bbox_object(bbox_coords)
        
    except Exception as e:
        logger.error(f'Error reading hgrid file: {str(e)}')
        sys.exit(1)


def parse_arguments():
    parser = argparse.ArgumentParser(description='Download HRRR/GFS data')
    
    # Required arguments
    parser.add_argument('--start-date', required=True,
                      help='Start date in format YYYY-mm-dd HH:MM')
    parser.add_argument('--rnday', type=float, required=True,
                      help='Number of days to download')
    
    # Model and format selection
    parser.add_argument('--model', choices=['hrrr', 'gfs'], default='hrrr',
                      help='Model to download from (default: hrrr)')
    parser.add_argument('--format', choices=['netcdf', 'grib2'], default='netcdf',
                      help='Output format (default: netcdf)')
    
    # Domain arguments
    parser.add_argument('--hgrid', type=str,
                      help='Path to hgrid.gr3 file for domain bounding box')
    parser.add_argument('--bbox', type=float, nargs=4,
                      metavar=('LON_MIN', 'LAT_MIN', 'LON_MAX', 'LAT_MAX'),
                      help='Bounding box coordinates')
    parser.add_argument('--region', choices=['conus', 'alaska'],
                      default='conus',
                      help='Region for HRRR (default: conus)')
    
    # Directory arguments
    parser.add_argument('--outdir', default='./sflux',
                      help='Output directory (default: ./sflux)')
    parser.add_argument('--pscr', 
                      help='Scratch directory for temporary files (default: system temp directory)')
    
    # Variable selection
    parser.add_argument('--var', nargs='+', 
                      choices=['uwind', 'vwind', 'prmsl', 'spfh', 'stmp', 'prate', 'dlwrf', 'dswrf'],
                      help='Variables to download (e.g., -var uwind vwind prmsl)')
    
    parser.add_argument('--level', type=int, default=2,
                      help='Output level (default: 2)')
    parser.add_argument('--debug', action='store_true',
                      help='Enable debug logging')
    
    args = parser.parse_args()
    
    # Set air, rad, prc flags based on variable selections
    if args.var:
        args.air = any(var in ['uwind', 'vwind', 'prmsl', 'spfh', 'stmp'] for var in args.var)
        args.rad = any(var in ['dlwrf', 'dswrf'] for var in args.var)
        args.prc = 'prate' in args.var
    else:
        # If no variables specified, download all
        args.var = ['uwind', 'vwind', 'prmsl', 'spfh', 'stmp', 'prate', 'dlwrf', 'dswrf']
        args.air = True
        args.rad = True
        args.prc = True
    
    return args


def main():
    args = parse_arguments()
    setup_logging(level=logging.DEBUG if args.debug else logging.INFO)

    try:
        start_date = datetime.strptime(args.start_date, '%Y-%m-%d %H:%M')
    except ValueError:
        logger.error('Invalid date format. Use YYYY-mm-dd HH:MM')
        sys.exit(1)

    bbox = None
    if args.hgrid:
        if not os.path.exists(args.hgrid):
            logger.error(f'hgrid file not found: {args.hgrid}')
            sys.exit(1)
        bbox = read_hgrid(args.hgrid)
    elif args.bbox:
        bbox = create_bbox_object(args.bbox)
    else:
        bbox = create_bbox_object(None)

    logger.info(f'Using bounding box: xmin={bbox.xmin}, ymin={bbox.ymin}, xmax={bbox.xmax}, ymax={bbox.ymax}')

    os.makedirs(args.outdir, exist_ok=True)

    if args.pscr:
        scratch_dir = args.pscr
    else:
        scratch_dir = tempfile.mkdtemp(prefix=f'{args.model}_')
        logger.info(f'Created temporary scratch directory: {scratch_dir}')

    logger.info(f"Model: {args.model.upper()}")

    try:
        if args.format == 'grib2':
            logger.info("Downloading in GRIB2 format")
            temp_dir = tempfile.mkdtemp(prefix=f'{args.model}_tmp_')
            
            try:
                # Create directory for model files
                model_dir = os.path.join(args.outdir, args.model)
                os.makedirs(model_dir, exist_ok=True)

                if args.model == 'hrrr':
                    inventory = HRRRGrib2Inventory(start_date, record=args.rnday, pscr=temp_dir)
                else:  # GFS
                    inventory = GFSGrib2Inventory(start_date, record=args.rnday, pscr=temp_dir)

                files = inventory.files
                if files:
                    import shutil
                    for file in files:
                        dest_file = os.path.join(model_dir, os.path.basename(file))
                        shutil.copy2(file, dest_file)
                        logger.info(f'Saved to: {dest_file}')
                    logger.info(f'Total files downloaded: {len(files)}')
                else:
                    logger.error("No files were downloaded")
            finally:
                # Clean up temporary directory used for downloads
                if os.path.exists(temp_dir):
                    import shutil
                    shutil.rmtree(temp_dir)

        else:  # NetCDF format
            logger.info("Downloading in NetCDF format")
            logger.info("Variables to download: " + ", ".join(args.var))

            if args.model == 'hrrr':
                model = HRRR(
                    level=args.level,
                    region=args.region,
                    bbox=bbox,
                    pscr=scratch_dir,
                    outdir=args.outdir
                )
            else:
                model = GFS(
                    level=args.level,
                    bbox=bbox,
                    pscr=scratch_dir,
                    outdir=args.outdir
                )

            model.write(
                start_date=start_date,
                rnday=args.rnday,
                air=args.air,
                prc=args.prc,
                rad=args.rad
            )

        logger.info('Download completed successfully')

    except Exception as e:
        logger.error(f'Error during download: {str(e)}')
        if args.debug:
            logger.exception(e)
        sys.exit(1)
    finally:
        if not args.pscr and os.path.exists(scratch_dir):
            try:
                import shutil
                shutil.rmtree(scratch_dir)
                logger.debug(f'Cleaned up temporary directory: {scratch_dir}')
            except Exception as e:
                logger.warning(f'Failed to clean up temporary directory: {e}')

if __name__ == '__main__':
    main()
