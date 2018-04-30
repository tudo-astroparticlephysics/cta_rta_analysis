import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.coordinates import EarthLocation, Latitude, Longitude
from astropy.coordinates.angle_utilities import angular_separation

from astropy import wcs

from dateutil import parser
import numpy as np


def select_points_in_ring(coords, center, inner_radius=0.32 * u.deg, outter_radius=1 * u.deg):
    d = coords.separation(center)
    m = (d > inner_radius) & (d < outter_radius)
    return m


def select_points_in_circle(coords, center, radius=0.32 * u.deg):
    d = coords.separation(center)
    m = d < radius
    return m


def split_on_off(gammas, protons, on_region_radius=0.32 * u.deg, off_region_radius=2.0 * u.deg, center=None):
    c_gamma = skyccords_from_dl3_table(gammas).icrs
    c_proton = skyccords_from_dl3_table(protons).icrs

    if center is None:
        center = horizontal_to_skycoord(70 * u.deg, 0 * u.deg)

    on_gammas = select_points_in_circle(c_gamma, center, radius=on_region_radius)
    on_protons = select_points_in_circle(c_proton, center, radius=on_region_radius)

    inner_ring_radius = 0.4 * u.deg
    off_gammas = select_points_in_ring(c_gamma, center, inner_radius=inner_ring_radius, outter_radius=off_region_radius)
    off_protons = select_points_in_ring(c_proton, center, inner_radius=inner_ring_radius, outter_radius=off_region_radius)

    off_events = gammas.iloc[off_gammas].append(protons.iloc[off_protons])
    on_events = gammas.iloc[on_gammas].append(protons.iloc[on_protons])

    off_area = np.pi * (off_region_radius**2 - inner_ring_radius**2)
    on_area = np.pi * on_region_radius**2
    alpha = (on_area / off_area).si.value
    return on_events, off_events, alpha


def build_standard_wcs(image_center, shape, naxis=2, fov=9 * u.deg):
    width, height = shape

    w = wcs.WCS(naxis=2)
    w.wcs.crpix = [width / 2 + 0.5, height / 2 + 0.5]
    w.wcs.cdelt = np.array([-fov.value / width, fov.value / height])
    w.wcs.crval = [image_center.ra.deg, image_center.dec.deg]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    w.wcs.radesys = 'ICRS'
    w.wcs.equinox = 2000.0
    w.wcs.cunit = ['deg', 'deg']
    w._naxis = [width, height]
    return w


def calculate_distance_theta(df, source_alt=70 * u.deg, source_az=0 * u.deg):
    source_az = Angle(source_az).wrap_at(180 * u.deg)
    source_alt = Angle(source_alt)

    az = Angle(df.az_prediction.values, unit=u.rad).wrap_at(180*u.deg)
    alt = Angle(df.alt_prediction.values, unit=u.rad)

    distance = angular_separation(source_az, source_alt, az, alt)
    return distance



def wrap_angles(alt, az):
    alt = alt.to('degree')
    az = Angle(az).wrap_at(180 * u.deg).degree * u.deg
    return alt, az


def skyccords_from_dl3_table(df):
    az = df.az_prediction.values * u.rad
    alt = df.alt_prediction.values * u.rad
    return horizontal_to_skycoord(alt, az)


@u.quantity_input(alt=u.rad, az=u.rad)
def horizontal_to_skycoord(alt, az, location=None, dt=None):

    alt, az = wrap_angles(alt, az)

    lat = Latitude((24, 37, 38), unit='deg')
    lon = Longitude((70, 34, 15), unit='deg')
    paranal = EarthLocation.from_geodetic(lon, lat, 2600)
    # paranal = EarthLocation.of_site('paranal')
    dt = parser.parse('2017-09-20 20:15')

    c = SkyCoord(
        alt=alt,
        az=az,
        obstime=dt,
        frame='altaz',
        location=paranal,
    )

    return c


@u.quantity_input(alt1=u.rad, az1=u.rad, alt2=u.rad, az2=u.rad,)
def distance_between_horizontal_coordinates(alt1, az1, alt2, az2):

    alt, az = wrap_angles(alt1, az1)
    mc_alt, mc_az = wrap_angles(alt2, az2)

    lat = Latitude((24, 37, 38), unit='deg')
    lon = Longitude((70, 34, 15), unit='deg')
    paranal = EarthLocation.from_geodetic(lon, lat, 2600)
    # paranal = EarthLocation.of_site('paranal')
    dt = parser.parse('2017-09-20 22:15')

    c = SkyCoord(
        alt=alt,
        az=az,
        obstime=dt,
        frame='altaz',
        location=paranal,
    )

    c_mc = SkyCoord(
        alt=mc_alt,
        az=mc_az,
        obstime=dt,
        frame='altaz',
        location=paranal,
    )

    return c.separation(c_mc)
