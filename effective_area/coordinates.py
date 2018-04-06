import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.coordinates import EarthLocation, Latitude, Longitude
from dateutil import parser


def wrap_angles(alt, az):
    alt = alt.to('degree')
    az = Angle(az).wrap_at(180 * u.deg).degree * u.deg

    return alt, az


@u.quantity_input(alt1=u.rad, az1=u.rad)
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
