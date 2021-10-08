#!/usr/bin/env python
""" Earth-Mars """
"""
	earth-mars.py - a program to create the Earth Mars distance calculation plus occultation
	Martin J Levy - W6LHI/G8LHI - https://github.com/mahtin/mars-earth
	Copyright (C) 2021 @mahtin - https://github.com/mahtin/earth-mars/blob/main/LICENSE
"""

import math
import time
import datetime

import ephem

class EarthMars(object):
	""" caculate Earth/Mars relationship based on Earth and Mars lat/long """

	def __init__(self, debug=False):
		""" Everything is based on a constellation name """

		self._debug = debug
		self._data = None
		self._earthobserver = None
		self._marsobserver = None

	def __call__(self):
		""" return whatever we have! """
		return self._data

	def earthobserver(self, lon=0.0, lat=0.0, amsl=0.0):
		""" set Earth observer location """

		if lon > 180 or lon < -180 or lat > 90 or lat < -90:
			raise ValueError

		if not self._earthobserver:
			self._earthobserver = ephem.Observer()
		self._earthobserver.lon = self._degrees_to_radians(lon)
		self._earthobserver.lat = self._degrees_to_radians(lat)
		self._earthobserver.elevation = amsl
		return True

	def sun(self, when=None):

		if not self._earthobserver:
			raise RuntimeError

		if when:
			self._earthobserver.date = when
		else:
			# All in UTC - https://rhodesmill.org/pyephem/date.html-
			self._earthobserver.date = datetime.datetime.utcnow()

		sun = ephem.Sun()
		sun.compute(self._earthobserver)

		r = {
			'name': sun.name,
			'mag': sun.mag,
			'radius': self._radians_to_degrees(sun.radius),
			'ra_dec': [self._radians_to_degrees(sun.ra), self._radians_to_degrees(sun.dec)],
			'alt_az': [self._radians_to_degrees(sun.alt), self._radians_to_degrees(sun.az)],
			'earth_distance_km': self._au_to_km(sun.earth_distance),
		}
		return r

	def mars(self, when=None):

		if not self._earthobserver:
			raise RuntimeError

		if when:
			self._earthobserver.date = when
		else:
			# All in UTC - https://rhodesmill.org/pyephem/date.html-
			self._earthobserver.date = datetime.datetime.utcnow()

		mars = ephem.Mars()
		mars.compute(self._earthobserver)

		previous_rising = self._earthobserver.previous_rising(mars).datetime()
		previous_setting = self._earthobserver.previous_setting(mars).datetime()
		next_rising = self._earthobserver.next_rising(mars).datetime()
		next_setting = self._earthobserver.next_setting(mars).datetime()

		r = {
			'name': mars.name,
			'mag': mars.mag,
			'radius': self._radians_to_degrees(mars.radius),
			'ra_dec': [self._radians_to_degrees(mars.ra), self._radians_to_degrees(mars.dec)],
			'alt_az': [self._radians_to_degrees(mars.alt), self._radians_to_degrees(mars.az)],
			'earth_distance_km': self._au_to_km(mars.earth_distance),
			'previous_rising': previous_rising,
			'previous_setting': previous_setting,
			'next_rising': next_rising,
			'next_setting': next_setting,
		}
		return r

	def now(self, when=None, horizon_angle=0.0, historical=False):
		""" return where the satellites are """
		if self._constellation and not self._data:
			self._update()
		now_utc = datetime.datetime.utcnow()
		too_old = datetime.date.today() - datetime.timedelta(21)
		r = {}
		r['satellite'] = {}
		r['directory'] = {}
		r['now'] = now_utc.replace(microsecond=0).isoformat(timespec='seconds')
		if not when:
			when = now_utc
		else:
			# is the provided 'when' value valid? - worry about later
			pass

		if self._observer:
			# All in UTC - https://rhodesmill.org/pyephem/date.html-
			self._observer.date = when
			arg = self._observer
		else:
			arg = when

		for satelliteId in self._data:
			s = self._data[satelliteId]
			if self._satellite:
				if isinstance(self._satellite, str) and self._satellite != s['name'].lower():
					continue;
				elif isinstance(self._satellite, int) and self._satellite != s['id']:
					continue;
			if 'valid' not in s:
				s['valid'] = True	# first time thru - no idea if satellite has valid tle
			if not s['valid']:
				continue
			if 'tle_rec' not in s or s['tle_rec'] == None:
				s['tle_rec'] = ephem.readtle(s['name'], s['line1'], s['line2'])

			# See documentation at https://www.celestrak.com/NORAD/documentation/tle-fmt.php

			tle_rec = s['tle_rec']
			try:
				tle_rec.compute(arg)
				# if arg:
				# 	tle_rec.compute(arg)
				# else:
				# 	tle_rec.compute()
				d = ephem.Date(s['date'].replace('T',' ')).datetime().date()
				if d < too_old:
					raise RuntimeError
				long_lat = [self._radians_to_degrees(tle_rec.sublong), self._radians_to_degrees(tle_rec.sublat)]
				elevation = tle_rec.elevation
				eclipsed = tle_rec.eclipsed
				try:
					surface_diameter = self._surface_diameter(elevation, horizon_angle)	## defaults to zero, but changable in call
				except:
					surface_diameter = None
				if surface_diameter != None:
					surface_diameters = {0: surface_diameter}
					for angle in [10, 15, 20, 25]:
						surface_diameters[angle] = self._surface_diameter(elevation, angle)	## always return more values - realistic viewing angle
				else:
					surface_diameters = None
			except RuntimeError:
				# RuntimeError: cannot compute the body's position at YYYY/MM/DD HH:MM:SS
				if self._debug:
					sys.stderr.write('now(): RuntimeError: %s\n' % (s['name']))
				s['valid'] = False
				long_lat = [None, None]
				elevation = None
				eclipsed = None
				surface_diameter = None
				surface_diameters = None
			except ValueError:
				# ValueError: TLE elements are valid for a few weeks around their epoch, but you are asking about a date ### days from the epoch
				if self._debug:
					sys.stderr.write('now(): ValueError %s\n' % (s['name']))
				s['valid'] = False
				long_lat = [None, None]
				elevation = None
				eclipsed = None
				surface_diameter = None
				surface_diameters = None
			sat_name = s['name']
			sat_id = s['id']
			tle = [s['line1'], s['line2']]
			if sat_name not in r['directory']:
				r['directory'][sat_name] = [sat_id]
			else:
				r['directory'][sat_name].append(sat_id)
			r['satellite'][sat_id] = {
					'id': sat_id,
					'name': sat_name,
					'long_lat': long_lat,
					'elevation': elevation,
					'eclipsed': eclipsed,
					'tle': tle
					}
			if surface_diameter != None:
				r['satellite'][sat_id]['surface_diameter'] = surface_diameter
			if surface_diameters != None:
				r['satellite'][sat_id]['surface_diameters'] = surface_diameters.copy()
			if s['valid']:
				r['satellite'][sat_id]['position'] = {}
								# https://rhodesmill.org/pyephem/radec.html
								# Astrometric Topocentric Position for the epoch of your Observer
								#'a_ra': self._hours_to_degrees(tle_rec.a_ra),
								#'a_dec': self._radians_to_degrees(tle_rec.a_dec),
								# Apparent Geocentric Position for the epoch-of-date
								#'g_ra': self._hours_to_degrees(tle_rec.g_ra),
								#'g_dec': self._radians_to_degrees(tle_rec.g_dec),
								# Apparent Topocentric Position for the epoch-of-date
								# 'ra': self._hours_to_degrees(tle_rec.ra),
								# 'dec': self._radians_to_degrees(tle_rec.dec)

				try:
					r['satellite'][sat_id]['position']['ra_dec'] = [tle_rec.ra, tle_rec.dec]
				except:
					pass

				try:
					r['satellite'][sat_id]['position']['alt_az'] = [self._radians_to_degrees(tle_rec.alt), self._radians_to_degrees(tle_rec.az)]
				except:
					pass

				try:
					r['satellite'][sat_id]['mag'] = tle_rec.mag
				except:
					r['satellite'][sat_id]['mag'] = None
				try:
					r['satellite'][sat_id]['neverup'] = tle_rec.neverup
				except:
					r['satellite'][sat_id]['neverup'] = None
					pass

			## now compute long_lat for +/- some number of seconds
			if historical and s['valid']:
				r['satellite'][sat_id]['historical'] = {}
				for timeskew in range(-480, 480+1, 30):
					skew = when + datetime.timedelta(seconds=timeskew)
					if self._observer:
						# All in UTC - https://rhodesmill.org/pyephem/date.html-
						self._observer.date = skew
						arg = self._observer
					else:
						arg = skew
					try:
						tle_rec.compute(arg)
						long_lat = [self._radians_to_degrees(tle_rec.sublong), self._radians_to_degrees(tle_rec.sublat)]
						elevation = tle_rec.elevation
						eclipsed = tle_rec.eclipsed
					except RuntimeError:
						# RuntimeError: cannot compute the body's position at YYYY/MM/DD HH:MM:SS
						if self._debug:
							sys.stderr.write('now(): RuntimeError: %s\n' % (s['name']))
						long_lat = [None, None]
						elevation = None
						eclipsed = None
					except ValueError:
						# ValueError: TLE elements are valid for a few weeks around their epoch, but you are asking about a date ### days from the epoch
						if self._debug:
							sys.stderr.write('now(): ValueError %s\n' % (s['name']))
						long_lat = [None, None]
						elevation = None
						eclipsed = None
					s_skew = skew.replace(microsecond=0).isoformat(timespec='seconds')
					r['satellite'][sat_id]['historical'][s_skew] = {'long_lat': long_lat.copy(), 'elevation': elevation, 'eclipsed': eclipsed}
					try:
						r['satellite'][sat_id]['historical'][s_skew]['ra_dec'] = [tle_rec.ra, tle_rec.dec]
					except:
						pass
					try:
						r['satellite'][sat_id]['historical'][s_skew]['alt_az'] = [self._radians_to_degrees(tle_rec.alt), self._radians_to_degrees(tle_rec.az)]
					except:
						pass

		return r

	def _au_to_km(self, a):
		""" I think in Km's - even if computers think in AU's """
		return a * 149598073.0

	def _hours_to_degrees(self, h):
		""" I think in degress - even if computers think in radians """
		return h * (360/15)

	def _radians_to_degrees(self, d):
		""" I think in degress - even if computers think in radians """
		return d * (180/math.pi)

	def _degrees_to_radians(self, a):
		""" I think in degress - even if computers think in radians """
		return a / (180/math.pi)

	def _surface_diameter(self, elevation, horizon_angle=0.0):
		""" calculate how much of the surface is seen """

		earth_radius = ephem.earth_radius
		earth_circumference = 2 * math.pi * earth_radius

		if horizon_angle == 0:
			# caculate the angle of the satellites view horizon
			angle = math.acos(earth_radius / (earth_radius + elevation))
		else:
			# more complex
			if horizon_angle < 0 or horizon_angle > 90:
				raise ValueError
			# Michael's math - glad someone was paying attention in class!
			theta = math.asin(math.sin(math.radians(90.0 + horizon_angle)) * earth_radius / (earth_radius + elevation))
			angle = math.radians(180.0 - 90.0 - horizon_angle) - theta
			if angle < 0.0:
				angle = 0.0

		# caculate the diameter on the earths surface
		percentage = (2 * angle) / (2 * math.pi)
		surface_diameter = earth_circumference * percentage
		# this is the distance (diameter) of the visible part of the earth surface (in meters)
		return surface_diameter

# https://www.techbeamers.com/python-float-range/
import decimal
def float_range(start, stop, step):
	while start < stop:
		yield float(start)
		start += decimal.Decimal(step)

def main():

	m = EarthMars()
	# Santa Cruz (36.9812°N, 122.0262°W)
	m.earthobserver(lon=-122.0262, lat=36.9812, amsl=40.0)
	mars = m.mars()

	now = datetime.datetime.utcnow()

	if False:
		all_times = {}
		all_times['previous_rising'] = mars['previous_rising'].replace(microsecond=0)
		all_times['previous_setting'] = mars['previous_setting'].replace(microsecond=0)
		all_times['now'] = now.replace(microsecond=0)
		all_times['next_rising'] = mars['next_rising'].replace(microsecond=0)
		all_times['next_setting'] = mars['next_setting'].replace(microsecond=0)

		for r in sorted(all_times, key=lambda y: all_times[y]):
			print("%s UTC ; %s" % (all_times[r], r))

	now = now.replace(hour=0, minute=0, second=0, microsecond=0)
	for timeskew in float_range(-27*4, 27*4*10+1, 0.25):
		delta_when = now + datetime.timedelta(weeks=timeskew)
		mars = m.mars(when=delta_when)

		mars_km = mars['earth_distance_km']
		secs = mars_km / (ephem.c / 1000.0)
		delay = datetime.datetime.fromtimestamp(secs)	## XXX maybe should not use datetime as a way of getting min/sec printed

		sun = m.sun(when=delta_when)
		sun_km = sun['earth_distance_km']

		occult_degrees = 4.0 * (sun['radius'] * 2)

		if abs(mars['ra_dec'][0] - sun['ra_dec'][0]) <= occult_degrees and abs(mars['ra_dec'][1] - sun['ra_dec'][1]) <= occult_degrees:
			if mars_km > sun_km:
				sun_mars = "M"
			else:
				sun_mars = "S"
		else:
			sun_mars = "-"

		print("%s %12.1f %5s | [%8.2f %8.2f] | [%8.2f %8.2f] | %s | %8.4f < %8.4f ; %6.1f" % (
				delta_when, mars_km, delay.strftime("%M:%S"),
				mars['ra_dec'][0], mars['ra_dec'][1],
				sun['ra_dec'][0], sun['ra_dec'][1],
				sun_mars,
				mars['radius'],
				sun['radius'],
				timeskew
				)
		)

if __name__ == '__main__':
	main()
