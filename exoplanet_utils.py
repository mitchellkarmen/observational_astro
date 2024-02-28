import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
from astropy import units as u
from astroplan.plots import plot_airmass,plot_sky
from astroplan import Observer, FixedTarget, time_grid_from_range, AirmassConstraint,AtNightConstraint,AltitudeConstraint,observability_table,SunSeparationConstraint,MoonSeparationConstraint, is_observable
from astropy.coordinates import EarthLocation, get_body, SkyCoord, Angle
from astroquery.simbad import Simbad
import datetime as dt
from datetime import timezone

apo_coords = ('105.8200d', '32.7805d') # longitude, lattitude
apo = Observer.at_site(('apache point observatory'), timezone = 'US/Mountain')

night_1_midnight = Time('2024-03-12 07:00:00', scale='utc', location=apo_coords)
night_2_midnight = Time('2024-03-13 07:00:00', scale='utc', location=apo_coords)
night_3_midnight = Time('2024-03-14 07:00:00', scale='utc', location=apo_coords)
midnights = [
	night_1_midnight,
	night_2_midnight,
	night_3_midnight
]

night1_sunset = apo.sun_set_time(night_1_midnight, which='nearest')
night2_sunset = apo.sun_set_time(night_2_midnight, which='nearest')
night3_sunset = apo.sun_set_time(night_3_midnight, which='nearest')
sunsets = [
	night1_sunset,
	night2_sunset,
	night3_sunset
]

night1_times = Time([night1_sunset.iso, "2024-03-12 07:00:00"])
night2_times = Time([night2_sunset.iso, "2024-03-13 07:00:00"])
night3_times = Time([night3_sunset.iso, "2024-03-14 07:00:00"])
obs_times = [
	night1_times,
	night2_times,
	night3_times
] # what time we are observing

default_constraints = [
	AltitudeConstraint((32.7805-45)*u.deg, (32.7805 + 45)*u.deg),
	AirmassConstraint(2),
	AtNightConstraint.twilight_civil(),
	MoonSeparationConstraint(30 * u.deg)
]

def get_exoplanets():
	"""
		returns a pandas dataframe of exoplanets, requires csv of planets.
		planets['coordinates'] are astropy SkyCoord of ra, dec for astroplan input
		planets['night{i}_transit'] is the datetime of the transit nearest midnight
		planets['night{i}_transit_delay'] is number of hours *after* midnight of nearest transit.  negative for before midnight transit.
		(we want this to be -6<x<0)

	"""
	sources = pd.read_csv('exoplanets.csv', skiprows=55)
	sources['coordinates']  = [SkyCoord(s['ra']*u.deg, s['dec']*u.deg, frame='icrs') for i, s in sources.iterrows()]

	for i, night in enumerate([night_1_midnight, night_2_midnight, night_3_midnight]):
	    time_diff = night.jd - sources['pl_tranmid'].values # days between observing night and transit midpoint date
	    n_transits = time_diff / sources['pl_orbper'].values # number of transits between observing night and transit midpoint date
	    nearest_transit_jd = (np.round(n_transits) * sources['pl_orbper'].values) + sources['pl_tranmid'].values # jd of nearest transit
	    transit_delay = (nearest_transit_jd - night.jd)*24 # hours after midnight of nearest transit
	    nearest_transit_jd[np.isnan(nearest_transit_jd)] = 0

	    nearest_transit_time = [Time(t, format='jd').iso for t in nearest_transit_jd]
	    sources[f'night{i+1}_transit'] = nearest_transit_time # iso date of nearest transit
	    sources[f'night{i+1}_transit_delay'] = transit_delay # hours after midnight of nearest transit

	return sources

def make_plot(targets, night, filepath=None):
	"""
	input:
		targets (pd.DataFrame): pandas dataframe of targets.  should have columns corresponding to given night

	returns figure
	"""
	ii=0
	i_night = night-1
	fig = plt.figure(dpi = 200)
	for i, target in targets.iterrows():
	    p = plot_airmass(targets=target['coordinates'],
	                 observer=apo,
	                 time=midnights[i_night],
	                 use_local_tz=False
	                 )
	    transit_airmass = apo.altaz(target[f'night{night}_transit'], target['coordinates']).secz
	    transit_time_utc = Time(target[f'night{night}_transit'], format='iso') #+ 7*u.hr
	    if (transit_airmass > 0) & (transit_airmass < 5):
	        plt.plot_date(transit_time_utc.datetime, transit_airmass, color=f'C{ii}')
	        plt.text(transit_time_utc.datetime, transit_airmass, target['pl_name'], fontsize=10)
	    ii+=1
	plt.xlim(obs_times[i_night][0].datetime, obs_times[i_night][1].datetime)
	plt.ylim(3, 0.98)
	if filepath:
		plt.savefig(filepath, dpi = 300)

	return fig


def get_best_planets(night, constraints=default_constraints, plotting=True):
	"""
	retrieves a list of planets with visible transits on a given observing night
	input:
		night (int): night 1, 2, or 3.  corresponds to night of March 11th, 12th, 13th.
		constraints (tuple): list of astroplan constraints.  if in doubt, use the ones above
		plotting (bool): whether or not to make a plot

	returns:
		targets (dataframe)
	"""
	i_night = night - 1
	all_planets = get_exoplanets()
	target_list = [FixedTarget(s['coordinates'], name=s['pl_name']) for i, s in all_planets.iterrows()]

	observable  = is_observable(constraints, apo, target_list, time_range=obs_times[i_night])
	targets  = all_planets.iloc[np.where(observable &
                                    (all_planets[f'night{night}_transit_delay'].values > -6) & # transit between 6pm and midnight
                                    (all_planets[f'night{night}_transit_delay'].values < 0))]


	if plotting:
		make_plot(targets, night=1)
		plt.show()

	return targets
