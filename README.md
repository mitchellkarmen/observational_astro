# observational_astro
some utils for class

## Requirements:
- pandas
- astropy
- astroplan
- astroquery

## How to use:
 `cd` into the directory where you cloned this.  you can now open python and `import exoplanet_utils as exutil`.
 Then, the most useful fn would be `night1_planets = exutil.get_best_planets(night=1)`.  If you automatically want to make a plot of those, you can specify the argument `plot_filepath='path/to/my/plot.png'`.  Otherwise, you can pick through the targets for your favorites, and then `exutil.make_plot(my_targets, night=2, filepath='path/to/my/plot.png'`.  To load all possible exoplanets with relevant (and irrelevant) info, `exutil.get_exoplanets()`.