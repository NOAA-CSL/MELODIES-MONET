---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.6
kernelspec:
  display_name: Python 3 (ipykernel)
  language: python
  name: python3
---

# Evaluate AirNow and WRF-Chem

Our first example will demonstrate the basics available in MELODIES MONET to 
compare WRF-Chem model results against AirNow surface observations for ozone, 
carbon monoxide, nitrogen oxide, nitrogen dioxide, and PM{subscript}`2.5`.

Please download the WRF-Chem example model dataset from TBA.

```{code-cell} ipython3
import sys; sys.path.append("../../")
# ^ sys.path-hacking needed until the package is installable

from melodies_monet import driver
```

## Driver class

Now lets create an instance of the analysis driver class, {class}`melodies_monet.driver.analysis`.
It consists of these main parts:
* model instances
* observation instances
* a paired instance of both

This will allow us to move things around the plotting function for
spatial and overlays and more complex plots.

```{code-cell} ipython3
an = driver.analysis()
```

## Control file

Now we set the YAML control file and begin by reading the file.

```{note}
Check out the {doc}`/appendix/yaml` for info on how to create
and modify these files.
```

```{code-cell} ipython3
:tags: [hide-output]

an.control = "control_wrfchem_mech-0905_2.yaml"
an.read_control()
an.control_dict
```

````{admonition} Note: This is the file that was loaded.
:class: dropdown

```{literalinclude} control_wrfchem_mech-0905_2.yaml
:caption:
:linenos:
```
````

## Load the model data

The driver will automatically loop through the "models" found in the `model` section
of the YAML file and create an instance of {class}`melodies_monet.driver.model` for each
that includes the
* label
* mapping information
* file names (can be expressed using a glob expression)
* xarray object

```{code-cell} ipython3
an.open_models()
```

Applying {meth}`~melodies_monet.driver.analysis.open_models`
populates the {attr}`~melodies_monet.driver.analysis.models` attribute.

```{code-cell} ipython3
an.models
```

We can access the underlying dataset with the
{attr}`~melodies_monet.driver.model.obj` attribute.

```{code-cell} ipython3
an.models['RACM_ESRL'].obj
```

## Load the observational data

As with the model data, the driver will loop through the "observations" found in
the `obs` section of the YAML file and create an instance of
{class}`melodies_monet.driver.observation` for each.

```{code-cell} ipython3
an.open_obs()
```

```{code-cell} ipython3
an.obs
```

```{code-cell} ipython3
an.obs['airnow'].obj
```

## Pair model and observational data

```{code-cell} ipython3
:tags: [hide-output]

an.pair_data()
```

```{code-cell} ipython3
an.paired
```

```{code-cell} ipython3
an.paired['airnow_RACM_ESRL']
```

## Plot

```{code-cell} ipython3
an.plotting()
```

## Statistics

```{code-cell} ipython3
an.stats()
```

The stats routine has produced six files, one of which is:
```{literalinclude} output/stats.OZONE.all.CONUS.2019-09-05_06.2019-09-06_06.csv
:caption:
```
