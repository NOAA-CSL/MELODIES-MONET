======================
Command Line Interface
======================

Installing :mod:`melodies_monet` provides a command line interface (CLI):
``melodies-monet``.

For example, you can use the CLI to run control files without writing
any Python code::

    melodies-monet run control.yml

**Subcommands**

* |run|_ -- run a control file
* |get-airnow|_ -- get AirNow data
* |get-aeronet|_ -- get AERONET data

.. |run| replace:: ``run``
.. _run: #melodies-monet-run

.. |get-airnow| replace:: ``get-airnow``
.. _get-airnow: #melodies-monet-get-airnow

.. |get-aeronet| replace:: ``get-aeronet``
.. _get-aeronet: #melodies-monet-get-aeronet


.. click:: melodies_monet._cli:_typer_click_object
   :prog: melodies-monet
   :nested: full
