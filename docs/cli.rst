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
* |get-aqs|_ -- get AQS data
* |get-ish|_ -- get ISH data
* |get-ish-lite|_ -- get ISH-Lite data

.. |run| replace:: ``run``
.. _run: #melodies-monet-run

.. |get-airnow| replace:: ``get-airnow``
.. _get-airnow: #melodies-monet-get-airnow

.. |get-aeronet| replace:: ``get-aeronet``
.. _get-aeronet: #melodies-monet-get-aeronet

.. |get-aqs| replace:: ``get-aqs``
.. _get-aqs: #melodies-monet-get-aqs

.. |get-ish| replace:: ``get-ish``
.. _get-ish: #melodies-monet-get-ish

.. |get-ish-lite| replace:: ``get-ish-lite``
.. _get-ish-lite: #melodies-monet-get-ish-lite

.. click:: melodies_monet._cli:_typer_click_object
   :prog: melodies-monet
   :nested: full
