Developer's Guide
=================

This whole page is under development.

Description of Branches
-----------------------

There are several primary repository branches
for specific development tracks.

main
____
This is the stable release branch.
It is updated from develop and tagged prior to each release.
The melodies-monet conda package is created from main.

develop
_______
This is the parent development branch in which
to consolidate the various development tracks.

develop_plots
_____________

develop_aircraft
________________

develop_satellite
_________________

develop_testsuite
_________________



How to incorporate updates to MELODIES MONET
--------------------------------------------


(a) Fork the GitHub repository to your own GitHub account:

    https://github.com/NOAA-CSL/MELODIES-MONET

    Note: you can pull updates from the main NOAA repository by using the “Fetch upstream” button on your fork. Alternatively::

    $ git pull upstream master
    $ git push origin master

(b) Navigate on cheyenne/casper to where you would like to keep the MELODIES MONET code (e.g. in your work location) and clone your fork to cheyenne::

    $ git clone https://github.com/$GitHubUsername/$ForkName.git

(c) Checkout the develop_plots branch - you need to do this with the remote branch as well as create a local tracking branch::

    $ git checkout origin/develop_plots
    $ git checkout develop_plots

    Then all development work will be in the monet_analysis folder::

    $ cd monet_analysis


User Guide Development
----------------------

If you add a component to MELODIES MONET, please follow the instructions below 
to update the readthedocs user guide. 

   * Add instructions.
   
Please see the `Documentation <https://github.com/NOAA-CSL/MELODIES-MONET/projects/2>`_ 
project on GitHub to learn about current and future development.   

Adding New Datasets
-------------------