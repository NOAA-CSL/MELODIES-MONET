In order to contribute to the readthedocs.io documents.

Download two conda packages:
conda install sphinx
conda install sphinx_rtd_theme

Edit the desired content.

Make sure you are in the docs folder
Create the html code by typing: make html
May need to use: 'make clean' prior to 'make html' to have all changes implenented.

Open the index.html file in _build/html using your favorite browser.

Once you are happy with your updates. Create a pull request on github 
with your changes. Once the changes are on github, readthedocs.io will
automatically include the updates online.

Note to make a pdf, you should be able to use: make latexpdf
