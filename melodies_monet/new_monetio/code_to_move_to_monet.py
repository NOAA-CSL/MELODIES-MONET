#These are functions in MONET that we are working on updating and then will move over to MONET.

def savefig(fname, loc=1, decorate=True, height=50.0, **kwargs):
    """Save figure and add the MELODIES-MONET logo.
    Parameters
    ----------
    fname : str
        output file name.
    loc : int
        the location for the monet logo.
    decorate : bool
        Include logo on plot (True) or not (False)
    **kwargs : dict
        kwargs for the matplotlib.pyplot.savefig function.
    Returns
    -------
    plots
        Saved plots optionally with MELODIES-MONET logo included
    """
    import io
    import os
    import sys
    from PIL import Image
    import matplotlib.pyplot as plt
    try:
        from pydecorate import DecoratorAGG
        pydecorate = True
    except ImportError:
        pydecorate = False
    plt.savefig(fname, **kwargs)
    if pydecorate and decorate:
        img = Image.open(fname)
        dc = DecoratorAGG(img)
        if loc == 1:
            dc.align_bottom()
        elif loc == 2:
            dc.align_bottom()
            dc.align_right()
        elif loc == 3:
            dc.align_right()
        elif loc == 4:
            dc.align_bottom()
            dc.align_left()
        # sys.argv[0])[-5] + 'data/MONET_logo.png'
        # print(os.path.basename(__file__))
        #logo = os.path.abspath(__file__)[:-17] + 'data/MONET-logo.png'
        logo = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data/MM_logo.png'))
        dc.add_logo(logo,height=height)
        if fname.split('.')[-1] == 'png':
            img.save(fname, "PNG")
        elif fname.split('.')[-1] == 'jpg':
            img.save(fname, "JPEG")
