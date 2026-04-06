import matplotlib.colors

def nb_colors():
    nonbinary_colors = [
        [0.0, "#000000"],  # Black 
        [0.33, "#9C59D1"], # Purple
        [0.66, "#FFFFFF"], # White
        [1.0, "#FFF430"],  # Yellow 
    ]

    return nonbinary_colors

def ace_colors(reverse=True):
    if reverse==False:
        ace_colors = [
            [0.0, "#000000"],  # Black 
            [0.33, "#A3A3A3"], # Purple
            [0.66, "#FFFFFF"], # White
            [1.0, "#800080"],  # Yellow 
        ]
    else:
        ace_colors = [
            [0.0, "#800080"],  # Black 
            [0.33, "#FFFFFF"], # Purple
            [0.66, "#A3A3A3"], # White
            [1.0, "#000000"],  # Yellow 
        ]
    return ace_colors

def mpl_ace():

    return matplotlib.colors.LinearSegmentedColormap.from_list("", ["#000000", "#A3A3A3", "#FFFFFF", "#800080"])

def mpl_nb():

    return matplotlib.colors.LinearSegmentedColormap.from_list("", ["#000000", "#9C59D1", "#FFFFFF", "#FFF430"])