import matplotlib.pyplot as plt
def plot_style():
   
    params = {'legend.fontsize': 20,
             'axes.labelsize': 20,
             'axes.titlesize':20,
             'xtick.labelsize':20,
             'ytick.labelsize':20,
             'xtick.major.size':5,
              'xtick.minor.size':2,
             'ytick.major.size':5,
              'ytick.minor.size':2,
             'figure.facecolor':'w',
             #'lines.linewidth' : 1.5,
              'xtick.major.width':1.0,
              'ytick.major.width':1.0,
              'xtick.minor.width':1.0,
              'ytick.minor.width':1.0,
              'axes.linewidth':1.5,
              'xtick.direction':'in',
              'ytick.direction':'in',
             'ytick.labelleft':True,
              'text.usetex' : False,
             'font.family': 'sans-serif'}
  
    plt.rcParams.update(params)