import matplotlib.pyplot as plt
def plot_style():
   
    params = {'legend.fontsize': 14,
             'axes.labelsize': 22,
             'axes.titlesize':20,
              'axes.titlepad':20,
             'xtick.labelsize':22,
             'ytick.labelsize':22,
             'xtick.major.size':5,
             'xtick.minor.size':2,
             'ytick.major.size':5,
             'ytick.minor.size':2,
             'figure.facecolor':'w',
             'lines.linewidth' : 1.5,
             'xtick.major.width':1.0,
             'ytick.major.width':1.0,
             'xtick.minor.width':1.0,
             'ytick.minor.width':1.0,
             'axes.linewidth':1.5,
             'xtick.direction':'in',
             'ytick.direction':'in',
              'xtick.major.pad':6,
              'ytick.major.pad':6,
             'ytick.labelleft':True,
             'text.usetex' : False,
              'lines.linewidth':2,
             'font.family': 'sans-serif'}
  
    plt.rcParams.update(params)
