
fig_ax_spine_width = 1.5
figure_font_size = 16
label_font_size = 11
num2char="abcdefghijk"

def get_standard_plt():

    import matplotlib
    import matplotlib.pyplot as plt

    matplotlib.rcParams['font.size'] = figure_font_size
    matplotlib.rcParams['legend.handlelength'] = 1
    plt.rcParams.update({'axes.labelsize': 'large'})
    return plt

def init_subplots(plt,ncols=1,nrows=1,fig_width=8,wid_high_ratio=0.63,more_options={},**kwargs):
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(fig_width * ncols,  wid_high_ratio * fig_width * nrows),constrained_layout=True,**kwargs)

    #TODO Set spine width
    #TODO Set label font size
    argkeys=more_options.keys()
    if "auto_title" in argkeys:
        title_font_size=30
        if more_options["auto_title"]:
            for i in range(len(axs)):
                for j in range(len(axs[i])):
                    index = j*nrows+i
                    axs[i][j].set_title("({})".format(num2char[index]),x=-0.12,y=1.15,size=title_font_size,pad=-title_font_size*2)


    if ncols==1 and nrows==1:
        return fig, [axs]
    else:
        return fig, axs