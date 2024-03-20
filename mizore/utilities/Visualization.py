def draw_x_y_line_relation(x_data, y_data, x_label, y_label, y_data2=None, y_label2=None, ylim=None, ylim2=None,
                           red_line=None, filename="Untitled"):
    import matplotlib
    import matplotlib.pyplot as plt
    _fig, ax = plt.subplots()
    plt.xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.plot(x_data, y_data, marker="o")

    if y_data2 != None:
        ax2 = ax.twinx()
        ax2.plot(x_data, y_data2, 'go-')
        if y_label2 != None:
            ax2.set_ylabel(y_label2)
        if ylim2 != None:
            ax2.set_ylim(ylim2)

    if ylim != None:
        ax.set_ylim(ylim)
    if red_line != None:
        ax.hlines(red_line, 0, x_data[len(x_data) - 1], colors="r", linestyles="dashed")

    if filename == None:
        plt.show()
    else:
        plt.savefig(filename + ".png", bbox_inches='tight')

    plt.close()
