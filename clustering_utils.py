import numpy as np 
import musclebeachtools_hlab.musclebeachtools as mbt
import os
import glob
import shutil
import pandas as pd
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
# import xlsxwriter
from io import BytesIO
from matplotlib.backends.backend_pdf import PdfPages
from pandas.plotting import table
import six
import seaborn as sns

h_neurons_clust_out = glob.glob('*/co/*neurons_group0.npy')

qualed_neurons = sorted(glob.glob('*/co/*qualed.npy'))

def render_mpl_table(data, col_width=3.0, row_height=0.625, font_size=14,
                     header_color='#40466e', row_colors=['#f1f1f2', 'w'], edge_color='w',
                     bbox=[0, 0, 1, 1], header_columns=0,
                     ax=None, **kwargs):
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=size)
        ax.axis('off')

    mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns, **kwargs)

    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)

    for k, cell in  six.iteritems(mpl_table._cells):
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        else:
            cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
    return fig, ax

def qual_cells(h_neurons_clust_out):
    '''
        goes through and quals each of the neuron lists given to it

        h_neurons_clust_out: paths to each file, make sure the paths are complete 
        or you're in a parent directory

        returns: nothing, just creates the qualed files in each directory
    '''
    for h_n in h_neurons_clust_out:
        cells = np.load(h_n,allow_pickle=True)
        for cell in cells:
            cell.remove_large_amplitude_spikes(4,lplot=False)
        mbt.autoqual(cells,'/media/HlabShare/models/xgb_model')
        front_letter = h_n.rfind('H_')
        end_letter = h_n.rfind('.npy')
        np.save(h_n[:end_letter] + 'qualed',cells)



def determine_shape(num_cells):
    if num_cells>15:
        if num_cells < 30:
            rows=3
            cols=10
        else:
            rows=int(np.ceil((num_cells/15)))
            cols=15
    else:
        cols=5
        rows=3
    
    return rows, cols

def plot_wfs(cells, quality, title_string):
    fig = plt.figure(constrained_layout=False, figsize=(11,8))
    fig.suptitle(f'Quality {quality}: {len(cells)} {title_string}')
    num_rows, num_cols = determine_shape(len(cells))
    print(f'num_cells: {len(cells)} rows: {num_rows} cols: {num_cols}')
    gs = fig.add_gridspec(ncols=num_cols,nrows=num_rows)

    print(f'plotting_{quality}s ....')
    for idx in range(len(cells)):

        ax = fig.add_subplot(gs[int(idx/num_cols), int(idx%num_cols) ])
        ax.set_title(f"CI: {cells[idx].clust_idx}")
        sns.despine(left=True, bottom=True)
        ax.get_xaxis().set_visible(False)
        ax.get_yaxis().set_visible(False)
        ax.plot(cells[idx].waveforms)

    plt.tight_layout(rect=[0,0,1,.95])

    return fig


def make_pdf(qualed_cell_paths, testing_params=True):
    
    if testing_params:
        pdf_title='evaluating_parameters.pdf'
        column_titles = ['thresh','cpz','probe','filt','total',
                            'ones','twos','threes','fours']

    else:
        pdf_title='evaluating_clustering_output.pdf'
        column_titles = ['total', 'ones','twos','threes','fours']


    big_panda = pd.DataFrame(columns = column_titles)   

    with PdfPages(pdf_title) as pdf: 
        for qn in qualed_cell_paths:

            cells = np.load(qn,allow_pickle=True)
            total = len(cells)


            if testing_params:

                first_thresh = qn.rfind('thresh')
                thresh = int(qn[first_thresh+6:first_thresh+7])
                first_czp = qn.rfind('czp')
                czp = float(qn[first_czp+3:first_czp+7])
                first_p = qn.rfind('5p')
                probe = int(qn[first_p+2:first_p+3])
                first_filt = qn.rfind('_1_')
                filt = 'average'

                title_string = f'Thresh: {thresh} CZP: {czp} TOTAL: {total}'
            else:
                title_string = f'TOTAL: {total}'


            ones = [cell for cell in cells if cell.quality==1]
            twos = [cell for cell in cells if cell.quality==2]
            threes = [cell for cell in cells if cell.quality==3]
            fours = [cell for cell in cells if cell.quality==4]

            all_len = [len(ones),len(twos),
                        len(threes),len(fours)]
                    
            max_len = max(all_len)

            imgdata_one = BytesIO()
            imgdata_two = BytesIO()
            imgdata_three = BytesIO()
            imgdata_four = BytesIO()
            
            tabledata = BytesIO()
            
        
            #FIGURE FOR ONES
            one_fig = plot_wfs(ones, 1, title_string)
            pdf.savefig()
            plt.close()

            #FIGURE FOR TWOS
            two_fig = plot_wfs(twos, 2, title_string)
            pdf.savefig()
            plt.close()

            #FIGURE FOR THREES
            three_fig = plot_wfs(threes, 3, title_string)
            pdf.savefig()
            plt.close()

            
            #FIGURE FOR FOURS
            four_fig = plot_wfs(fours, 4, title_string)
            pdf.savefig()
            plt.close()

            if testing_params:
                small_panda = pd.DataFrame({'thresh':[thresh,
                                                    np.NaN],
                                            'cpz':[czp,
                                                    np.NaN,],
                                            'probe':[probe,
                                                    np.NaN],
                                            'filt':[filt,
                                                    np.NaN],
                                            'total':[total,
                                                    'percentages:'],
                                            'ones':[len(ones),
                                                    round(len(ones)*100/total,1)],
                                            'twos':[len(twos),
                                                    round(len(twos)*100/total,1)],
                                            'threes':[len(threes),
                                                    round(len(threes)*100/total,1)],
                                            'fours':[len(fours),
                                                    round(len(fours)*100/total,1)]})
            else:
                small_panda = pd.DataFrame({
                                            'total':[total,
                                                    'percentages:'],
                                            'ones':[len(ones),
                                                    round(len(ones)*100/total,1)],
                                            'twos':[len(twos),
                                                    round(len(twos)*100/total,1)],
                                            'threes':[len(threes),
                                                    round(len(threes)*100/total,1)],
                                            'fours':[len(fours),
                                                    round(len(fours)*100/total,1)]})

            small_panda_fig,small_panda_ax = render_mpl_table(small_panda,col_width=2)
            pdf.savefig()
            plt.close()
            # small_panda_fig.savefig(tabledata, format='png')
            big_panda = big_panda.append(small_panda)
        big_panda_fig,big_panda_ax = render_mpl_table(big_panda,col_width=2)
        pdf.savefig()
        plt.close()
        big_panda.to_pickle('caf22_average_big_panda.pkl')

