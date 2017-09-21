# -*-coding:Utf-8 -*
# ----------------------------------------------------------------------
# Align data without model
# ----------------------------------------------------------------------
#   External libraries
# ----------------------------------------------------------------------
import sys
import os
# Full path of this file
full_path_here = os.path.realpath(__file__)
text = full_path_here.split('/')
a = ''
i = 0
while i < len(text)-1:
   a = a + text[i] + '/'
   i = i + 1
full_path = a

#filename = full_path + '../' + '.pythonexternallibpath'
#file = open(filename, 'r')
#for line in file:
#    path_lib_ext=line
#file.close()
#if path_lib_ext != 'None':
#    sys.path.insert(0, path_lib_ext[:-1])
# ----------------------------------------------------------------------
#    Packages
# ----------------------------------------------------------------------
import os
import sys
import glob
import numpy as np
import pandas as pd
from scipy import stats
import ConfigParser as cp
import bokeh.layouts as blyt
import statsmodels.api as sm
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import bokeh.plotting as bplt
from bokeh.models import HoverTool, TapTool, ColumnDataSource, OpenURL
from bokeh.models.widgets import DateFormatter, NumberFormatter, DataTable, TableColumn
import bokeh.io as io
#
# ======================================================================
#   To align the data without model.
# ======================================================================
def help():
    text = "Non parametric alignement of the data."
    return text

def communicate(cfg, verbose, text):
    if cfg.getint('Modelling', 'Verbose') >= verbose:
        print text
# ----------------------------------------------------------------------
def plot(cfgsetup=False, models=False, model_param=False, time_serie=False,\
         obs_properties=False, options=False, interpol_method=False):
    """
    To align the data without model.
    """
    # ----------------------------------------------------------------------
    #   Arnaud's method
    # ----------------------------------------------------------------------
    def align(photometryData, tel_ref, obs_name, cfg):

        cond = (photometryData['reject'] == 0) &\
               (photometryData['obs'] == tel_ref)

        ref_date = photometryData[cond]['dates'].values
        ref_mag = photometryData[cond]['magnitude'].values
        ref_flux = np.power(10,-ref_mag/2.5)

        # Read lc to calibrate
        cond = (photometryData['reject'] == 0) &\
               (photometryData['obs'] == obs_name)
        obs_date = photometryData[cond]['dates'].values
        obs_mag = photometryData[cond]['magnitude'].values
        obs_flux = np.power(10,-obs_mag/2.5)

        # Plot LC
    #    fig1 = plt.figure(1)
    #    plt.subplot(211)
    #    plt.plot(ref_date,ref_flux,'+r')
    #    plt.subplot(212)
    #    plt.plot(obs_date,obs_flux,'+g')
    #    plt.savefig("LC.pdf")
    #    plt.close(1)

        ## Interpolation inter-lc
        datmin = np.max([np.min(ref_date), np.min(obs_date)])
        datmax = np.min([np.max(ref_date), np.max(obs_date)])
    #    ref_ind = np.where((ref_date >= datmin) & (ref_date <= datmax))
        obs_ind = np.where((obs_date >= datmin) & (obs_date <= datmax))
    #    ref_fselec = ref_flux[ref_ind]
        obs_fselec = obs_flux[obs_ind]
        ##
    #    t_com = np.concatenate([ref_date[ref_ind],obs_date[obs_ind]])
        t_com = obs_date[obs_ind]
        ##
        ref_int = interp1d(ref_date,ref_flux)
        ref_fint = ref_int(t_com)
    #    obs_int = interp1d(obs_date,obs_flux)
    #    obs_fint = obs_int(t_com)
        if len(t_com) < 3:
            text = "     %3d data: WARNING no data plotted" % len(t_com)
            communicate(cfg, 1, text)
            return [ref_date,ref_flux,obs_date,np.zeros_like(ref_fint)]
        else:
            text = "     %3d data" % len(t_com)
            communicate(cfg, 1, text)
        ## Plot LC interpol
    #    fig1 = plt.figure(3)
    #    plt.subplot(211)
    #    plt.title('Reference interpolee')
    #    plt.plot(t_com,ref_fint,'+r')
    #    plt.subplot(212)
    #    plt.title('Observatoire interpole')
    #    plt.plot(t_com,obs_fint,'+g')
    #    plt.savefig("LCinterp.pdf")
    #    plt.close(3)

        ## Calcul des regressions
        fig2 = plt.figure(2)
        xreg = np.linspace(np.amin(ref_fint),np.amax(ref_fint),2)
        ##
        ## Methode regression lineaire (event. avec sigma-clipping)
        ##
    #    a1, b1, r_value, p_value, std_err = stats.linregress(ref_fint,obs_fint) # reglin 1 fois
        a1, b1, r_value, p_value, std_err = stats.linregress(ref_fint,obs_fselec) # reglin 1 fois
        ## Visualisation resultat
        plt.subplot(211)
        plt.title('Least square')
        plt.plot(xreg,a1*xreg+b1,'-k')
    #    plt.plot(ref_fint,obs_fint,'or')
        plt.plot(ref_fint,obs_fselec,'or')
        ##
        ## Methode RLM
        ##
        rlmx = sm.add_constant(ref_fint)
    #    rlmy = obs_fint
        rlmy = obs_fselec
        rlm_model = sm.RLM(rlmy,rlmx,M=sm.robust.norms.HuberT())
        rlm_results = rlm_model.fit()
        [b2, a2] = rlm_results.params
        ## Visualisation resulta
        plt.subplot(212)
        plt.title('RLM')
        plt.plot(xreg,a2*xreg+b2,'-k')
    #    plt.plot(ref_fint,obs_fint,'or')
        plt.plot(ref_fint,obs_fselec,'or')
        filename = cfgsetup.get('FullPaths', 'Event')\
                + cfgsetup.get('RelativePaths', 'Plots')\
                + obs_name + '_reglin.pdf'
        plt.savefig(filename)
        plt.close(2)

        ## Trace des CL alignees pour toutes les donnees
        obs_fluxali1 = (obs_flux-b1)/a1 # methode 1
        obs_fluxali2 = (obs_flux-b2)/a2 # methode 2
        ##
    #    fig1 = plt.figure(4)
    #    ax1 = plt.subplot(211)
    #    plt.title('Least square')
    #    plt.plot(ref_date,ref_flux,'or',obs_date,obs_fluxali1,'og')
    #    ax2 = plt.subplot(212)
    #    plt.title('RLM')
    #    plt.plot(ref_date,ref_flux,'or',obs_date,obs_fluxali2,'og')
    #    plt.savefig("LCalign.pdf")
    #    plt.close(4)

        cond = (photometryData['reject'] == 0) &\
               (photometryData['obs'] == obs_name)
        photometryData.loc[cond, 'mag_align'] = -2.5*np.log10(obs_fluxali2)

        cond2 = pd.isnull(photometryData['mag_align'])
        if np.any(np.array(cond2)):
            photometryData.loc[cond, 'mag_align'] = int(0)

        cond = (photometryData['reject'] == 0) &\
               (photometryData['obs'] == tel_ref)
        photometryData.loc[cond, 'mag_align'] = -2.5*np.log10(ref_flux)

        return 0

    # ----------------------------------------------------------------------
    #   Loading data full information.
    # ----------------------------------------------------------------------

    # Call function to align data - WRAP
    # ----------------------------------------------------------------------
    tel_list = np.array(obs_properties['key'])
    tel_col = np.array(['#' + a for a in obs_properties['colour']])
    tel_ref = np.where(tel_list == cfgsetup.get('Observatories', 'Reference'))[0][0]

    text = " Telescope de reference : " + obs_properties['name'][tel_ref]
    communicate(cfgsetup, 1, text)
    obs_list = np.delete(tel_list, tel_ref)
    tel_col = np.delete(tel_col, tel_ref)

    photometryData = pd.DataFrame({})
    for key in time_serie:
        photometryData[key]=time_serie[key]
    photometryData['mag_align'] = 0
    photometryData['reject'] = 0
    photometryData['color'] = '#000000'
    photometryData['obs_leg'] = '?'
    photometryData['alpha'] = 0.7

    # Trace interactif RLM
    for i in xrange(len(obs_list)):
        obs_name = obs_list[i]
        text = " %s" % obs_name
        communicate(cfgsetup, 1, text)
        align(photometryData, tel_list[tel_ref], obs_name, cfgsetup)

    # ----------------------------------------------------------------------
    #   Create an html webpage
    # ----------------------------------------------------------------------
    filename = cfgsetup.get('FullPaths', 'Event')\
            + cfgsetup.get('RelativePaths', 'Plots')\
            + cfgsetup.get('Controls', 'Archive') + '-datalign.html'

    bplt.output_file(filename)
    fig = np.array([])

    # ..................................................................
    #   Plot 0
    # ..................................................................

    # cond = photometryData['obs'] == 'moa-i'
    # photometryData = photometryData.loc[cond]

    observatories = np.unique(photometryData['obs'])
    for i in xrange(len(observatories)):
        cond = photometryData['obs'] == observatories[i]

        photometryData.loc[cond, 'color'] = '#' +\
                obs_properties['colour'][np.where(np.array(obs_properties['key']) == observatories[i])[0][0]]

        photometryData.loc[cond, 'obs_leg'] = \
                obs_properties['name'][np.where(np.array(obs_properties['key']) == observatories[i])[0][0]]

    cond = photometryData['reject'] == 0
    photometryData = photometryData[cond]

    source = ColumnDataSource(photometryData)

    hover0 = HoverTool(
            tooltips=[
                ("ID", "@id{int}"),
                ("Obs", "@obs"),
                ("Date", "@dates{1.11}"),
                ("Mag", "@mag_align{1.111}"),
                ("Err", "@err_magn{1.111}"),
                ("Seeing", "@seeing{1.11}"),
                ("Bkg", "@background{1.11}"),
            ]
        )

    # Plot the data
    tools = ["save", "pan", "box_zoom", "wheel_zoom", "reset", "tap", hover0]
    cond = (photometryData['reject'] == 0)
    # xmin = np.min(photometryData[cond]['dates'].values)
    # xmax = np.max(photometryData[cond]['dates'].values)
    xmin = float(options.split('-')[0].strip())
    xmax = float(options.split('-')[1].strip())
    ymin = np.min(photometryData[cond].mag_align.values)
    ymax = np.max(photometryData[cond].mag_align.values)
    fig = np.append(fig,\
            bplt.figure(toolbar_location="above", plot_width=1200, plot_height=600, x_range=(xmin, xmax), y_range=(ymax, ymin),\
            title=None, min_border=10, min_border_left=50, tools=tools))

    fig[0].circle('dates', 'mag_align', size=8, color='color', alpha='alpha', source=source)

    # Plot the errorbars
    for i in xrange(len(observatories)):
        cond = photometryData['obs'] == observatories[i]

        fig[0].circle(0, 0, size=8, color=photometryData[cond]['color'].values[0], alpha=0.8, legend=photometryData[cond]['obs_leg'].values[0])

        x = photometryData[cond]['dates']
        y = photometryData[cond]['mag_align']
        color = "#" + obs_properties['colour'][np.where(np.array(obs_properties['key']) == observatories[i])[0][0]]

        err_xs = []
        err_ys = []
        err_xs2 = []
        err_ys2 = []
        err_alpha = []
        for x, y, yerr, alpha_c in zip(x, y, photometryData[cond]['err_magn'], photometryData[cond]['alpha']):
            err_xs.append((x, x))
            err_xs2.append((x-5, x+5))
            err_ys.append((y - yerr, y + yerr))
            err_ys2.append((y - yerr, y - yerr))
            err_alpha.append(alpha_c)

        fig[0].multi_line(err_xs, err_ys, color=color, alpha=err_alpha)

    # Layout
    # ^^^^^^
    fig[0].xaxis.axis_label = 'HJD - 2,450,000'
    fig[0].yaxis.axis_label = 'Magnitude'
    fig[0].xaxis.axis_label_text_font = 'helvetica'
    fig[0].yaxis.axis_label_text_font = 'helvetica'
    fig[0].xaxis.axis_label_text_font_size = '10pt'
    fig[0].yaxis.axis_label_text_font_size = '10pt'

    fig[0].min_border_top = 10
    fig[0].min_border_bottom = 0
    fig[0].min_border_left = 0

    fig[0].xgrid.grid_line_color = None
    fig[0].ygrid.grid_line_color = None


    # ..................................................................
    #   Plot 1
    # ..................................................................

    # Plot the residus
    hover1 = HoverTool(
            tooltips=[
                ("ID", "@ID{int}"),
                ("Obs", "@obs"),
                ("Date", "@dates{1.11}"),
                ("Mag", "@mag_align{1.111}"),
                ("Err", "@err_magn{1.111}"),
                ("Seeing", "@seeing{1.11}"),
                ("Bkg", "@background{int}"),
            ]
        )

    tools = ["save", "pan", "box_zoom", "wheel_zoom", "reset", "tap", hover1]
    fig = np.append(fig,\
            bplt.figure(toolbar_location="above", plot_width=fig[0].plot_width, plot_height=400, x_range=fig[0].x_range,
            title=None, min_border=10, min_border_left=50, tools=tools))

    fig[1].circle('dates', 'seeing', size=8, color='color', alpha='alpha', source=source)


    # Layout
    # ^^^^^^
    fig[1].xaxis.axis_label = 'HJD - 2,450,000'
    fig[1].yaxis.axis_label = 'Seeing'
    fig[1].xaxis.axis_label_text_font = 'helvetica'
    fig[1].yaxis.axis_label_text_font = 'helvetica'
    fig[1].xaxis.axis_label_text_font_size = '10pt'
    fig[1].yaxis.axis_label_text_font_size = '10pt'

    fig[1].min_border_top = 10
    fig[1].min_border_bottom = 0
    fig[1].min_border_left = 0

    fig[1].xgrid.grid_line_color = None
    fig[1].ygrid.grid_line_color = None

    # ..................................................................
    #   Plot 2
    # ..................................................................

    # Plot the residus
    hover2 = HoverTool(
            tooltips=[
                ("ID", "@ID{int}"),
                ("Obs", "@obs"),
                ("Date", "@dates{1.11}"),
                ("Mag", "@mag_align{1.111}"),
                ("Err", "@err_magn{1.111}"),
                ("Seeing", "@seeing{1.11}"),
                ("Bkg", "@background{int}"),
            ]
        )

    tools = ["save", "pan", "box_zoom", "wheel_zoom", "reset", "tap", hover2]
    fig = np.append(fig,\
            bplt.figure(toolbar_location="above", plot_width=fig[0].plot_width, plot_height=400, x_range=fig[0].x_range,
            title=None, min_border=10, min_border_left=50, tools=tools))

    fig_curr = fig[2]

    fig_curr.circle('dates', 'background', size=8, color='color', alpha='alpha', source=source)

    # Layout
    # ^^^^^^
    fig_curr.xaxis.axis_label = 'HJD - 2,450,000'
    fig_curr.yaxis.axis_label = 'Background'
    fig_curr.xaxis.axis_label_text_font = 'helvetica'
    fig_curr.yaxis.axis_label_text_font = 'helvetica'
    fig_curr.xaxis.axis_label_text_font_size = '10pt'
    fig_curr.yaxis.axis_label_text_font_size = '10pt'

    fig_curr.min_border_top = 10
    fig_curr.min_border_bottom = 0
    fig_curr.min_border_left = 0

    fig_curr.xgrid.grid_line_color = None
    fig_curr.ygrid.grid_line_color = None

    # ..................................................................
    #   Table
    # ..................................................................

    columns = [
            TableColumn(field="id", title="Data ID", width=50),
            TableColumn(field="obs", title="Observatory", width=50),
            TableColumn(field="dates", title="Date", width=100),
            TableColumn(field="mag_align", title="Magnitude (Output code)", width=100),
            TableColumn(field="err_magn", title="Err_Magnitude", width=100),
            TableColumn(field="seeing", title="Seeing", width=100),
            TableColumn(field="background", title="Background", width=100),
        ]

    data_table = DataTable(source=source, columns=columns, width=1200, height=280)

    # Save
    # ^^^^
    final = blyt.column(fig[0], fig[1], fig[2], blyt.WidgetBox(data_table))
    bplt.save(final)

    # ------------------------------------------------------------------
    #   Modify the html page
    # ------------------------------------------------------------------
    filename = cfgsetup.get('FullPaths', 'Event')\
            + cfgsetup.get('RelativePaths', 'Plots')\
            + cfgsetup.get('Controls', 'Archive') + '-datalign.html'

    title = cfgsetup.get('EventDescription',
                             'Name') + ' - Datalign'

    file = open(filename, 'r')
    file_new = ''
    for line in file:
        # print line.strip()[:7]
        if line.strip()[:7] == '<title>':
            file_new = file_new \
                       + '        <style type="text/css">\n' \
                       + '        p.p1 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 43.0px; font: 36.0px "Lucida Grande"; color: #000000; -webkit-text-stroke: #000000}\n' \
                       + '        p.p2 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 21.0px; font: 18.0px "Lucida Grande"; color: #000000; -webkit-text-stroke: #000000}\n' \
                       + '        p.p3 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 15.0px; font: 12.0px "Lucida Grande"; color: #000000; -webkit-text-stroke: #000000}\n' \
                       + '        p.p4 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 17.0px; font: 14.0px "Lucida Grande"; color: #000000; -webkit-text-stroke: #000000; min-height: 17.0px}\n' \
                       + '        p.p5 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 14.0px; font: 12.0px Times; color: #000000; -webkit-text-stroke: #000000; min-height: 14.0px}\n' \
                       + '        p.p6 {margin: 0.0px 0.0px 12.0px 0.0px; line-height: 14.0px; font: 12.0px Times; color: #000000; -webkit-text-stroke: #000000; min-height: 14.0px}\n' \
                       + '        p.p7 {margin: 0.0px 0.0px 0.0px 0.0px; line-height: 17.0px; font: 14.0px "Lucida Grande"; color: #000000; -webkit-text-stroke: #000000}\n' \
                       + '        span.s1 {font-kerning: none}\n' \
                       + '        span.s10 {font: 14.0px "Lucida Grande"; color: #585858}\n' \
                       + '        hr {\n' \
                       + '        display: block;\n' \
                       + '        margin-top: 0.5em;\n' \
                       + '        margin-bottom: 0.5em;\n' \
                       + '        margin-left: auto;\n' \
                       + '        margin-right: auto;\n' \
                       + '        border-style: inset;\n' \
                       + '        border-width: 1px;\n' \
                       + '        }\n' \
                       + '        </style>\n' \
                       + '        <title>' + 'muLAn ' + cfgsetup.get('EventDescription', 'Name')[4:] + '- Datalign</title>\n' \
                       + '        <meta name="Author" content="Clement Ranc">\n'
        elif line.strip()[:6] == '<body>':
            file_new = file_new \
                       + '    <body>\n\n' \
                       + '<p class="p1"><span class="s1"><b>' + title + '</b></span></p>\n' \
                       + '<p class="p2"><span class="s1"><br>\n' \
                       + '</span></p>\n'
        elif line.strip()[:7] == '</body>':
            file_new = file_new \
                       + '        <BR>\n' \
                       + '        <hr>\n' \
                       + '        <BR>\n' \
                       + '        <footer>\n' \
                       + '        <p class="p7"><span class="s10">Modelling and page by muLAn (MicroLensing Analysis software).</span></p>\n' \
                       + '        <BR>\n' \
                       + '        <BR>\n' \
                       + '        <BR>\n' \
                       + '        <BR>\n' \
                       + '        <BR>\n' \
                       + '        </footer>\n' \
                       + '        </body>\n'
        else:
            file_new = file_new + line
    file.close()

    file = open(filename, 'w')
    file.write(file_new)
    file.close()

#
# ======================================================================
#   Main
# ======================================================================
if (__name__ == "__main__"):

    sys.exit("Option not supported yet.")








