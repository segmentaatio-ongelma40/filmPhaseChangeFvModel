import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

FONT_SIZE=30
LINE_WIDTH=4.0

# start-of-injection (CFD cases)
SOI = 0.5 # ms

# ------------------------------- FUNCTIONS --------------------------------------- #

def calculate_evap_rate(df):
    new_df = pd.DataFrame()

    new_df["dt"] = df["t"].diff().iloc[1:]
    new_df["dm"] = df["m"].diff().iloc[1:]
    new_df["mdot"] = -new_df["dm"]/(new_df["dt"]*0.001) # convert to mg/s

    new_df["t"] = df["t"].to_numpy()[:-1] + new_df["dt"].to_numpy()/2

    return new_df

def calculate_absolute_error(measured_val, correct_val):
    return measured_val - correct_val

def calculate_relative_error(measured_val, correct_val):
    return calculate_absolute_error(measured_val, correct_val) / correct_val

def plot_mass(data, ax, style='-', label='CFD', color='black', markevery=None):
    ax.plot(
            data["t"][1:], 
            data["m"][1:], 
            style, 
            label=label, 
            linewidth=LINE_WIDTH,
            markersize=LINE_WIDTH*5,
            markerfacecolor='none',
            markeredgewidth=LINE_WIDTH*0.66,
            color=color,
            markevery=markevery
            )
    ax.set_ylabel("film mass [mg]", fontsize=FONT_SIZE)
    ax.set_xlabel("$t$ [ms aSOI]", fontsize=FONT_SIZE)
    ax.tick_params(labelsize=FONT_SIZE)
    ax.grid(True)
    ax.set_xlim([-25, 625])
    ax.set_ylim([0, 0.45])

def plot_rate(data, ax, style='-', label='CFD', color='black', markevery=None):
    ax.plot(
            data["t"], 
            data["mdot"], 
            style, 
            label=label, 
            linewidth=LINE_WIDTH,
            markersize=LINE_WIDTH*5,
            markerfacecolor='none',
            markeredgewidth=LINE_WIDTH*0.66,
            color=color,
            markevery=markevery
            )
    ax.set_ylabel("film evaporation rate [mg/s]", fontsize=FONT_SIZE)
    ax.set_xlabel("$t$ [ms aSOI]", fontsize=FONT_SIZE)
    ax.tick_params(labelsize=FONT_SIZE)
    ax.grid(True)
    ax.set_xlim([-25, 625])
    ax.set_ylim([0.15, 0.75])

def plot_ref(case_path, axs):
    print("Printing reference data")

    refdata_path = case_path + "/experimental/jungst-et-al_2021_fig8_lif.csv"
    refdata = pd.read_table(
            refdata_path,
            comment='#',
            names=["t", "m"],
            sep=','
            )

    plot_mass(
            refdata, 
            axs[0], 
            style='k:+', 
            label='Jüngst et. al 2021',
            )

    # clean up spike in the beginning
    refdata_short = refdata.drop(refdata[refdata.t < 150].index)
    plot_rate(
            calculate_evap_rate(refdata_short), 
            axs[1], 
            style='k:+', 
            label='Jüngst et. al 2021'
            )

def plot_cfd(case_path, axs, style='-', label='CFD', color='black', markevery=None):

    print("Printing case ", case_path)

    data_path = case_path + "/postProcessing/film/volumeIntegral/0/volFieldValue.dat"
    data = pd.read_table(
            data_path, 
            comment='#', 
            names=["t", "m"], 
            sep='\s+'
            )

    # correct data...
    data["t"] = data["t"] * 10**(3) # convert time to ms
    data["m"] = data["m"] * 10**(6) # convert mass to mg

    data["t"] = data["t"] - SOI

    plot_mass(
            data, 
            axs[0], 
            style=style, 
            label=label,
            color=color,
            markevery=markevery
            )

    # clean up spike in the beginning
    data_short = data.drop(data[data.t < 50].index)
    plot_rate(
            calculate_evap_rate(data_short), 
            axs[1], 
            style=style, 
            label=label,
            color=color,
            markevery=markevery
            )


# - - - - - - - - - - - - - - - END FUNCTIONS - - - - - - - - - - - - - - - - - - - #

# ------------------------------- PLOT MASS --------------------------------------- #

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["DejaVu Serif"]

cases = [
          # path,  label,   style
          ("..", "1.830 mm", '-')
        ]

idx = 0
cmap=plt.get_cmap('tab10')
idxs=np.linspace(0.1, 0.9, num=len(cases))
colors = cmap(idxs)

fig, axs = plt.subplots(1,2)

for c in cases:
    color = colors[idx]
    plot_cfd(c[0], axs, label=c[1], style=c[2], color=color, markevery=0.1)
    idx = idx + 1

plot_ref("..", axs)

fig.set_size_inches(21.0, 10.5)
fig.tight_layout()
plt.legend(loc="upper right", fontsize=FONT_SIZE)
plt.savefig('masses.png')
plt.show()

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ##


