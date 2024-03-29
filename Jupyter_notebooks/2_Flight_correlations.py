#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os, sys
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import spearmanr
from scipy.stats import linregress
from sklearn.metrics import r2_score
import plotly.express as px
import plotly.graph_objects as go
import xarray as xr


# In[2]:


data_dir = "../data/"
#bring in manual phenotype BLUES
manual_phenos = pd.read_csv(data_dir+"Pheno_blues.csv", index_col=[0])
manual_phenos.index.name="Pedigree"
manual_phenos.index = manual_phenos.index.str.replace("/"," X ")
#add calculated manual phenos ASI and EHT/PHT
manual_phenos["ASI"] = manual_phenos["DTS"] - manual_phenos["DTA"]
manual_phenos["EHT/PHT"] = manual_phenos["EHT"]/manual_phenos["PHT"]

#bring in temporal vegitation index data
temp_VIs = pd.read_csv(data_dir+"Temporal_VIs.csv", index_col=[0]) #Phenomic data Washburn.csv
#create 2-level column index with first level being VI and second level DAP
temp_VIs.columns = pd.MultiIndex.from_tuples([x.split(" ") for x in temp_VIs.columns.tolist()])

#bring in temporal plant hieght data (weilbull fit)
temp_ht = pd.read_csv(data_dir+"Tpht.csv", index_col=[0])
temp_ht.columns = pd.MultiIndex.from_tuples([x.split(".") for x in temp_ht.columns.tolist()])
#fix pedegree name
temp_ht.index = temp_ht.index.str.replace("W10004_0019X PHZ51", "W10004_0019 X PHZ51")

#combine hieght data with the VIs
all_phenos = pd.concat([temp_VIs, temp_ht], axis=1)

#bring in gdds
gdd = pd.read_csv(data_dir+"GDDs_DAPs.csv", index_col=[0], parse_dates=[0])

#flowering windows
flower_DAP = (51.19470817271569, 57.2200534210139)
flower_GDD = (750.1523934725153, 838.2223942479286)


# In[3]:


#create all phenots with GDD instead of DAP
#create DAP to GDD dictionary
dap_to_gdd_dict = gdd.copy()
#dap_to_gdd_dict[dap_to_gdd_dict["DAP"].duplicated()]
dap_to_gdd_dict.index = dap_to_gdd_dict["DAP"]
dap_to_gdd_dict = dap_to_gdd_dict["Cum_GDD [C]"].to_dict()

tmp_index = all_phenos.columns.to_frame().copy()
tmp_index[1] = tmp_index[1].astype(int)
tmp_index[1] = tmp_index[1].replace(dap_to_gdd_dict)
#tmp_index[1] = tmp_index[1].round(0).astype(int).astype(str)
#tmp_index[1].max(), tmp_index[1].min()
tmp_index = pd.MultiIndex.from_frame(tmp_index)

all_phenos_GDD = all_phenos.copy()
all_phenos_GDD.columns = tmp_index


# In[4]:


#custom helper functions
#custom rounding function with any base
def custom_round(pdSeries, base=1):
    return (pdSeries / base).round().astype(type(base)) * base

def get_overall_stats(obs_values, pred_values):
    stats_dict = {}
    #calculates statistics on complete dataset given with no spliting by environment or averageing
    #print(obs_values.name, pred_values.name)
    stats_dict["pearson_r"] = pearsonr(np.array(obs_values), np.array(pred_values))[0]
    #stats_dict["pearson_r_p"] = pearsonr(np.array(obs_values), np.array(pred_values))[1]
    return stats_dict

def calc_metrics(df):
    x = df.dropna().iloc[:,0] #.tolist()
    y = df.dropna().iloc[:,1] #.tolist()
    stats_dict = get_overall_stats(x, y)
    #r, pVal = pearsonr(x,y)
    return stats_dict

def calc_all_days(sngl_VI):
    sngl_VI = sngl_VI.copy()
    statsDF = []
    done = []
    for firstDay in sngl_VI.columns.sort_values(ascending=True).tolist():
        for scndDay in sngl_VI.columns.sort_values(ascending=True).tolist():
            distance = abs(scndDay-firstDay)
            if set([firstDay, scndDay]) in done:
                continue
            else:
                done.append(set([firstDay, scndDay]))
                mtrcs_dict = calc_metrics(sngl_VI[[firstDay, scndDay]].copy())
                statsDF.append([firstDay, scndDay, distance] + list(mtrcs_dict.values()))
    statsDF = pd.DataFrame(statsDF, columns=["First Date [DAP]", "Second Date [DAP]", "Distance [DAP]"]+list(mtrcs_dict.keys()))
    return statsDF


# In[5]:


#UNCOMENT IF NEED TO REGENERATE
#calculate correlations and days apart for all VIs
'''
all_VIs_stats = []
for vi in all_phenos.columns.to_frame()[0].unique():
    print(vi)
    sngl_VI = all_phenos[vi].copy()
    sngl_VI.columns = sngl_VI.columns.astype(int)
    statsDF = calc_all_days(sngl_VI)
    statsDF.index = pd.MultiIndex.from_frame(statsDF[["First Date [DAP]", "Second Date [DAP]", "Distance [DAP]"]])
    statsDF = statsDF.drop(columns=["First Date [DAP]", "Second Date [DAP]", "Distance [DAP]"])
    all_VIs_stats.append(statsDF.rename(columns={"pearson_r":vi}).copy())
all_VIs_stats = pd.concat(all_VIs_stats, axis=1)
#save for future easy use.
all_VIs_stats.to_csv("../data/all_VIs_stats.csv")
'''


# In[6]:


#Figure ploting functions

def correlation_heatmap(all_phenos, vi, x_y_label = "Days After Planting", value_col = "Pearson r", title="",
                        save_html_path="", save_svg_path="", fig_width=800, fig_height=400):
    for_corr = all_phenos[vi].copy()
    for_corr = for_corr.corr()
    
    fig = px.imshow(for_corr, origin="lower", #height=fig_height, width=fig_width,
                    labels=dict(x=x_y_label, y=x_y_label, color=value_col))
    fig.update_layout(title_text=title, title_x=0.5)
    fig.update_layout(margin=dict(l=0, r=1, t=30, b=1))
    
    if save_html_path != "":
        fig.write_html(save_html_path)
    
    if save_svg_path != "":
        fig.update_layout(autosize=False, width=fig_width, height=fig_height)
        fig.write_image(save_svg_path)#, scale=scale)
    fig.show()
    
    return for_corr

def box_plot_w_trendline(all_VIs_stacked, vi, x_col="Distance [DAP]", value_col = "Pearson r", title="",
                         x_title = "Days Between Flights", distance_cutoff = 120, save_html_path="", round_base=0,
                         save_svg_path="", fig_width=800, fig_height=400):
    for_box_plot = all_VIs_stacked[all_VIs_stacked["Vegetation Index"]==vi].copy()
    for_box_plot = for_box_plot[for_box_plot[x_col] < distance_cutoff]
    
    #round x_col values to round_base
    if round_base != 0:
        for_box_plot[x_col] = custom_round(for_box_plot[x_col], base=round_base)
    
    fig = px.box(for_box_plot, x = x_col, y=value_col)

    #draw trenlind over the box plot
    trendline = for_box_plot.groupby([x_col])[value_col].mean()

    x1,y1 = list(trendline.index), trendline.values
    #stds = for_box_plot.groupby(["Distance [DAP]"]).value.std().values
    #n_vals = for_box_plot.groupby(["Distance [DAP]"])["Distance [DAP]"].count().values

    fig.add_trace(go.Scatter(x=x1, y=y1, line_shape="spline", line_color="blue", name="trendline"))
    fig.update_layout(title_text=title, title_x=0.5)
    fig.update_layout(xaxis_title=x_title)
    fig.update_layout(margin=dict(l=0, r=1, t=30, b=1))
    fig.update_layout(showlegend=False)
    
    if save_html_path != "":
        fig.write_html(save_html_path)
    
    if save_svg_path != "":
        fig.update_layout(autosize=False, width=fig_width, height=fig_height)
        fig.write_image(save_svg_path)#, scale=scale)
    fig.show()
    
def round_2nd_index_values(all_phenos_GDD, base=1):
    #round second index values
    return_phenos_GDD = all_phenos_GDD.copy()
    tmp_index = return_phenos_GDD.columns.to_frame().copy()
    tmp_index[1] = custom_round(tmp_index[1], base).astype(int).astype(str)
    tmp_index = pd.MultiIndex.from_frame(tmp_index)
    return_phenos_GDD.columns = tmp_index
    return return_phenos_GDD


# In[7]:


#read VI stats from file
all_VIs_stats = pd.read_csv("../data/all_VIs_stats.csv", index_col=[0,1,2])
#stack the data for easier processing
all_VIs_stacked = all_VIs_stats.stack().reset_index().copy()
all_VIs_stacked = all_VIs_stacked.rename(columns={"level_3":"Vegetation Index", 0:"Pearson r"})

#add in GDD distance
all_VIs_stacked = all_VIs_stacked.merge(gdd[["DAP","Cum_GDD [C]"]], left_on="First Date [DAP]", right_on=["DAP"], how="left")
all_VIs_stacked = all_VIs_stacked.drop(columns="DAP").rename(columns={"Cum_GDD [C]":"1st Cum_GDD [C]"})
all_VIs_stacked = all_VIs_stacked.merge(gdd[["DAP","Cum_GDD [C]"]], left_on="Second Date [DAP]", right_on=["DAP"], how="left")
all_VIs_stacked = all_VIs_stacked.drop(columns="DAP").rename(columns={"Cum_GDD [C]":"2nd Cum_GDD [C]"})
all_VIs_stacked["Distance [GDD [C]]"] = (all_VIs_stacked["2nd Cum_GDD [C]"] - all_VIs_stacked["1st Cum_GDD [C]"]).abs()


# In[8]:


#FIGURE 1 A
vi = "RmB"
_ = correlation_heatmap(all_phenos, vi, 
                        x_y_label= "Days After Planting",
                        save_html_path="../Figures/Fig1A.html",
                        save_svg_path="../Figures/Fig1A.svg",
                        fig_width=450
                       )


# In[9]:


#FIGURE 1 B
vi = "RmB"
box_plot_w_trendline(all_VIs_stacked, vi, value_col="Pearson r", distance_cutoff=110,
                     save_html_path="../Figures/Fig1B.html",
                     save_svg_path="../Figures/Fig1B.svg",
                     fig_width=450
                    )


# In[ ]:





# In[10]:


#SUP FIGURE 2 A
vi = "RmB"
_ = correlation_heatmap(round_2nd_index_values(all_phenos_GDD, base=10),
                        vi,
                        x_y_label= "Cumulative Growing Degree Days [C]",
                        save_html_path="../Figures/SupFig2A.html",
                        save_svg_path="../Figures/SupFig2A.svg",
                        fig_width=450
                       )


# In[11]:


#SUP FIGURE 2 B
vi = "RmB"
box_plot_w_trendline(all_VIs_stacked, vi, x_col= "Distance [GDD [C]]", value_col="Pearson r",
                     x_title = "Growthing Degree Days [C] Between Flights", distance_cutoff=2000, round_base=15,
                     save_html_path="../Figures/SupFig2B.html",
                     save_svg_path="../Figures/SupFig2B.svg",
                     fig_width=450
                    )


# In[ ]:





# In[12]:


#SUP FIGURE 3
#make correlations for all VIs
drop_VIs = ["Height"]
for_corr = {}
for vi in all_phenos.columns.to_frame()[0].unique():
    if vi in drop_VIs: continue
    df = all_phenos[vi].copy() 
    df = df.corr()
    df.index = df.index.astype(int)
    df.columns = df.columns.astype(int)
    df = df.copy().stack()#.reset_index()
    df = df.to_xarray()
    for_corr[vi] = df
    #for_corr.append(df)
for_corr_plot = xr.concat([df for df in for_corr.values()], dim="concat_dim")#, dim=list(for_corr.keys()))#, name= "dataset")
fig = px.imshow(for_corr_plot, facet_col="concat_dim", origin="lower", facet_col_wrap=4, facet_row_spacing=0.02,
                labels=dict(x="DAP", y="DAP", color="Pearson r"))
fig.update_layout(margin=dict(l=0, r=1, t=30, b=1))
width = 800
height = 1000
fig.update_layout(autosize=False, width=width,height=height)
#fix facet labels
fig.for_each_annotation(lambda a: a.update(text=list(for_corr.keys())[int(a.text.split("=")[-1])]))
fig.write_html("../Figures/SupFig3.html")
fig.write_image("../Figures/SupFig3.svg")
fig.show()


# In[ ]:





# In[13]:


#SUP FIGURE 4
drop_VIs = ["Height"]
avg_by_dist = all_VIs_stats[[x for x in all_VIs_stats.columns if x not in drop_VIs]].copy()
avg_by_dist = avg_by_dist.reset_index().pivot_table(index=["Distance [DAP]"], values = [x for x in avg_by_dist.columns], aggfunc=np.mean).copy()
avg_by_dist = avg_by_dist.stack().reset_index().rename(columns={"level_1":"Vegetation Index", 0:"mean"})

for_box_facet_plot = all_VIs_stacked.copy()
for_box_facet_plot = for_box_facet_plot[for_box_facet_plot["Vegetation Index"].isin(drop_VIs)==False]
vis = for_box_facet_plot["Vegetation Index"].unique().tolist()
print("number of VIs:", len(vis))
cols=4
rows = int(len(vis)/4)
print("Rows: ", rows, "Cols: ", cols)

fig = px.line(avg_by_dist[avg_by_dist["Vegetation Index"].isin(drop_VIs)==False], x = "Distance [DAP]", y="mean",
              facet_col="Vegetation Index", facet_col_wrap=4, 
              facet_row_spacing=0.02)
fig.update_layout(margin=dict(l=0, r=1, t=30, b=1))
width = 800
height = 1000
fig.update_layout(autosize=False, width=width,height=height)
fig.write_html("../Figures/SupFig4.html")
fig.write_image("../Figures/SupFig4.svg")
fig.show()


# In[ ]:





# In[14]:


#find correlations between manual phenotypes and UAV phenotypes
def get_stats_single_vi_man_pheno(all_phenos, manual_phenos, man_pheno, vi):
    sngl_man_pheno = manual_phenos[man_pheno].copy()
    sngl_VI = all_phenos[vi].copy()
    man_uas_stats = []
    for day in sngl_VI.columns:
        #print(day)
        tmp_day = pd.concat([sngl_VI[day].copy(), sngl_man_pheno.copy()], axis=1)
        mtrcs_dict = calc_metrics(tmp_day)
        man_uas_stats.append([vi, day, man_pheno] + list(mtrcs_dict.values()))
    man_uas_stats = pd.DataFrame(man_uas_stats, columns=["Vegetation Index", "DAP", "Manual Phenotype"] + list(mtrcs_dict.keys()))
    return man_uas_stats
    
#man_uas_stats = get_stats_single_vi_man_pheno(all_phenos, manual_phenos, man_pheno = "Yield", vi = "Red")

man_uas_stats = []
for man_pheno in manual_phenos.columns:
    for vi in all_phenos.columns.to_frame()[0].unique():
        print(man_pheno, vi)
        man_uas_stats.append(get_stats_single_vi_man_pheno(all_phenos, manual_phenos, man_pheno, vi))
man_uas_stats = pd.concat(man_uas_stats).reset_index(drop=True)

#add GDD
print(len(man_uas_stats))
man_uas_stats["DAP"] = man_uas_stats["DAP"].astype(int)
man_uas_stats = man_uas_stats.merge(gdd[["DAP","Cum_GDD [C]"]].reset_index(drop=True), on="DAP", how="left")
#remove ASI and EHT/PHT
man_uas_stats_rm = man_uas_stats[man_uas_stats["Manual Phenotype"].isin(["ASI","EHT/PHT"])==False].copy()
print(len(man_uas_stats))


# In[15]:


def line_plot_manual_vi_cor(man_uas_stats, man_phenos, vi, flower_DAP=flower_DAP, x="DAP", save_html_path="", save_svg_path="",
                            xaxis_title="Days After Planting",
                            fig_width=800, fig_height=400):
    if len(man_phenos) == 0:
        man_phenos = man_uas_stats["Manual Phenotype"].unique().tolist()
    tmp_for_line_plot = man_uas_stats[(man_uas_stats["Manual Phenotype"].isin(man_phenos)) & 
                                      (man_uas_stats["Vegetation Index"]==vi)].copy()
    tmp_for_line_plot[x] = tmp_for_line_plot[x].astype(int)
    tmp_for_line_plot["pearson_r"].max()
    fig = px.line(tmp_for_line_plot, x=x, y="pearson_r", color="Manual Phenotype",
                  #line_dash=x,
                  #line_dash_sequence=["dashdot", "dash"] #["solid","dash","dashdot","dot","longdash"]
                 )    
    fig["data"][1]['line']['dash'] = "dash"
    fig["data"][2]['line']['dash'] = "dashdot"
    fig["data"][3]['line']['dash'] = "dot"
    fig["data"][4]['line']['dash'] = "longdashdot"
    #fig.update_layout(title_text=vi, title_x=0.5)
    fig.update_layout(xaxis_title=xaxis_title)
    fig.update_layout(margin=dict(l=0, r=1, t=30, b=1))
    fig.update_layout(yaxis_title="Pearson r")
    min_range = -0.12
    max_range= 0.55
    fig.update_yaxes(range=[min_range, max_range])
    fig.add_shape(type="rect", x0=flower_DAP[0], y0=min_range, x1=flower_DAP[1], y1=max_range,
                  line=dict(
                      color="black",
                      width=3,
                  ), fillcolor="black", opacity=0.25)
    fig.update_layout(legend_title_text='Phenotype')
    fig.update_layout(legend=dict(yanchor="top",
                                  y=0.99,
                                  xanchor="right",
                                  x=0.99))
    
    if save_html_path != "":
        fig.write_html(save_html_path)
    if save_svg_path != "":
        fig.update_layout(autosize=False, width=fig_width, height=fig_height)
        fig.write_image(save_svg_path)#, scale=scale)
    fig.show()
    return fig


# In[16]:


#FIGURE 2
vi="GmR"
fig = line_plot_manual_vi_cor(man_uas_stats, man_phenos=[], vi=vi,
                              save_html_path="../Figures/Fig2.html",
                              save_svg_path="../Figures/Fig2.svg",
                             )


# In[17]:


#SUP FIGURE 5
vi="GmR"
fig = line_plot_manual_vi_cor(man_uas_stats, man_phenos=[], vi=vi,
                              flower_DAP=flower_GDD,
                              x="Cum_GDD [C]",
                              xaxis_title="Cumulative Growing Degree Days [C]",
                              save_html_path="../Figures/SupFig5.html",
                              save_svg_path="../Figures/SupFig5.svg",
                             )


# In[18]:


#SUP FIGURE 6
#create plots for all VIs
drop_VIs = ["Height"]
tmp_for_plot = man_uas_stats.copy()
tmp_for_plot["DAP"] = tmp_for_plot["DAP"].astype(int)
tmp_for_plot = tmp_for_plot.rename(columns={"Vegetation Index":"VI"})
fig = px.line(tmp_for_plot[tmp_for_plot["VI"].isin(drop_VIs)==False], x = "DAP", y="pearson_r", color="Manual Phenotype",
              facet_col="VI", facet_col_wrap=4, 
              facet_row_spacing=0.02)
fig.update_layout(margin=dict(l=0, r=1, t=30, b=1))
width = 900
height = 1100
fig.update_layout(autosize=False, width=width,height=height)
fig.update_layout(legend=dict(yanchor="bottom",
                              y=-0.08,
                              xanchor="center",
                              x=0.5,
                              orientation="h",
                             ))

fig.write_html("../Figures/SupFig6.html")
fig.write_image("../Figures/SupFig6.svg")
fig.show()


# In[ ]:





# In[ ]:




