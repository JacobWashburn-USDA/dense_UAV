#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os, sys
import pandas as pd
import numpy as np
import plotly.express as px
import io

from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd


# In[2]:


#bring in data
data_dir = "../data/"
sng_long_pred = pd.read_csv("../data/Prediction ability.csv")
sng_long_pred_noPHT = pd.read_csv("../data/Prediction ability no hieght.csv")

temp_VIs = pd.read_csv(data_dir+"Temporal_VIs.csv", index_col=[0])


# In[3]:


def process_labels(sng_long_pred):
    #add M number as int for sorting
    sng_long_pred["M_num"] = sng_long_pred["Model"].str.split("M",expand=True)[1].astype(int)
    #add genomic verses phenomic
    sng_long_pred["Type"] = ""
    sng_long_pred.loc[sng_long_pred["Model"]=="M1", "Type"] = "Genomic"
    sng_long_pred.loc[sng_long_pred["Model"]!="M1", "Type"] = "Phenomic"
    return sng_long_pred


# In[4]:


def add_DAP(temp_VIs, sng_long_pred):
    #get and add DAP to labels
    temp_VIs=temp_VIs.copy()
    temp_VIs.columns = pd.MultiIndex.from_tuples([x.split(" ") for x in temp_VIs.columns.tolist()])
    daps = [0] + temp_VIs["Red"].columns.tolist()
    daps = pd.DataFrame(daps, columns=["DAP"])
    daps["Model"]=["M"+str(x) for x in range(1, len(daps)+1)]

    labels = []
    for ind in daps.index:
        label = daps.loc[ind, "Model"] + ": "
        if ind == 0:
            label = "M1: Genomic"
        elif ind == 1:
            label += daps.loc[ind,"DAP"] + " DAP"
        else:
            label += daps.loc[ind-1, "Model"] + " + " + daps.loc[ind,"DAP"] + " DAP"
        labels.append(label)
        #print(label)
    daps["Label"] = labels
    #long_pred = long_pred.merge(daps[["Model","Label"]], on="Model")
    sng_long_pred = sng_long_pred.merge(daps[["Model","Label"]], on="Model")
    return sng_long_pred


# In[5]:


def anova_tukey(df, group_column="Model", alpha = 0.05):
    print(len(df))
    models = []
    for model in df[group_column].unique():
        models.append(df.loc[df[group_column]==model, "cor"].tolist())
    anova_res = f_oneway(*models)
    print(anova_res)
    
    if anova_res[1] < alpha:
        print("Preforming Tukey's HSD test...")
    
        tukey = pairwise_tukeyhsd(endog=df['cor'],
                                  groups=df[group_column],
                                  alpha=alpha)
        return tukey
    else: return None 


# In[6]:


def run_tukey(sng_long_pred):
    tukey_res = anova_tukey(sng_long_pred[sng_long_pred["CV"]=="CV1"]) #, group_column="M_num")
    #get single confidence interval and mean data
    tukey_res._simultaneous_ci()
    means = tukey_res._multicomp.groupstats.groupmean
    #minrange = [means[i] - tukey_res.halfwidths[i] for i in range(len(means))]
    #maxrange = [means[i] + tukey_res.halfwidths[i] for i in range(len(means))]
    sim_CIs = pd.DataFrame([tukey_res.groupsunique, means, tukey_res.halfwidths], index=["Model","Mean","Halfwidth"]).T
    sim_CIs["M_num"] = sim_CIs["Model"].str[1:].astype(int)
    sim_CIs = sim_CIs.sort_values("M_num").reset_index(drop=True)
    sim_CIs = sim_CIs.merge(sng_long_pred[["Model","Label"]].drop_duplicates(), on="Model", how="left")
    return sim_CIs


# In[7]:


sng_long_pred = process_labels(sng_long_pred)
sng_long_pred = add_DAP(temp_VIs, sng_long_pred)
sim_CIs = run_tukey(sng_long_pred) #run anova


# In[8]:


#FIGURE 5
for_fig = sng_long_pred[sng_long_pred["CV"]=="CV1"].copy()
for_fig = for_fig[for_fig["Model"].isin(["M1","M44","M15"])]
for_fig["Model"] = for_fig["Model"].replace({"M1":"Genomic<br>(M1)",
                                             "M44":"Phenomic:<br>All Season<br>(M2)", 
                                             "M15":"Phenomic:<br>Pre-Flowering<br>(M3)"})
for_fig = for_fig.sort_values("Model")
fig = px.box(for_fig, x="Model", y="cor", color="Model")
fig.update_layout(yaxis_title="Prediction Ability (pearson r)", xaxis_title="")
fig.update_layout(showlegend=False)

width = 500
height = 500
fig.update_layout(autosize=False, width=width,height=height)
#fig.write_html("../Figures/Fig5.html")
#fig.write_image("../Figures/Fig5.svg")
#fig.write_image("../Figures/Fig5.png")
fig.write_image("../Figures/Fig5.pdf")
fig.show()

#ALT FIGURE 5
#mean = for_fig.pivot_table(index=["Model"]).reset_index()
#std = for_fig.pivot_table(index=["Model"], aggfunc="std").reset_index()
#fig = px.bar(mean, x="Model", y="cor", color="Model",
#             text=mean["cor"].round(2).to_list(),
#             error_y= std["cor"].tolist())
#fig.show()


# In[ ]:





# In[9]:


def tukey_simultanious_fig(sim_CIs, save_html_path="", save_svg_path="", width=900, height=750):
    fig = px.scatter(sim_CIs, x="Label", y="Mean", error_y="Halfwidth")
    fig.update_layout(yaxis_title="Prediction Ability (pearson r)", xaxis_title="")
    fig.update_layout(margin=dict(l=1, r=1, t=1, b=1))

    #add hline 1
    graph_line = sim_CIs[sim_CIs["Model"]=="M44"].copy().iloc[0]
    graph_line = graph_line["Mean"]-graph_line["Halfwidth"]
    fig.add_hline(y=graph_line, line_dash="dot", line_width=1)

    #add hline 2
    graph_line = sim_CIs[sim_CIs["Model"]=="M2"].copy().iloc[0]
    graph_line = graph_line["Mean"]+graph_line["Halfwidth"]
    fig.add_hline(y=graph_line, line_dash="dot", line_width=1)

    #add flowering time box
    fig.add_vrect(x0=15.5, x1=18.5)
    fig.update_xaxes(range=(-1,44))

    fig.update_layout(autosize=False, width=width,height=height)
    
    if save_html_path != "":
        fig.write_html(save_html_path)
    
    if save_svg_path != "":
        fig.write_image(save_svg_path)#, scale=scale)

    fig.show()


# In[10]:


#FIGURE 6
tukey_simultanious_fig(sim_CIs,
                       #save_html_path="../Figures/Fig6.html",
                       #save_svg_path="../Figures/Fig6_noMod.svg"
                       #save_svg_path="../Figures/Fig6_noMod.pdf"
                      )


# In[11]:


sng_long_pred_noPHT = process_labels(sng_long_pred_noPHT)
sng_long_pred_noPHT = add_DAP(temp_VIs, sng_long_pred_noPHT)
sim_CIs_noPHT = run_tukey(sng_long_pred_noPHT) #run anova


# In[12]:


#SUP FIGURE 9
tukey_simultanious_fig(sim_CIs_noPHT,
                       #save_html_path="../Figures/SupFig8.html",
                       #save_svg_path="../Figures/SupFig8.svg"
                      )


# In[ ]:




