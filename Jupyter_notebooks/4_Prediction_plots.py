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
#sng_long_pred = pd.read_csv("../data/Prediction ability.csv")
sng_long_pred = pd.read_csv("../data/Prediction ability no hieght.csv")

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
    daps["Label"] = labels
    #long_pred = long_pred.merge(daps[["Model","Label"]], on="Model")
    sng_long_pred = sng_long_pred.merge(daps[["Model","Label"]], on="Model")

    #fix labels
    tmp = sng_long_pred.loc[sng_long_pred["Label"]!="M1: Genomic", "Label"].copy()
    sng_long_pred.loc[sng_long_pred["Label"]!="M1: Genomic", "Label"] = tmp.str.replace(" M\d+ \+", "", regex=True)#.unique()
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
sim_CIs = run_tukey(sng_long_pred)


# In[9]:


#FIGURE 6
fig = px.scatter(sim_CIs, x="Label", y="Mean", error_y="Halfwidth")
fig.update_layout(yaxis_title="Prediction Ability (pearson r)", xaxis_title="")
fig.update_layout(margin=dict(l=1, r=1, t=1, b=1))

#fig2 = px.line(sim_CIs, x="Label", y="Mean")
#fig2.update_traces(line_color='black')#, line_width=5)
#fig.add_traces(list(fig2.select_traces()))


width = 900
height = 750
fig.update_layout(autosize=False, width=width,height=height)
fig.write_html("../Figures/Fig6.html")
fig.write_image("../Figures/Fig6.svg")

fig.show()
#tukey_res.plot_simultaneous()


# In[ ]:





# In[16]:


#FIGURE 5
for_fig = sng_long_pred[sng_long_pred["CV"]=="CV1"].copy()
for_fig = for_fig[for_fig["Model"].isin(["M1","M44","M17"])]
for_fig.pivot_table(index=["Model"])


# In[15]:


fig = px.box(for_fig, x="Model", y="cor", color="Model")
fig.show()


# In[11]:


cv_cors


# In[ ]:





# In[ ]:





# In[ ]:





# In[2]:



long_pred = pd.read_csv(data_dir+"Longitudinal Prediction ability no hieght.csv", index_col=[0])

#add M number as int for sorting
long_pred["M_num"] = long_pred["Model"].str.split("M",expand=True)[1].astype(int)
#add genomic verses phenomic
long_pred["Type"] = ""
long_pred.loc[long_pred["Model"]=="M1", "Type"] = "Genomic"
long_pred.loc[long_pred["Model"]!="M1", "Type"] = "Phenomic"

sng_long_pred = pd.read_csv(data_dir+"single_flight_Prediction ability.csv", index_col=[0])

#add M number as int for sorting
sng_long_pred["M_num"] = sng_long_pred["Model"].str.split("M",expand=True)[1].astype(int)
#add genomic verses phenomic
sng_long_pred["Type"] = ""
sng_long_pred.loc[sng_long_pred["Model"]=="M1", "Type"] = "Genomic"
sng_long_pred.loc[sng_long_pred["Model"]!="M1", "Type"] = "Phenomic"


# In[ ]:




