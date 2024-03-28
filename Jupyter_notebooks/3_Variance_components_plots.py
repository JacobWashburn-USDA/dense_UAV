#!/usr/bin/env python
# coding: utf-8

# In[1]:


import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from plotly.subplots import make_subplots


# In[2]:


#bring in data
data_dir = "../data/"

#bring in temporal plant hieght data (weilbull fit) from Variance Components r script
temp_ht = pd.read_csv(data_dir+"Tpht_subset.csv", index_col=[0])

#bring in gdds
gdd = pd.read_csv("../data/GDDs_DAPs.csv", index_col=[0], parse_dates=[0])

#add in GDDS
temp_ht = temp_ht.merge(gdd[["DAP","Cum_GDD [C]"]], on="DAP")

#flowering windows
flower_DAP = (51.19470817271569, 57.2200534210139)
flower_GDD = (750.1523934725153, 838.2223942479286)

man_VCs = pd.read_csv("../data/Var_comps_manual_phenos_Tpht.csv", index_col=[0])
#fix naming
man_VCs["grp"] = man_VCs["grp"].str.replace("Pedigree","Hybrid")
man_VCs.rename(columns={"grp":"Variance Components", "Percent":"Explaned Percent Variation (%)",
                        "Heritability":"Repetability"}, inplace=True)


# In[ ]:





# In[3]:


#plot heights by date and value
def plot_heights(df_stacked, flower_DAP=flower_DAP, x_col="DAP", value_col="Height", color="Tester",
                 testers_to_keep=['PHK76','PHP02','PHZ51'], xaxis_title = "Days After Planting (DAP)",
                 save_html_path="", save_svg_path="", fig_width=2000, fig_height=800):
    #filter testers
    for_fig = df_stacked[df_stacked["Tester"].isin(testers_to_keep)].copy()
    print("Number of Hybrids: ", len(df_stacked["Pedigree"].unique()), len(for_fig["Pedigree"].unique()))
    fig = px.box(for_fig, x = x_col, y=value_col, color=color)
    
    fig.update_layout(margin=dict(l=0, r=1, t=30, b=1))
    fig.update_layout(xaxis_title=xaxis_title)
    fig.update_layout(yaxis_title="Plant Height")
    
    #add flowering time box
    min_range = -0.4
    max_range= 2.6
    fig.update_yaxes(range=[min_range, max_range])
    fig.add_shape(type="rect", x0=flower_DAP[0], y0=min_range, x1=flower_DAP[1], y1=max_range,
                  line=dict(
                      color="black",
                      width=3,
                  ), fillcolor="black", opacity=0.25)
    fig.update_layout(legend=dict(yanchor="bottom",
                                  y=0.05,
                                  xanchor="right",
                                  x=0.99))
    fig.show()
    fig.update_layout(autosize=False, width=fig_width, height=fig_height)
    
    if save_html_path != "":
        fig.write_html(save_html_path)
    
    if save_svg_path != "":
        fig.write_image(save_svg_path)#, scale=scale)
    
    
    #return df_stacked


# In[4]:


#FIGURE 3A
plot_heights(temp_ht, x_col="DAP",
             save_html_path="../Figures/Fig3A.html",
             save_svg_path="../Figures/Fig3A.svg",
            )


# In[5]:


#SUP FIG 7
plot_heights(temp_ht, flower_DAP=flower_GDD, x_col="Cum_GDD [C]", xaxis_title="Cumulative Growing Degree Days [C]",
             save_html_path="../Figures/SupFig7.html",
             save_svg_path="../Figures/SupFig7.svg",
            )


# In[17]:


#FIGURE 3B
#fix ordering for figure
man_VCs["Trait_order"] = man_VCs["Trait"]
man_VCs["Trait_order"] = man_VCs["Trait_order"].replace({
    'Asymptote':1,
    'Growth rate':2, 
    'Inflection point':3,
    'YLD':4,
    'DTA':5,
    'DTS':6,
    'ASI':7,
    'PHT':8, 
    'EHT':9,
    'ASI':10,
    'Temporal plant height':11})
man_VCs["VC_order"] = man_VCs["Variance Components"]
man_VCs["VC_order"] = man_VCs["VC_order"].replace({'Hybrid':1,
                                                   'Row':5,
                                                   'Range':4,
                                                   'Rep':6,
                                                   'Residual':7,
                                                   'Hybrid:Flight':3,
                                                   'Flight':2})
man_VCs = man_VCs.sort_values(["VC_order", "Trait_order"], ascending=[False,True])

#setup metrics df
#CV_R2_Rep = man_VCs[["Repetability","Rsquared","CV"]].copy()
#CV_R2_Rep.index = man_VCs["Trait"]
#CV_R2_Rep.drop_duplicates(inplace=True)
#CV_R2_Rep = CV_R2_Rep.stack().reset_index().rename(columns={"level_1":"Metric",0:"Value"})
#CV_R2_Rep["Value"] = CV_R2_Rep["Value"]*100
CV_R2_Rep = man_VCs[["Trait","Repetability","Rsquared","CV"]].copy()
CV_R2_Rep.drop_duplicates(inplace=True)

#create figure
fig1 = px.bar(man_VCs, x="Trait", y="Explaned Percent Variation (%)",
             color="Variance Components", #pattern_shape="Variance Components",
             )
#fig2 = px.scatter(CV_R2_Rep, x="Trait", y="Value", symbol="Metric")

fig = make_subplots(specs=[[{"secondary_y": True}]])
cv = go.Scatter(x=CV_R2_Rep.Trait, y=CV_R2_Rep.CV, name="CV", mode='markers',
                marker=dict(size=12, color = 'white', symbol='square',
                            line=dict(width=1, color='black')
                           ),
                legendgroup="Metrics",
                legendgrouptitle_text="Metrics",
               )
rep = go.Scatter(x=CV_R2_Rep.Trait, y=CV_R2_Rep.Repetability, name="Repetability", mode='markers',
                marker=dict(size=12, color = 'white', symbol='triangle-down',
                            line=dict(width=1, color='black')
                           ),
                legendgroup="Metrics",
               )
r2 = go.Scatter(x=CV_R2_Rep.Trait, y=CV_R2_Rep.Rsquared, name="Rsquared", mode='markers',
                marker=dict(size=12, color = 'white', symbol='triangle-up',
                            line=dict(width=1, color='black')
                           ),
                legendgroup="Metrics",
               )
fig.add_traces([cv, rep, r2], secondary_ys=[True,True,True])
fig.add_traces(fig1.data)
fig.update_layout(barmode='stack')
#seperate legends into two groups
for trace_num in range(0, len(fig.data)):
    if fig.data[trace_num].name in man_VCs["Variance Components"].unique():
        fig.data[trace_num].legendgroup = "VC"
        fig.data[trace_num].legendgrouptitle = dict(text="Variance Components")
    #print(fig.data[trace_num].legendgroup)
    #print("\t",fig.data[trace_num].name)
    #print("\t", fig.data[trace_num].legendgrouptitle)
#fig.update_layout(legend_traceorder="reversed")

fig.update_yaxes(title_text="Explaned Percent Variation (%)", secondary_y=False)
fig.update_yaxes(title_text="Repeatability, Rsquard, and CV", secondary_y=True)
fig.update_layout(margin=dict(l=0, r=1, t=30, b=1))
fig.update_xaxes(tickangle=15)
fig.update_layout(autosize=False, width=800, height=400)

fig.write_html("../Figures/Fig3C.html")
fig.write_image("../Figures/FigC.svg")#, scale=scale)

fig.show()


# In[ ]:




