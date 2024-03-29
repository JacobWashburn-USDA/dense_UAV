#!/usr/bin/env python
# coding: utf-8

# In[1]:


import plotly.express as px
import plotly.graph_objects as go
import pandas as pd
from plotly.subplots import make_subplots


# In[20]:


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
man_calc_BLUPS = pd.read_csv("../data/blups_manual_phenos.csv", index_col=[0])

vcVARI = pd.read_csv("../data/VARI.varcomp.csv")
vcRmB = pd.read_csv("../data/RmB.varcomp.csv")


# In[3]:


triColors = ['#636efa','#EF553B','#00cc96']


# In[ ]:





# In[16]:


#plot heights by date and value
def plot_heights(df_stacked, flower_DAP=flower_DAP, x_col="DAP", value_col="Height", color="Tester",
                 testers_to_keep=['PHK76','PHP02','PHZ51'], xaxis_title = "Days After Planting (DAP)",
                 save_html_path="", save_svg_path="", fig_width=1300, fig_height=600):
    #filter testers
    for_fig = df_stacked[df_stacked["Tester"].isin(testers_to_keep)].copy()
    print("Number of Hybrids: ", len(df_stacked["Pedigree"].unique()), len(for_fig["Pedigree"].unique()))
    fig = px.box(for_fig, x = x_col, y=value_col, color=color,
                 color_discrete_sequence=triColors,
                 notched=True,
                )
    #fig.data[0].fillcolor = triColors[0]
    #fig.data[0].line = {"color":"black"}
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
    fig.update_layout(autosize=False, width=fig_width, height=fig_height, 
                      #font=dict(size=20)
                     )
    
    if save_html_path != "":
        fig.write_html(save_html_path)
    
    if save_svg_path != "":
        fig.write_image(save_svg_path)#, scale=scale)
    
    
    return fig


# In[17]:


#FIGURE 3A
fig3A = plot_heights(temp_ht, x_col="DAP",
                   save_html_path="../Figures/Fig3A.html",
                   save_svg_path="../Figures/Fig3A.svg",
                  )


# In[18]:


#SUP FIG 7
fig = plot_heights(temp_ht, flower_DAP=flower_GDD, x_col="Cum_GDD [C]", xaxis_title="Cumulative Growing Degree Days [C]",
             save_html_path="../Figures/SupFig7.html",
             save_svg_path="../Figures/SupFig7.svg",
            )


# In[14]:


#FIGURE 3B
import plotly.figure_factory as ff

def dist_plots():
    #create each individual figures
    figs_dict = {}
    for trait in man_calc_BLUPS["Trait"].unique():
        for_plot = []
        group_labels = []
        for tester in ["PHK76","PHP02","PHZ51"]: #man_calc_BLUPS["Tester"].unique():
            tmp = man_calc_BLUPS[(man_calc_BLUPS["Trait"]==trait) & (man_calc_BLUPS["Tester"]==tester)].copy()
            for_plot.append(tmp.Data.tolist())
            group_labels.append(tmp["Tester"].unique()[0])

        fig = ff.create_distplot(for_plot, group_labels, show_hist=False, show_rug=False, colors=triColors)
        fig["data"][0]['line']['dash'] = "dash"
        fig["data"][1]['line']['dash'] = "dot"
        fig["data"][2]['line']['dash'] = "dashdot"
        fig.update_layout(title=tmp["Trait"].unique()[0])
        figs_dict[trait] = fig

    #group figures
    fig = make_subplots(rows=3, cols=3,
                        vertical_spacing=0.09,
                        subplot_titles=list(figs_dict.keys()))
    row_col = [1,1]
    for key in figs_dict.keys():
        #print(row_col)
        #fix legend issues
        if row_col != [1,1]:
            for x in range(0,len(figs_dict[key].data)):
                figs_dict[key].data[x].showlegend=False

        fig.add_traces(figs_dict[key].data, rows=[row_col[0]]*3, cols=[row_col[1]]*3)
        #fig.add_trace(figs_dict[key].data[0], row=[row_col[0]], col=[row_col[1]])
        if row_col[1] < 3:
            row_col[1] = row_col[1]+1
        else:
            row_col[0] = row_col[0]+1
            row_col[1] = 1
    fig.data[0].legendgrouptitle = dict(text="Tester")

    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))
    width = 500
    height = 450
    fig.update_layout(autosize=False, width=width,height=height)
    fig.update_layout(legend=dict(yanchor="top",
                                  y=1.12,
                                  xanchor="left",
                                  x=0,
                                  #entrywidthmode="fraction",
                                  #entrywidth=1,
                                  orientation="h",
                                  #tracegroupgap=0,
                                  #itemwidth=40,
                                  #font_size=5,
                                 ))

    fig.write_html("../Figures/Fig3B.html")
    fig.write_image("../Figures/Fig3B.svg")#, scale=scale)

    fig.show()
    
    return fig
fig3B = dist_plots()


# In[ ]:





# In[15]:


#FIGURE 3C

def vc_plot(man_VCs):
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
                  color="Variance Components",
                  #pattern_shape="Variance Components",
                  color_discrete_sequence=px.colors.qualitative.Dark2,
                  #color_discrete_sequence=["#E69F00", "#56B4E9", "#009E73",
                  #                         "#F0E442", "#0072B2", "#D55E00", "#CC79A7"],
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
    fig.update_layout(autosize=False, width=800, height=450)

    fig.write_html("../Figures/Fig3C.html")
    fig.write_image("../Figures/Fig3C.svg")#, scale=scale)

    fig.show()
fig3C = vc_plot(man_VCs)


# In[21]:


#FIGURE 4
vcRmB


# In[ ]:





# In[ ]:




