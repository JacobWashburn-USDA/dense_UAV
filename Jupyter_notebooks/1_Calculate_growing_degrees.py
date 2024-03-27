#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd


# In[2]:


#units = "F"
#up_thresh = 86 #F #30 #C
#dn_thresh = 50 #F #10 #C

units = "C"
up_thresh = 30 #C
dn_thresh = 10 #C

field_data = pd.read_csv( "../data/Masternote Washburn.csv")
sow_date = pd.to_datetime(field_data['Date Plot Planted [MM/DD/YY]'].unique()[0])
met_data = pd.read_csv("../data/2020_Bradford_MOH1_weather.csv")


# In[3]:


#calculate GDD Based on https://ndawn.ndsu.nodak.edu/help-corn-growing-degree-days.html

#convert to F
#(Celsius * 1.8) + 32
temp_data = met_data[["Date_key","Temperature [C]"]].dropna().reset_index(drop=True).copy()
if units == "F":
    temp_data["Temperature [F]"] = ((temp_data["Temperature [C]"]*1.8) + 32)
    temp_data.drop(columns="Temperature [C]", inplace=True)
    
#get max and min temps
temp_data["Date_time"] = pd.to_datetime(temp_data["Date_key"])
temp_data["Date"] = pd.to_datetime(temp_data["Date_time"].dt.date)
gdd = temp_data.pivot_table(index=["Date"], values = "Temperature ["+units+"]", aggfunc=['min',"max"])
gdd.columns = [x[0]+" "+x[1] for x in gdd.columns]

#if min and/or max < dn_thresh set to dn_thresh
gdd.loc[gdd["min Temperature ["+units+"]"] < dn_thresh, "min Temperature ["+units+"]"] = dn_thresh
gdd.loc[gdd["max Temperature ["+units+"]"] < dn_thresh, "max Temperature ["+units+"]"] = dn_thresh

#if min and/or max > up_thresh set to up_thresh
gdd.loc[gdd["min Temperature ["+units+"]"] > up_thresh, "min Temperature ["+units+"]"] = up_thresh
gdd.loc[gdd["max Temperature ["+units+"]"] > up_thresh, "max Temperature ["+units+"]"] = up_thresh

#calculate mean temp and subtract dn_thresh from it
gdd["mean temperature ["+units+"]"] = (gdd["max Temperature ["+units+"]"] + gdd["min Temperature ["+units+"]"])/2
gdd["GDD ["+units+"]"] = gdd["mean temperature ["+units+"]"] - dn_thresh

#calculate cummulative GDD
gdd = gdd.sort_index() # make sure in correct order
gdd["DAP"] = 1 #add DAP - days after planting
gdd["DAP"] = gdd["DAP"].cumsum()
gdd["Cum_GDD ["+units+"]"] = gdd["GDD ["+units+"]"].cumsum()

#save to file
gdd.to_csv("../data/GDDs_DAPs.csv")
gdd


# In[ ]:





# In[4]:


#calculate GDD to silk/pollen
flowering = field_data[['Pedigree', 'Plot', 'Pollen Date [MM/DD/YY]', 'Silk Date [MM/DD/YY]']].copy()
flowering["Pollen Date"] = pd.to_datetime(flowering["Pollen Date [MM/DD/YY]"])
flowering['Silk Date'] = pd.to_datetime(flowering['Silk Date [MM/DD/YY]'])
flowering.drop(columns=['Pollen Date [MM/DD/YY]', 'Silk Date [MM/DD/YY]'], inplace=True)
#remove NaN dates
flowering = flowering.dropna()
print(len(flowering))
#add in DAP and GDD
flowering = flowering.merge(gdd[["DAP", "Cum_GDD ["+units+"]"]], left_on="Pollen Date", right_index=True, how="left")
flowering = flowering.rename(columns={"DAP":"Pollen_DAP", "Cum_GDD ["+units+"]":"Pollen_Cum_GDD ["+units+"]"})
flowering = flowering.merge(gdd[["DAP", "Cum_GDD ["+units+"]"]], left_on="Silk Date", right_index=True, how="left")
flowering = flowering.rename(columns={"DAP":"Silk_DAP", "Cum_GDD ["+units+"]":"Silk_Cum_GDD ["+units+"]"})
print(len(flowering))

flowering.to_csv("../data/flowering_time_GDD_DAP.csv")
print(flowering[['Pollen_DAP',"Pollen_Cum_GDD ["+units+"]", 'Silk_DAP', "Silk_Cum_GDD ["+units+"]"]].min())
print(flowering[['Pollen_DAP',"Pollen_Cum_GDD ["+units+"]", 'Silk_DAP', "Silk_Cum_GDD ["+units+"]"]].max())
flowering


# In[5]:


mean_DAP = flowering[["Pollen_DAP", "Silk_DAP"]].values.flatten().mean()
std_DAP = flowering[["Pollen_DAP", "Silk_DAP"]].values.flatten().std()
print("Mean DAP: ", mean_DAP)
print("Flowering DAP, Mean +/- Std: ", mean_DAP-std_DAP,",", mean_DAP+std_DAP)

mean_GDD = flowering[["Pollen_Cum_GDD ["+units+"]", "Silk_Cum_GDD ["+units+"]"]].values.flatten().mean()
std_GDD = flowering[["Pollen_Cum_GDD ["+units+"]", "Silk_Cum_GDD ["+units+"]"]].values.flatten().std()
print("Mean Cum_GDD ["+units+"]: ", mean_GDD)
print("Flowering Cum_GDD ["+units+"], Mean +/- Std: ", mean_GDD-std_GDD,",", mean_GDD+std_GDD)


# In[ ]:




