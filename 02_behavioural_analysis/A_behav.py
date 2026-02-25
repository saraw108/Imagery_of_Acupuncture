#!/usr/bin/env python
# coding: utf-8

# # IMACU
# bahavioural data analysis in python 
# 
# ## 1. install, load and assign

# In[1]:


import numpy as np
from numpy import genfromtxt
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import ttest_rel
import os
import statsmodels.api as sm
from statsmodels.formula.api import ols


# In[5]:


# "D:\IMACU\behav\IMACU_behav_final.xlsx"

df_raw_Acu = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/ExpectationStimulation_session1_Acu.csv')#get the values for a given file
df_raw_C1 = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/ExpectationStimulation_session1_C1.csv')#get the values for a given file
df_raw_C2 = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/ExpectationStimulation_session1_C2.csv')#get the values for a given file


# In[7]:


# Sess 1, before, Stim

FORMAT = ['soreness', 'aching', 'deep_pressure', 'heaviness', 'fullness', 'tingling', 'numbness', 'sharp_pain', 'dull_pain', 'warmness', 'cold', 'throbbing']
df_Acu = df_raw_Acu[FORMAT]
df_C1 = df_raw_C1[FORMAT]
df_C2 = df_raw_C2[FORMAT]


# In[8]:


df_Acu_clean = np.delete(df_Acu.values, [10], 0)
df_C1_clean = np.delete(df_C1.values, [10], 0)
df_C2_clean = np.delete(df_C2.values, [10], 0)

stim_acu = np.nan_to_num(df_Acu_clean) # sj 6 has some NaNs
stim_C1 = np.nan_to_num(df_C1_clean) 
stim_C2 = np.nan_to_num(df_C2_clean) 

plt.bar(np.arange(12)*2, np.mean(stim_acu,0), label = 'ExStimAcu', color =
'#7f0000')
plt.bar(np.arange(12)*2+0.3, np.mean(stim_C1,0), label = 'ExStimC1', color =
'#ff0000')
plt.bar(np.arange(12)*2+0.6, np.mean(stim_C2,0), label = 'ExStimC2', color =
'#ff7f7f')
plt.legend()

plt.ylim(0, 6)


# In[9]:


flat_data = np.append(stim_acu, stim_C1)
flat_data = np.append(flat_data, stim_C2)

dataframe = pd.DataFrame({'Points': np.repeat(['Acu', 'C1', 'C2'], 25*12),
                          'Items': np.tile(['soreness', 'aching', 'deep_pressure', 'heaviness', 'fullness', 'tingling', 'numbness', 'sharp_pain', 'dull_pain', 'warmness', 'cold', 'throbbing'], 25*3),
                          'intensity': flat_data})

model = ols('intensity ~ C(Points) + C(Items) + C(Points):C(Items)', data=dataframe).fit()
anova_result = sm.stats.anova_lm(model, type=2)
anova_result


# In[10]:


plt.boxplot((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1))) 
ex_stim = np.mean((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),0)
plt.ylim(0, 7)
plt.title("MASS expectation stim, collapsed over items")

save_three = np.round((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),2)
np.savetxt('stim_pre.txt', save_three)
save_three_mean = np.mean(np.round((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),2),1)
np.savetxt('MEAN_stim_pre.txt', save_three_mean)
weights_save = np.divide(save_three.T, save_three_mean)
np.savetxt('WEIGHTS_stim_pre.txt', weights_save.T)


# In[11]:


# Sess 1, before, Imag
df_raw_Acu = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/ExpectationImagery_session1_Acu.csv')#get the values for a given file
df_raw_C1 = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/ExpectationImagery_session1_C1.csv')#get the values for a given file
df_raw_C2 = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/ExpectationImagery_session1_C2.csv')#get the values for a given file

df_Acu = df_raw_Acu[FORMAT]
df_C1 = df_raw_C1[FORMAT]
df_C2 = df_raw_C2[FORMAT]


# In[12]:


df_Acu_clean = np.delete(df_Acu.values, [10], 0)
df_C1_clean = np.delete(df_C1.values, [10], 0)
df_C2_clean = np.delete(df_C2.values, [10], 0)

stim_acu = np.nan_to_num(df_Acu_clean) 
stim_C1 = np.nan_to_num(df_C1_clean) 
stim_C2 = np.nan_to_num(df_C2_clean) 


plt.bar(np.arange(12)*2, np.mean(stim_acu,0), label = 'ExImagAcu', color =
'#004471')
plt.bar(np.arange(12)*2+0.3, np.mean(stim_C1,0), label = 'ExImagC1', color =
'#0072bd')
plt.bar(np.arange(12)*2+0.6, np.mean(stim_C2,0), label = 'ExImagC2', color =
'#66aad7')
plt.legend()
plt.ylim(0, 6)


# In[13]:


flat_data = np.append(stim_acu, stim_C1)
flat_data = np.append(flat_data, stim_C2)

dataframe = pd.DataFrame({'Points': np.repeat(['Acu', 'C1', 'C2'], 25*12),
                          'Items': np.tile(['soreness', 'aching', 'deep_pressure', 'heaviness', 'fullness', 'tingling', 'numbness', 'sharp_pain', 'dull_pain', 'warmness', 'cold', 'throbbing'], 25*3),
                          'intensity': flat_data})

model = ols('intensity ~ C(Points) + C(Items) + C(Points):C(Items)', data=dataframe).fit()
anova_result = sm.stats.anova_lm(model, type=2)
anova_result


# In[14]:


plt.boxplot((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1))) 
ex_imag = np.mean((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),0)

save_three = np.round((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),2)
np.savetxt('imags_pre.txt', save_three)
save_three_mean = np.mean(np.round((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),2),1)
np.savetxt('MEAN_imags_pre.txt', save_three_mean)
weights_save = np.divide(save_three.T, save_three_mean)
np.savetxt('WEIGHTS_imags_pre.txt', weights_save.T)

plt.ylim(0, 7)
plt.title("MASS expectation imag, collapsed over items")


# In[15]:


# Sess 1, after, Stim

df_raw_Acu = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/Stimulation_session1_Acu.csv')#get the values for a given file
df_raw_C1 = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/Stimulation_session1_C1.csv')#get the values for a given file
df_raw_C2 = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/Stimulation_session1_C2.csv')#get the values for a given file

df_Acu = df_raw_Acu[FORMAT]
df_C1 = df_raw_C1[FORMAT]
df_C2 = df_raw_C2[FORMAT]


# In[16]:


df_Acu_clean = np.delete(df_Acu.values, [10], 0)
df_C1_clean = np.delete(df_C1.values, [10], 0)
df_C2_clean = np.delete(df_C2.values, [10], 0)

stim_acu = np.nan_to_num(df_Acu_clean) # sj 6 has some NaNs
stim_C1 = np.nan_to_num(df_C1_clean) 
stim_C2 = np.nan_to_num(df_C2_clean) 

plt.bar(np.arange(12)*2, np.mean(stim_acu,0), label = 'StimAcu', color =
'#7f0000')
plt.bar(np.arange(12)*2+0.3, np.mean(stim_C1,0), label = 'StimC1', color =
'#ff0000')
plt.bar(np.arange(12)*2+0.6, np.mean(stim_C2,0), label = 'StimC2', color =
'#ff7f7f')
plt.legend()
plt.ylim(0, 6)


# In[17]:


flat_data = np.append(stim_acu, stim_C1)
flat_data = np.append(flat_data, stim_C2)

dataframe = pd.DataFrame({'Points': np.repeat(['Acu', 'C1', 'C2'], 25*12),
                          'Items': np.tile(['soreness', 'aching', 'deep_pressure', 'heaviness', 'fullness', 'tingling', 'numbness', 'sharp_pain', 'dull_pain', 'warmness', 'cold', 'throbbing'], 25*3),
                          'intensity': flat_data})

model = ols('intensity ~ C(Points) + C(Items) + C(Points):C(Items)', data=dataframe).fit()
anova_result = sm.stats.anova_lm(model, type=2)
anova_result


# In[18]:


test1= ttest_rel(np.mean(stim_acu,1), np.mean(stim_C1,1))
print('t(24) value of acu vs c1: ', np.round(test1.statistic,2), 'p value of acu vs c1: ', np.round(test1.pvalue,2))
test2= ttest_rel(np.mean(stim_acu,1), np.mean(stim_C2,1))
print('t(24) value of acu vs c2: ', np.round(test2.statistic,2), 'p value of acu vs c2: ', np.round(test2.pvalue,2))
test3= ttest_rel(np.mean(stim_C1,1), np.mean(stim_C2,1))
print('t(24) value of c1 vs c2: ', np.round(test3.statistic,2), 'p value of c1 vs c2: ', np.round(test3.pvalue,2))


# In[19]:


plt.boxplot((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1))) 
stim = np.mean((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),0)

save_three = np.round((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),2)
np.savetxt('stim_post.txt', save_three)
save_three_mean = np.mean(np.round((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),2),1)
np.savetxt('MEAN_stim_post.txt', save_three_mean)
weights_save = np.divide(save_three.T, save_three_mean)
np.savetxt('WEIGHTS_stim_post.txt', weights_save.T)

plt.ylim(0, 7)
plt.title("MASS Stim, collapsed over items")


# In[20]:


# Sess 1, after, Imag

df_raw_Acu = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/Imagery_session1_Acu.csv')#get the values for a given file
df_raw_C1 = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/Imagery_session1_C1.csv')#get the values for a given file
df_raw_C2 = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/Imagery_session1_C2.csv')#get the values for a given file

df_Acu = df_raw_Acu[FORMAT]
df_C1 = df_raw_C1[FORMAT]
df_C2 = df_raw_C2[FORMAT]


# In[21]:


df_Acu_clean = np.delete(df_Acu.values, [10], 0)
df_C1_clean = np.delete(df_C1.values, [10], 0)
df_C2_clean = np.delete(df_C2.values, [10], 0)

stim_acu = np.nan_to_num(df_Acu_clean) 
stim_C1 = np.nan_to_num(df_C1_clean) 
stim_C2 = np.nan_to_num(df_C2_clean) 

plt.bar(np.arange(12)*2, np.mean(stim_acu,0), label = 'ImagAcu', color =
'#004471')
plt.bar(np.arange(12)*2+0.3, np.mean(stim_C1,0), label = 'ImagC1', color =
'#0072bd')
plt.bar(np.arange(12)*2+0.6, np.mean(stim_C2,0), label = 'ImagC2', color =
'#66aad7')
plt.legend()
plt.ylim(0, 6)


# In[22]:


flat_data = np.append(stim_acu, stim_C1)
flat_data = np.append(flat_data, stim_C2)

dataframe = pd.DataFrame({'Points': np.repeat(['Acu', 'C1', 'C2'], 25*12),
                          'Items': np.tile(['soreness', 'aching', 'deep_pressure', 'heaviness', 'fullness', 'tingling', 'numbness', 'sharp_pain', 'dull_pain', 'warmness', 'cold', 'throbbing'], 25*3),
                          'intensity': flat_data})

model = ols('intensity ~ C(Points) + C(Items) + C(Points):C(Items)', data=dataframe).fit()
anova_result = sm.stats.anova_lm(model, type=2)
anova_result


# In[23]:


plt.boxplot((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1))) 
imag = np.mean((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),0)

save_three = np.round((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),2)
np.savetxt('imags_post.txt', save_three)
save_three_mean = np.mean(np.round((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),2),1)
np.savetxt('MEAN_imags_post.txt', save_three_mean)
weights_save = np.divide(save_three.T, save_three_mean)
np.savetxt('WEIGHTS_imags_post.txt', weights_save.T)

plt.ylim(0, 7)
plt.title("MASS Imag, collapsed over items")


# In[24]:


# Sess 2, before, Imag

df_raw_Acu = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/ExpectationImagery_session2_Acu.csv')#get the values for a given file
df_raw_C1 = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/ExpectationImagery_session2_C1.csv')#get the values for a given file
df_raw_C2 = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/ExpectationImagery_session2_C2.csv')#get the values for a given file

df_Acu = df_raw_Acu[FORMAT]
df_C1 = df_raw_C1[FORMAT]
df_C2 = df_raw_C2[FORMAT]


# In[25]:


df_Acu_clean = np.delete(df_Acu.values, [10], 0)
df_C1_clean = np.delete(df_C1.values, [10], 0)
df_C2_clean = np.delete(df_C2.values, [10], 0)
print(df_Acu_clean.shape)

stim_acu = np.nan_to_num(df_Acu_clean) 
stim_C1 = np.nan_to_num(df_C1_clean) 
stim_C2 = np.nan_to_num(df_C2_clean) 

print(np.mean(stim_C2,0))

plt.bar(np.arange(12)*2, np.mean(stim_acu,0), label = 'Ses2 ExImagAcu', color =
'#004471')
plt.bar(np.arange(12)*2+0.3, np.mean(stim_C1,0), label = 'Ses2 ExImagC1', color =
'#0072bd')
plt.bar(np.arange(12)*2+0.6, np.mean(stim_C2,0), label = 'Ses2 ExImagC2', color =
'#66aad7')
plt.legend()
plt.ylim(0, 6)


# In[26]:


flat_data = np.append(stim_acu, stim_C1)
flat_data = np.append(flat_data, stim_C2)
print(25*12*3)
print(flat_data.shape)

dataframe = pd.DataFrame({'Points': np.repeat(['Acu', 'C1', 'C2'], 25*12),
                          'Items': np.tile(['soreness', 'aching', 'deep_pressure', 'heaviness', 'fullness', 'tingling', 'numbness', 'sharp_pain', 'dull_pain', 'warmness', 'cold', 'throbbing'], 25*3),
                          'intensity': flat_data})

model = ols('intensity ~ C(Points) + C(Items) + C(Points):C(Items)', data=dataframe).fit()
anova_result = sm.stats.anova_lm(model, type=2)
anova_result


# In[27]:


plt.boxplot((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1))) 
ex_imag2 = np.mean((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),0)

save_three = np.round((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),2)
np.savetxt('imags2_pre.txt', save_three)
save_three_mean = np.mean(np.round((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),2),1)
np.savetxt('MEAN_imags2_pre.txt', save_three_mean)
weights_save = np.divide(save_three.T, save_three_mean)
np.savetxt('WEIGHTS_imags2_pre.txt', weights_save.T)

plt.ylim(0, 7)
plt.title("MASS Expectation Imag Sess 2, collapsed over items")


# In[28]:


# Sess 2, before, Imag

df_raw_Acu = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/Imagery_session2_Acu.csv')#get the values for a given file
df_raw_C1 = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/Imagery_session2_C1.csv')#get the values for a given file
df_raw_C2 = pd.read_csv('C:/Users/mcnbf/Desktop/sara/IMACU/MASS/Imagery_session2_C2.csv')#get the values for a given file

df_Acu = df_raw_Acu[FORMAT]
df_C1 = df_raw_C1[FORMAT]
df_C2 = df_raw_C2[FORMAT]


# In[29]:


df_Acu_clean = np.delete(df_Acu.values, [10], 0)
df_C1_clean = np.delete(df_C1.values, [10], 0)
df_C2_clean = np.delete(df_C2.values, [10], 0)
print(df_Acu_clean.shape)

stim_acu = np.nan_to_num(df_Acu_clean) 
stim_C1 = np.nan_to_num(df_C1_clean) 
stim_C2 = np.nan_to_num(df_C2_clean) 

print(np.mean(stim_C2,0))

plt.bar(np.arange(12)*2, np.mean(stim_acu,0), label = 'Ses2 ImagAcu', color =
'#004471')
plt.bar(np.arange(12)*2+0.3, np.mean(stim_C1,0), label = 'Ses2 ImagC1', color =
'#0072bd')
plt.bar(np.arange(12)*2+0.6, np.mean(stim_C2,0), label = 'Ses2 ImagC2', color =
'#66aad7')
plt.legend()

plt.ylim(0, 6)


# In[30]:


flat_data = np.append(stim_acu, stim_C1)
flat_data = np.append(flat_data, stim_C2)
print(25*12*3)
print(flat_data.shape)

dataframe = pd.DataFrame({'Points': np.repeat(['Acu', 'C1', 'C2'], 25*12),
                          'Items': np.tile(['soreness', 'aching', 'deep_pressure', 'heaviness', 'fullness', 'tingling', 'numbness', 'sharp_pain', 'dull_pain', 'warmness', 'cold', 'throbbing'], 25*3),
                          'intensity': flat_data})

model = ols('intensity ~ C(Points) + C(Items) + C(Points):C(Items)', data=dataframe).fit()
anova_result = sm.stats.anova_lm(model, type=2)
anova_result


# In[31]:


plt.boxplot((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1))) 
imag2 = np.mean((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),0)

save_three = np.round((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),2)
np.savetxt('imags2_post.txt', save_three)
save_three_mean = np.mean(np.round((np.mean(stim_acu,1), np.mean(stim_C1,1), np.mean(stim_C2,1)),2),1)
np.savetxt('MEAN_imags2_post.txt', save_three_mean)
weights_save = np.divide(save_three.T, save_three_mean)
np.savetxt('WEIGHTS_imags2_post.txt', weights_save.T)

plt.ylim(0, 7)
plt.title("MASS Imag Sess 2, collapsed over items")


# In[32]:


plt.boxplot((ex_stim, ex_imag, stim, imag, ex_imag2, imag2)) 
plt.ylim(0, 7)
plt.title("All MASSes, collapsed over items and points")
to_save = np.round((ex_stim, ex_imag, stim, imag, ex_imag2, imag2),2)
np.savetxt('exstim_eximag_stim_imag_eximag2_imag2.txt', to_save)
to_save_mean = np.mean((ex_stim, ex_imag, stim, imag, ex_imag2, imag2),1)
np.savetxt('MEANs_exstim_eximag_stim_imag_eximag2_imag2.txt', to_save_mean)

weights_save = np.divide(to_save.T, to_save_mean)
np.savetxt('WEIGHTS_exstim_eximag_stim_imag_eximag2_imag2.txt', weights_save.T)


# ## lets look at the behav data from the log files

# In[36]:


answers = np.zeros([24,7,6,5])
path = 'D:/IMACU/Logs/logs_sorted/'

counter = 0
for s in np.arange(26):

    if s == 4 or s == 10:
        continue
    else:
        for r in np.arange(7):

            files = []
            for i in os.listdir(path):
                if os.path.isfile(os.path.join(path,i)) and 'newLogFile_subject_' + str(s+1) + '_sess_1_run_' + str(r+1) in i:
                    files.append(i)        

            log_file = pd.read_csv ('D:/IMACU/Logs/logs_sorted/' + files[0], sep = '\t')

            answers[counter,r,0,:] = log_file[log_file['TrialNr'].values == 1]['Response'].values
            answers[counter,r,1,:] = log_file[log_file['TrialNr'].values == 2]['Response'].values
            answers[counter,r,2,:] = log_file[log_file['TrialNr'].values == 3]['Response'].values
            answers[counter,r,3,:] = log_file[log_file['TrialNr'].values == 4]['Response'].values
            answers[counter,r,4,:] = log_file[log_file['TrialNr'].values == 5]['Response'].values
            answers[counter,r,5,:] = log_file[log_file['TrialNr'].values == 6]['Response'].values

        counter += 1


# In[37]:


print(np.shape(answers))

# mean over trials
first_mean = np.mean(answers,3, where = answers>0) # exclude zeros
print(np.shape(first_mean))

# mean over sessions
second_mean = np.mean(first_mean,1)
print(np.shape(second_mean))

# mean over subjects?
third_mean = np.mean(second_mean,0)
print(np.shape(third_mean))

sj_counts = np.zeros([24,6,5])
for j in np.arange(24):
    for c in np.arange(6):
        unique, counts = np.unique(answers[j,:,c,:], return_counts=True)
        sj_counts[j,c,unique.astype(np.int64)] = counts

plt.boxplot(second_mean) 
plt.ylim(0.5, 5)
plt.title("Intensity/Vividness as per button presses")


# In[38]:


flat_data = np.ndarray.flatten(first_mean[:,:,0:3])
print(24*7*3)
print(flat_data.shape)

dataframe = pd.DataFrame({'Conditions': np.tile(['StimAcu', 'StimC1', 'StimC2'], 24*7),
                          'Runs': np.tile(np.repeat(['run1', 'run2', 'run3', 'run4', 'run5', 'run6', 'run7'], 3), 24),
                          'rating': flat_data})

model = ols('rating ~ C(Conditions) + C(Runs) + C(Conditions):C(Runs)', data=dataframe).fit()
anova_result = sm.stats.anova_lm(model, type=2)
anova_result


# In[39]:


flat_data = np.ndarray.flatten(first_mean[:,:,3:7])
print(24*7*3)
print(flat_data.shape)

dataframe = pd.DataFrame({'Conditions': np.tile(['ImagAcu', 'ImagC1', 'ImagC2'], 24*7),
                          'Runs': np.tile(np.repeat(['run1', 'run2', 'run3', 'run4', 'run5', 'run6', 'run7'], 3), 24),
                          'rating': flat_data})

model = ols('rating ~ C(Conditions) + C(Runs) + C(Conditions):C(Runs)', data=dataframe).fit()
anova_result = sm.stats.anova_lm(model, type=2)
anova_result


# In[40]:


mean_stim_int = np.mean(second_mean[:,:3], 1)
mean_imag_int = np.mean(second_mean[:,3:], 1)

plt.boxplot((mean_stim_int,mean_imag_int)) 
plt.ylim(0.5, 5)
plt.title("Intensity/Vividness as per button presses, collapsed")


# In[41]:


np.mean(sj_counts[:,:,0],axis=0)

ans_opt = np.arange(6)
weight_counts = {
        "no answer": np.mean(sj_counts[:,:,0],axis=0),
        "1": np.mean(sj_counts[:,:,1],axis=0),
        "2": np.mean(sj_counts[:,:,2],axis=0),
        "3": np.mean(sj_counts[:,:,3],axis=0),
        "4": np.mean(sj_counts[:,:,4],axis=0),
}

width = 0.5

fig, ax = plt.subplots()
bottom = np.zeros(6)

for boolean, ans_count in weight_counts.items():
    p = ax.bar(ans_opt, ans_count, width, label=boolean, bottom=bottom)
    bottom += ans_count

ax.set_title("Number of answers per option, all subjects")
ax.legend(loc="upper right")

plt.show()


# In[42]:


collapsed = np.mean(np.mean(sj_counts[:,:3,:],axis=1),0)
print(collapsed.shape)
collapsed = np.array([collapsed, np.mean(np.mean(sj_counts[:,3:,:],axis=1),0)])
print(collapsed.shape)

ans_opt = np.arange(2)
weight_counts = {
        "no answer": collapsed[:,0],
        "1": collapsed[:,1],
        "2": collapsed[:,2],
        "3": collapsed[:,3],
        "4": collapsed[:,4],
}

width = 0.5

fig, ax = plt.subplots()
bottom = np.zeros(2)

for boolean, ans_count in weight_counts.items():
    p = ax.bar(ans_opt, ans_count, width, label=boolean, bottom=bottom)
    bottom += ans_count

ax.set_title("Number of answers per option, all subjects")
ax.legend(loc="upper right")

plt.show()


# In[43]:


for s in np.arange(24):
    print(sj_counts[s,:,:])

    ans_opt = np.arange(6)
    weight_counts = {
        "no answer": sj_counts[s,:,0],
        "1": sj_counts[s,:,1],
        "2": sj_counts[s,:,2],
        "3": sj_counts[s,:,3],
        "4": sj_counts[s,:,4],
    }

    width = 0.5

    fig, ax = plt.subplots()
    bottom = np.zeros(6)

    for boolean, ans_count in weight_counts.items():
        p = ax.bar(ans_opt, ans_count, width, label=boolean, bottom=bottom)
        bottom += ans_count

    ax.set_title("Number of answers per option")
    ax.legend(loc="upper right")

    plt.show()

