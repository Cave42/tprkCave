import pandas as pd
import sys
import os
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('all_summary_stats')

parser.add_argument('metadata')

args = parser.parse_args()

meta = pd.read_csv(args.metadata)

df = pd.read_csv(args.all_summary_stats)

#df = pd.read_csv('./all_summary_stats.csv')

df['V_region_summary'] = ''
df['% input'] = ''

sampleNum = int(len(df.index)-1)

regionSum = [0]*sampleNum
percentInput = [0]*sampleNum

for i in range(len(df.index)-1):
    for j in range(7):
        regionSum[i] = regionSum[i] + df.iloc[i,2+j]
    percentInput[i] = df.iloc[i, 1] / regionSum[i]

for q in range(len(df.index)-1):
    df.iloc[q, 9] = regionSum[q]
    df.iloc[q, 10] = percentInput[q]

df.to_csv('More_summary_stats.csv', index = False, header = True)

df2 = pd.read_csv('./all_summary_stats.csv')

sampleNum = int(len(df2.index)-1)

regionSum = [0]*sampleNum
percentInput = [0]*sampleNum

for i in range(len(df.index)-1):
    for j in range(7):
        if(df2.iloc[i,2+j] < 5000):
            print(df2.iloc[i,2+j])
        else: df2.iloc[i,2+j] = ""

df2.to_csv('less_than_5000_all_summary_stats.csv', index = False, header = True)

#meta = pd.read_csv(args[2])

df['Rep Ratio'] = ''

sampleNum = int(len(meta.index))

for i in range(sampleNum):
    if(meta.iloc[i,4] == "A"):

        id = "Ill_" + (meta.iloc[i, 0])

        test = df[df['Sample'].str.contains(id)].index

        for j in range(sampleNum):
            if(meta.iloc[j,3] == meta.iloc[i,3] and meta.iloc[j,4] == 'B'):
                id2 = "Ill_" + (meta.iloc[j, 0])

                test2 = df[df['Sample'].str.contains(id2)].index

        ratio = int(df.iloc[test, 9])/int(df.iloc[test2, 9])

        df.iloc[test, 11] = ratio

        df.iloc[test2, 11] = ratio

df.to_csv('technicalRep_all_summary_stats.csv', index = False, header = True)

for k in range(sampleNum):
        if(df.iloc[k,11] < 0.8 or df.iloc[k,11] > 1.25):
            print("")
        else:
            for l in range(12):
                df.iloc[k,l] = ""

df.to_csv('technicalRepCheck_all_summary_stats.csv', index = False, header = True)
