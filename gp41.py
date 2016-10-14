
# coding: utf-8

# In[26]:

import numpy as np
import pandas as pd
import re
import sys



##########################################################
# read ALL_GP41 and filter only R13 Sample Time
##########################################################
def read_GP41():
    df = pd.read_csv('ALL_GP41_FINAL_Genbank.csv')
    print df.shape
    #print df.head(5)
    return df

##########################################################
# read ALL_GP41 and filter only R13 Sample Time
##########################################################
def filter_SampleTime(df, SampleTime):
    df_r = df[df["Sample Time"]==SampleTime]
    print df_r.shape

    #df_r.to_csv("GP41_" + SampleTime + ".csv", index=False)
    return df_r

def read_R13():
    df = pd.read_csv("GP41_R13.csv")
    print df.shape
    # print df['Sequence Name'].describe()
    # print df['Accession Number'].describe()
    # print df[df['Accession Number'].isnull()].shape
    return df

##########################################################
# filter records with couples info (Union ID notnull)
# find couple id and add couple column
##########################################################
def filter_with_couple(df, name):

    def find_my_couple(row):
        #print "find my couple"
        me = row["Sequence Name"]
        #print me
        ustr = row["Union ID 1"]
        #print ustr
        f = re.findall(r"([A-Z]{1}[0-9]{6})", ustr)
        #print f
        for item in f:
            #print item
            if(item != me):
                #print "my couple found!"
                return item

    df_c = df[pd.notnull(df["Union ID 1"])]
    df_c["Couple"] = df_c.apply(find_my_couple, axis=1)
    print df_c.shape
    df_c.to_csv(name+"With_Couple.csv", index=False)



##########################################################
#read the couples data
#filter out region and subtype from id, add columns
##########################################################
def read_couples():
    df_c = pd.read_csv('couples_gp41.csv', header=None)
    df_c.columns=["c1","c2"]
    print df_c.shape

    def split_id(row):
        c1str = row.c1
        c2str = row.c2
        c1f = re.findall(r"([A-Z]{1}[0-9]{6})_([0-9]+)_([A-Z])", c1str)[0]
        c2f = re.findall(r"([A-Z]{1}[0-9]{6})_([0-9]+)_([A-Z])", c2str)[0]
        return pd.Series([c1f[0],c1f[1],c1f[2],c2f[0],c2f[1],c2f[2]])


    newcols = df_c.apply(split_id, axis=1)
    newcols.columns = ["c1id","c1region","c1subtype","c2id","c2region","c2subtype"]

    newdf_c = df_c.join(newcols)

    newdf_c.to_csv("GP41_R13_couples.csv", index=False)


##########################################################
# read cluster Data
# filter out region and subtype from id, add columns
##########################################################
def read_cluster(clusterid):
    print clusterid
    clustername = "R13_gp41_georegion_subtype_hivtraceout"
    #clusterid = "02"
    clusterfilename = clustername + "_" + clusterid +".out"

    df_cls = pd.read_csv(clusterfilename)
    print df_cls.shape

    def split_seqid(row):
        c1str = row.SequenceID
        c1f = re.findall(r"([A-Z]{1}[0-9]{6})R13_([0-9]+)_([A-Z])", c1str)[0]
        return pd.Series([c1f[0],c1f[1],c1f[2]])


    newcols = df_cls.apply(split_seqid, axis=1)
    newcols.columns = ["seqid","region","subtype"]
    newdf_cls = df_cls.join(newcols)

    newfilename = "GP41_R13_cluster_" + clusterid + ".csv"
    newdf_cls.to_csv(newfilename, index=False)


##########################################################
# examine GP41_All DataFrame
# filter out records with unique Accession Number
##########################################################

def filter_multi_sequence_name(df):
    func_name = sys._getframe().f_code.co_name
    print func_name

    print df['Sequence Name'].describe()
    # NanSeqName = df[df['Sequence Name'].isnull()]
    # print NanSeqName.shape
    print '\n'


    MultiSample = df.groupby("Sequence Name").filter(lambda x: len(x) > 1)
    print MultiSample.shape
    print MultiSample["Sequence Name"].value_counts()

    # UniqueSample = df.groupby("Sequence Name").filter(lambda x: len(x) == 1)
    # print UniqueSample.shape


    MultiSample.reset_index(drop=True, inplace=True)
    MultiSample.to_csv("GP41_MultiSample.csv", index=True)


##########################################################
# pivot sequence
##########################################################
pivot_seq_df = pd.DataFrame({})#columns=['id','seq-pos','seq'])

def pivot_sequence(sample):

    f = open('GP41_R13_pivot_seq.csv', 'w')

    def make_pivot(arow):
        global pivot_seq_df
        aseq = arow['Sequence']
        #print len(aseq)
        new_df = pd.DataFrame({'id':arow['Sequence Name'],'seq-pos':range(len(aseq)) , 'seq':list(aseq)})

        #this solution write to file many times but consume less memory
        #change the file open mode to append "a" instead of "w"
        #new_df.to_csv(f, index=False, header=False)

        pivot_seq_df = pivot_seq_df.append(new_df)
        #return pivot_seq_df.shape


    print sample.shape
    sample.apply(make_pivot, axis=1)

    pivot_seq_df.to_csv(f, index=False)
    f.close()

    nrec = pivot_seq_df.shape[0]/390
    print nrec
