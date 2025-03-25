import os
import argparse
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
pd.options.mode.chained_assignment = None
import scipy.stats
from statistics import mean
from scipy import stats
from collections import Counter
import umap
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from sklearn.cluster import DBSCAN



def get_index(lst=None, item=''):
    return [index for (index,value) in enumerate(lst) if value == item]

def Normal_SHAPE(inlist,windowsize=200,step=5): #采用icSHAPE-pipe的策略，在Windows内，top5%为1，最低5%为0
    """
    inlist                      -- A list of rate scores
    windowsize                  -- The sliding window size.
    step                        -- The sliding steps.
    Normalize the SHAPE reactivity scores 
    """    
    shape=[[] for i in range(len(inlist))]
    finalshape=[0 for i in range(len(inlist))]
    start=0
    while start+windowsize< len(inlist):
        sublist=np.array(inlist[start:start+windowsize])
        s95=stats.scoreatpercentile(sublist[~np.isnan(sublist)], 95)
        s5=stats.scoreatpercentile(sublist[~np.isnan(sublist)], 5)
        subshape=[1 if i>s95 else (0 if i<s5 else (i-s5)/(s95-s5)) for i in sublist]
        for i in range(windowsize):
            shape[start+i].append(subshape[i])
        start+=step
    #处理最后一个window
    endsublist=np.array(inlist[start:])
    s95=stats.scoreatpercentile(endsublist[~np.isnan(endsublist)], 95)
    s5=stats.scoreatpercentile(endsublist[~np.isnan(endsublist)], 5)
    subshape=[1 if i>s95 else (0 if i<s5 else (i-s5)/(s95-s5)) for i in endsublist]
    for i in range(len(endsublist)):
        shape[start+i].append(subshape[i])
    #对每个位置的值求平均
    for i in range(len(inlist)):
        finalshape[i]=mean(shape[i])
    return finalshape


def plotScore(shapescore,poslist,outfile):
    """
    shapescore                     -- A list of SHAPE reactivity score
    poslist                        -- A list of the positions.
    outfile                        -- The output figure.
    Plot the SHAPE reactivity scores 
    """    
    plt.figure(figsize=(10,4))
    plt.bar(poslist, shapescore, color =color_SHAPE(shapescore))
    plt.ylabel("Reactivity scores",fontsize=15) 
    plt.xlabel("Position",fontsize=15)
    plt.xlim(0,len(poslist))
    plt.savefig(outfile)
    plt.close()

def plotUMAP(UMAPdata,clust_ids,samplelabel,cluser_ratio,outfile):
    """
    UMAPdata                      --The UMAP dataframe in which include two columns(UMAP1 and UMAP2)
    clust_ids                     --A list of cluster labels of each reads. e.g [0,0,1,2,0,1,2]
    cluser_ratio                  -- The reads ratio of each cluster.
    outfile                       -- The output figure.
    Plot the UMAP clustering result 
    """    
    cluster_color = { 0: 'red', 1: '#FF7F00', 2: '#FF00FF',3:'#0000FF',4:'yellow',5:'black'}
    fig = plt.figure(figsize = (6,6))
    ax = fig.add_subplot(111)
    plt.scatter(UMAPdata['PC1'], UMAPdata['PC2'], c = [cluster_color[x]  if x in samplelabel else 'gray' for x in clust_ids], s = 20, alpha = 0.50 )
    ax.set_xlabel('UMAP1', fontsize = 15)
    ax.set_ylabel('UMAP2', fontsize = 15)
    N=len(cluser_ratio)
    for i in range(N):
        ax.text(0.05,0.95-0.05*i,'Cluter{}: {:.2%}'.format(i,cluser_ratio[i]),c=cluster_color[i], fontsize = 12, transform=ax.transAxes)
    ax.set_aspect('equal', 'box')
    plt.axis('equal')
    plt.savefig(outfile)
    plt.close()

def color_SHAPE(shape_list, cutoff=[0.3, 0.5, 0.7]):
    """
    shape_list              -- A list of SHAPE scores
    cutoff                  -- Cutoff of SHAPE color boundaries.
    
    Transform SHAPE values to color blocks
    """
    color_blocks=[]
    for value in shape_list:
        if value == 'nan':
            color_blocks.append('lightgray')
        else:
            shape = float(value)
            if shape < cutoff[0]:
                color_blocks.append('black')
            elif shape < cutoff[1]:
                color_blocks.append('blue')
            elif shape < cutoff[2]:
                color_blocks.append('orange')
            else:
                color_blocks.append('red')
    return color_blocks


def args():
    """
    Calculate and plot the reactivity scores according to different methods
    """ 
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--Mod_vector", type=str, required=True,      help='The Mod_Profile vector file')
    parser.add_argument("-u", "--Unmod_vector", type=str, required=True,     help='The Unmod_Profile vector file.')
    parser.add_argument("-o", "--output",      type=str, required=True,      help='The otput reactivity score file.')
    parser.add_argument("-s", "--method",      type=str, default='mean',     help='The reactivity score calculating methods. mean: calculate the mean score as one single structure; heter: calculate the alternative comformations separately. The default is (mean).')


    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = args()

    if not os.path.exists(args.Mod_vector):
        print(f'Mod_Profile vector file not found: {args.Mod_vector}')
        exit()
    if not os.path.exists(args.Unmod_vector):
        print(f'Unmod_Profile vector file not found: {args.Unmod_vector}')
        exit()
    assert args.method in ('mean', 'heter'), "method should be one of mean/heter"

    #read in the Profile vectors
    Modvect=pd.read_csv(args.Mod_vector,sep='\t',index_col=0)
    Unmodvect=pd.read_csv(args.Unmod_vector,sep='\t',index_col=0)

    #get the intermidiate length
    assert Modvect.shape[1]==Unmodvect.shape[1], "The lengths of modified and unmodifired vectors are not the same"
    poslist=[i for i in range(1,Modvect.shape[1]+1)]
    
    #calculate the scores according to the methods

    if args.method=='mean':
        NAI_one=np.sum(Modvect==1,axis=0)
        NAI_zero=np.sum(Modvect==0,axis=0)
        NAI_two=np.sum(Modvect==2,axis=0)

        DMSO_two=np.sum(Unmodvect==2,axis=0)
        DMSO_zero=np.sum(Unmodvect==0,axis=0)

        delta_rates=list(map(lambda x: x[0]-x[1] if x[0]>x[1] else 0, zip(NAI_two/(NAI_one+NAI_two+NAI_zero),DMSO_two/(DMSO_two+DMSO_zero))))
        rates=NAI_one/(NAI_one+NAI_zero)+delta_rates
        shapescore=Normal_SHAPE(rates)
        #Plot the PCA and Clustering results into PDF
        plotScore(shapescore,poslist,args.output+'.pdf')
        #Write the scores to outfile
        shapescoredit=pd.DataFrame(shapescore,index=poslist,columns=['score'])
        #shapescoredit.fillna('Null').to_csv(args.output+'.score',sep='\t')
        shapescoredit.fillna('NaN').to_csv(args.output+'.score',sep='\t')
        print("The mean Reactivity scores were calculated!")


    if args.method=='heter':
        # Classify reads using UMAP
        UP = umap(n_components=2)
        UP_out=pd.DataFrame(UP.fit_transform(Modvect),columns = ['UMAP1', 'UMAP2'])
        
        plt.scatter(UP_out['UMAP1'], UP_out['UMAP2'])
        plt.savefig(args.output+'.UMAP.pdf')
        plt.close()

        # Check the UMAP.PDF file and assign the cluster number and method
        N=input("Please check the UMAP.pdf file and then assign the Number of Clusers:  ")
        print(f"Number of Clusers is assigned:{N}!\n")

        M=input("Please choose the Clustering method, G for GaussianMixture and D for DBSCAN:  ")
        #assert M in ('G', 'D'), "method should be one of G/D"
        print(f"The cluser method is assigned:{M}!\n")

        if M=='G':
            clustering =GaussianMixture(n_components = int(N),covariance_type='full',random_state=1).fit(UP_out)
            cluster_ids=clustering.predict(UP_out).tolist()
            lists=Counter(cluster_ids).most_common(int(N))
            samplelabel={int(k):v for k, v in lists}.keys() 

        if M=='D':
            clustering =DBSCAN(eps=0.1, min_samples=int(N)).fit(UP_out)
            cluster_ids =clustering.labels_.tolist()
            lists=Counter(cluster_ids).most_common(int(N))
            samplelabel={int(k):v for k, v in lists}.keys() 

        #Calculate the scores of each cluster
        DMSO_two=np.sum(Unmodvect==2,axis=0)
        DMSO_zero=np.sum(Unmodvect==0,axis=0)

        cluser_ratio=[]
        cluser_score=[]
        cluser_name=[]
        for i in range(int(N)):
            cluser_ratio.append(cluster_ids.count(i)/len(cluster_ids))
            cluser_vect=Modvect.iloc[get_index(cluster_ids,i),:]

            NAI_one=np.sum(cluser_vect==1,axis=0)
            NAI_zero=np.sum(cluser_vect==0,axis=0)
            NAI_two=np.sum(cluser_vect==2,axis=0)
            delta_rates=list(map(lambda x: x[0]-x[1] if x[0]>x[1] else 0, zip(NAI_two/(NAI_one+NAI_two+NAI_zero),DMSO_two/(DMSO_two+DMSO_zero))))
            rates=NAI_one/(NAI_one+NAI_zero)+delta_rates
            shapescore=Normal_SHAPE(rates)            
            cluser_score.append(shapescore)
            cluser_name.append('Cluster'+str(i+1))
        print("Scores calculate successfully!")
        #Plot the PCA and Clustering results into PDF
        plotUMAP(UP_out,cluster_ids,samplelabel,cluser_ratio,args.output+'.UMAP.pdf')
        #Write the scores to outfile
        shapescoredit=pd.DataFrame(cluser_score,index=cluser_name,columns=poslist).T
        shapescoredit.fillna('Null').to_csv(args.output+'.scors',sep='\t')
        print("The heterogeneous reactivity scores were calculated!")

