import os
import argparse
from sklearn.svm import OneClassSVM
import pandas as pd
import numpy as np
from collections import Counter
pd.options.mode.chained_assignment = None


def SVM(testevent,trainevent,poslist,KN,GA,NU,cols):
    Staticfile={}
    TestwithPred=pd.DataFrame()
    for b in poslist: #对每个坐标，提取出所有reads在该位点的值，判断两个群之间的差异
        Train=trainevent[trainevent['position']==b]
        Test=testevent[testevent['position']==b]
        if Test.shape[0] >10 and Train.shape[0] >10:
            #training of unmodified
            clf = OneClassSVM(kernel=KN,gamma=GA,nu=NU).fit(np.array(Train.iloc[:,cols]))
            #Predicting of modified
            Testout=clf.predict(np.array(Test.iloc[:,cols]))
            Test['Predict']=Testout
            TestwithPred=pd.concat([TestwithPred,Test],axis=0)
            #统计-1和1的个数
            d = Counter(Testout)
            d_s = sorted(d.items(),key=lambda x:x[1],reverse=True)
            #将结果输出到Staticfile字典
            if len(d_s)==2:
                Staticfile[b]=[d_s[1][1],Train['reference_kmer'].tolist()[1][0:1],Test.shape[0],Train.shape[0]]
            if len(d_s)==1:
                Staticfile[b]=[0,Train['reference_kmer'].tolist()[1][0:1],Test.shape[0],Train.shape[0]]
        else:
            Staticfile[b]=[0,'NULL',Test.shape[0],Train.shape[0]]
    return Staticfile,TestwithPred
def fill_ends(row):
    first_valid_index = int(row.first_valid_index())
    last_valid_index = int(row.last_valid_index())
    row.iloc[:first_valid_index] = 0
    row.iloc[last_valid_index + 1:] = 0
    return row
def transPredictedEventToMatrix(TestwithPred,poslist): #SVM步骤中产生的profile文件含有信号值为1和-1，miscalled碱基没有信号值，可以将这些值作为表达矩阵，用于reads按照结构的聚类
    """
    TestwithPred                -- Test event file with the prediction labels. Generated from SVM
    poslist                     -- The position list of reference.
    Translate the event file to vector. the mod bases were set to 1, unmod bases to 0, and indel bases to nan
    """    
    readnames=np.unique(TestwithPred['read_name']).tolist()
    matrix={}
    print(f'There are {len(readnames)} modified reads mapped on reference and generate the digital vectors ... ')
    for i in readnames:
        dfi=TestwithPred[TestwithPred['read_name']==i]
        unmodpos=dfi[dfi['Predict']==1]['position'].tolist()
        modpos=dfi[dfi['Predict']==-1]['position'].tolist()
        matrix[i]=[1 if (k in modpos) else (0 if (k in unmodpos) else np.nan) for k in poslist] 
    return pd.DataFrame(matrix).T 

def transUnmod_EventToMatrix(DMSOevent,poslist):
    """
    DMSOevent                -- Test event file with the prediction labels. Generated from SVM
    poslist                     -- The position list of reference.
    Translate the event file to vector. the called bases were set to 0, and the indel bases to nan
    """    
    matrix = DMSOevent.groupby(by=['read_name']).position.apply(list).to_dict()
    print(f'There are {len(matrix)} Unmodified reads mapped on reference and generate the digital vectors ...')
    vect={}
    for key,value in matrix.items():
        #if len(value)>len(poslist)*0.8:
        vect[key]=[0 if (k in value) else np.nan for k in poslist] 
    return pd.DataFrame(vect).T
def args():
    """
    predict modification of each bases by comparing the event file from Mod and Unmod group

    """ 
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--Mod_Event",   type=str, required=True,      help='The Modified collapsed event file.')
    parser.add_argument("-u", "--Unmod_Event", type=str, required=True,      help='The Unmodified collapsed event file.')
    parser.add_argument("-l", "--length",      type=int, required=True,      help='The RNA intermidiate length.')
    parser.add_argument("-c", "--columns",     type=int, nargs='+', default=[4,5,6],   help='the features used for SVM, 4: current, 5: current stdv, 6: dwell time. The default is [4,5,6].')
    parser.add_argument("-KN", "--kernel",     type=str, default='rbf',      help='The SVM parameter: kernel ，the default is (rbf).')
    parser.add_argument("-GA", "--gamma",      type=float, default=0.01,     help='The SVM parameter: gamma ，the default is (0.01).')
    parser.add_argument("-NU", "--nu",         type=float, default=0.01,     help='The SVM parameter: nu ，the default is (0.01).')
    args = parser.parse_args()
    return args


if __name__ == "__main__":

    args = args()

    if not os.path.exists(args.Mod_Event):
        print(f'Modified event file not found: {args.Mod_Event}')
        exit()
    if not os.path.exists(args.Unmod_Event):
        print(f'Unmodified event file not found: {args.Unmod_Event}')
        exit()
    assert args.columns in ([4,5], [4,5,6],[4,6]), "method should be one of [4,5], [4,5,6], [4,6]."

    #get the positions of specific sequence 
    poslist=[i for i in range(1,args.length+1)]
    # read in the Event files
    Mod_Event= pd.read_csv(args.Mod_Event,sep='\t',index_col=0)
    Unmod_Event= pd.read_csv(args.Unmod_Event,sep='\t',index_col=0)
    # predict every bases, add the results to the Mod_event file as a new column
    print(args.columns)
    Staticfile,TestwithPred=SVM(Mod_Event,Unmod_Event,poslist,args.kernel,args.gamma,args.nu,args.columns)
    #transform the profile to digital vectors,
    Modvect=transPredictedEventToMatrix(TestwithPred,poslist).apply(fill_ends, axis=1).fillna(2)
    print('Modified profile vectors generating succeed！')
    Unmodvect=transUnmod_EventToMatrix(Unmod_Event,poslist).apply(fill_ends, axis=1).fillna(2)
    print('Unmodified profile vectors generating succeed！')
    # write the results to files
    Statis_out=pd.DataFrame.from_dict(Staticfile, orient='index',columns=['Modifired_NO.','Kmer','Test_NO.','Train_NO.'])
    Statis_out.to_csv(args.Mod_Event+'.static',sep='\t',index=True)

    Modvect_out=args.Mod_Event+'.profile.vect'
    Modvect.to_csv(Modvect_out,sep='\t',header=True, index=True)

    Unmodvect_out=args.Unmod_Event+'.profile.vect'
    Unmodvect.to_csv(Unmodvect_out,sep='\t',header=True, index=True)

