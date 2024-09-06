

#%% gene symbol convert
import pandas as pd

def EnsemblID2geneSymbol(EnsemblID: str = "") -> str:
    EnsemblID = EnsemblID.split(".")[0]  # 去掉尾部类似.1 .2 等的version号
    reference_dict = pd.read_csv(REFERENCE_DICT_DIR, sep='\t')
    df = reference_dict[reference_dict["Gene stable ID"] == EnsemblID]
    try:
        geneSymbol = df["Gene name"].tolist()[0]
        return geneSymbol
    except IndexError:  # 查找不到则返回None
        return None


def geneSymbol2EnsemblID(geneSymbol: str = "") -> str:
    reference_dict = pd.read_csv(REFERENCE_DICT_DIR, sep='\t')
    df = reference_dict[reference_dict["Gene name"] == geneSymbol]
    try:
        EnsemblID = df["Gene stable ID"].tolist()[0]
        return EnsemblID
    except IndexError:  # 查找不到则返回None
        return None

#%% data preprogress

def KNN_standardilze(X,z_score=False):
    #input dataframe (n_samples, n_features)
    import numpy as np
    import pandas as pd
    from sklearn.impute import KNNImputer
    from sklearn.preprocessing import StandardScaler
    X = pd.DataFrame(data=np.nan_to_num(X, copy=True,posinf=np.nan,neginf=np.nan,nan=np.nan),index=X.index,columns=X.columns)

    X1 = KNNImputer(n_neighbors=3).fit_transform(X)
    if z_score:
        scaler = StandardScaler().fit(X) 
        X1 = scaler.transform(X)
    X = pd.DataFrame(data=X1,index=X.index,columns=X.columns)
    return X

#%% univarante analysis

def Hedge_g(aval_ck, aval_case):
    MD = np.mean(aval_case) - np.mean(aval_ck)
    n1 = len(aval_ck)
    n2 = len(aval_case)
    SD1 = np.std(aval_ck)
    SD2 = np.std(aval_case)
    SDp = np.sqrt(((n1 - 1) * SD1**2 + (n2 - 1) * SD2**2) / (n1 + n2 - 2))
    g = (MD / SDp) * (1 - (3 / (4 * (n1 + n2) - 9)))
    return g


    #usage: uni_analysis(sam,injodr,'Sample.Type','Control','Case','t-test')
    import numpy as np
    import pandas as pd
    from scipy import stats
    from sklearn.metrics import roc_auc_score
    from statsmodels.stats.multitest import fdrcorrection
    
    ck_sample = list(injodr[injodr[group_col].isin([ck_name])].index.ravel())
    case_sample = list(injodr[injodr[group_col].isin([case_name])].index.ravel())
    DE_table = pd.DataFrame(data=None,index=sam.index,columns=['pval_median','pval_95CI-','pval_95CI+',
                                                               'FC_median','FC_95CI-','FC_95CI+',
                                                               'baseMean','log2FC','logPval',
                                                               'AUC','AUC_95CI-','AUC_95CI+','nar'])
   
    sam = sam.replace([np.inf, -np.inf], np.nan)
    for i,r in sam.iterrows():
        aval_ck = r[ck_sample].dropna().values
        aval_case = r[case_sample].dropna().values
        if len(aval_ck) >=3 and len(aval_case) >=3:
            if stat_func == 't-test':
                pval, pval_95CI_lb,pval_95CI_ub = uni_bootstrap_CI(aval_ck, aval_case,p_ttest)
            else:
                pval, pval_95CI_lb,pval_95CI_ub = uni_bootstrap_CI(aval_ck, aval_case,p_utest)
            if np.nanmean(aval_ck) != 0:
                FC,FC_95CI_lb,FC_95CI_ub = uni_bootstrap_CI(aval_ck, aval_case,fold_change)
                if FC > 0:
                    log2FC = np.log2(FC)
                else:
                    log2FC = 0
                
            else:
                FC = 0
                log2FC = 0
            baseMean = np.nanmean(r)

            logPval = -np.log10(pval)
            
            nar = (len(r[ck_sample+case_sample])-len(list(aval_ck)+list(aval_case)))/len(r[ck_sample+case_sample])
 
            y_true = np.array(list(np.ones((len(aval_case))))+list(np.zeros((len(aval_ck)))))
            y_pred = np.array(list(aval_case)+list(aval_ck))
            auc_mean = auc_roc(y_true,y_pred)
   
            auc_ci_lb, auc_ci_ub,auc_se = bootstrap_CI(y_true,y_pred, auc_roc)

            DE_table.loc[i,:]=[pval,pval_95CI_lb,pval_95CI_ub,
                               FC,FC_95CI_lb,FC_95CI_ub,
                               baseMean,log2FC,logPval,
                               auc_mean,auc_ci_lb, auc_ci_ub,nar]
    
    DE_table['p.adj'] = fdrcorrection(DE_table['pval_median'])[1]
    DE_table = DE_table[['baseMean','pval_median','pval_95CI-','pval_95CI+','p.adj','logPval',
                         'FC_median','FC_95CI-','FC_95CI+','log2FC','AUC','AUC_95CI-','AUC_95CI+','nar']]
    return DE_table

    #usage: pathway_AUCweight_sam(sam,injodr,'Sample.Type','Control','Case','p.adj','<',0.01,KEGG2ID,KEGG_info)
    import pandas as pd
    from scipy import stats
    from statsmodels.stats.multitest import fdrcorrection
    pathway_enrich_stat = pd.DataFrame(data=None,columns=['PathwayInfo','#Involved','#Total',
                                          'EnrichRatio','EnrichPvalue',
                                          'IncludeProtein'])
    for pw in pw_map:
        all_genes = pw_map[pw]
        intersection_genes = set(all_genes).intersection(set(used_genes))
        if len(intersection_genes) > 0:
            enrich = len(intersection_genes)/len(used_genes)/len(all_genes)*len(id2pw.keys())
            p = stats.fisher_exact([[len(id2pw.keys()),len(all_genes)],[len(used_genes),len(intersection_genes)]])[1]
            incl_pro = ';'.join(intersection_genes & set(used_genes))
            pathway_enrich_stat.loc[pw,:]=[pw_anno[pw],len(intersection_genes),len(all_genes),enrich,p,incl_pro]
    pathway_enrich_stat['Enrich_p.adj'] = fdrcorrection(pathway_enrich_stat['EnrichPvalue'])[1]
    pathway_enrich_stat = pathway_enrich_stat[['PathwayInfo','#Involved','#Total',
                                          'EnrichRatio','EnrichPvalue','Enrich_p.adj',
                                          'IncludeProtein']]
    return pathway_enrich_stat

#%% validation

from sklearn.metrics import confusion_matrix, classification_report, roc_auc_score, f1_score
import numpy as np
def bootstrap_CI(y_true,y_pred, func):
    from scipy.stats import bootstrap
    from scipy.stats import norm
    try:
        res = bootstrap((y_true,y_pred), func,vectorized=False, paired=True,n_resamples=1000)
        mean = np.nanmean(res.bootstrap_distribution)
        std = np.nanstd(res.bootstrap_distribution, ddof=1)
        ci_lower, ci_upper = norm.interval(0.95, loc=mean, scale=std)
        se = np.nanstd(res.bootstrap_distribution)
    except:
        ci_lower = ci_upper = se = np.nan
    return ci_lower, ci_upper, se


    from scipy.stats import bootstrap
    try:
        bootstrap = bootstrap((aval_ck, aval_case), func,vectorized=False,n_resamples=1000)
        median = np.median(bootstrap.bootstrap_distribution)
        ci_lower, ci_upper = bootstrap.confidence_interval
        se = bootstrap.standard_error
    except:
        median = ci_lower = ci_upper = se = np.nan  
    return median,ci_lower, ci_upper

def p_ttest(aval_ck, aval_case):
    from scipy import stats
    t, pval = stats.ttest_ind(aval_ck, aval_case)
    return pval
    
def p_utest(aval_ck, aval_case): 
    from scipy import stats
    t, pval = stats.mannwhitneyu(aval_ck, aval_case)
    return pval

def fold_change(aval_ck, aval_case):
    import numpy as np
    return np.mean(aval_case)/np.mean(aval_ck)

def sensitivity(y_true, y_pred):
    from sklearn.metrics import roc_curve
    fpr,tpr,thresholds = roc_curve(y_true,y_pred)
    best_thr = thresholds[np.argmax(np.abs(tpr-fpr))]
    if best_thr >= 1 or best_thr<= 0:
        best_thr = 0.5
    y_pred = [1 if i>=best_thr else 0 for i in y_pred]
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    if (tp + fn) != 0:
        return tp / (tp + fn)
    else:
        return 0
def specificity(y_true, y_pred):
    from sklearn.metrics import roc_curve
    fpr,tpr,thresholds = roc_curve(y_true,y_pred)
    best_thr = thresholds[np.argmax(np.abs(tpr-fpr))]
    if best_thr >= 1 or best_thr<= 0:
        best_thr = 0.5
    y_pred = [1 if i>=best_thr else 0 for i in y_pred]
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    if (tn + fp) != 0:
        return tn / (tn + fp)
    else:
        return 0
def ppv(y_true, y_pred):
    from sklearn.metrics import roc_curve
    fpr,tpr,thresholds = roc_curve(y_true,y_pred)
    best_thr = thresholds[np.argmax(np.abs(tpr-fpr))]
    if best_thr >= 1 or best_thr<= 0:
        best_thr = 0.5
    y_pred = [1 if i>=best_thr else 0 for i in y_pred]
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    if (tp + fp) != 0:
        return tp / (tp + fp)
    else:
        return 0
def npv(y_true, y_pred):
    from sklearn.metrics import roc_curve
    fpr,tpr,thresholds = roc_curve(y_true,y_pred)
    best_thr = thresholds[np.argmax(np.abs(tpr-fpr))]
    if best_thr >= 1 or best_thr<= 0:
        best_thr = 0.5
    y_pred = [1 if i>=best_thr else 0 for i in y_pred]
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    if (tn + fn) != 0:
        return tn / (tn + fn)
    else:
        return 0
def auc_roc(y_true, y_pred):
    auc = None
    try:
        auc = roc_auc_score(y_true, y_pred)
        if auc<0.5:
            y_true = 1-y_true
            auc = roc_auc_score(y_true, y_pred)
    except:
        auc = None
    return auc
def f1(y_true, y_pred):
    from sklearn.metrics import roc_curve
    fpr,tpr,thresholds = roc_curve(y_true,y_pred)
    best_thr =thresholds[np.argmax(np.abs(tpr-fpr))]
    if best_thr >= 1 or best_thr<= 0:
        best_thr = 0.5
    y_pred = [1 if i>=best_thr else 0 for i in y_pred]
    return f1_score(y_true, y_pred)
def explained_variance(y_true, y_pred):
    from sklearn.metrics import explained_variance_score
    return explained_variance_score(y_true, y_pred)

def LOO_CV(X,y,model,pred_mth='predict_proba'):
    import numpy as np
    from sklearn.model_selection import LeaveOneOut
    from sklearn.metrics import roc_curve
    loo = LeaveOneOut()
    loo.get_n_splits(X)
    y_pred = []
    for train_index, test_index in loo.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train = np.array(y)[train_index]
        model.fit(X_train, y_train)
        if pred_mth == 'predict':
            x_test_pred = model.predict(X_test)[0]
        elif pred_mth == 'predict_proba':
            x_test_pred = model.predict_proba(X_test)[0,1]
        y_pred.append(x_test_pred) 
    y_pred = np.array(y_pred) # list to array
    fpr, tpr, threshold = roc_curve(np.array(y), y_pred) 
    roc_auc = auc_roc(np.array(y), y_pred)  
 
    return {'AUC':roc_auc,
            'FPR_TPR_TH':(fpr,tpr,threshold),
            'Ytrue_Ypred':(np.array(y),y_pred)}

def LOO_CV_coef(X,y,X_col,model_name,model,pred_mth='predict_proba'):
    import pandas as pd
    import numpy as np
    from sklearn.model_selection import LeaveOneOut
    from sklearn.metrics import roc_curve
    loo = LeaveOneOut()
    loo.get_n_splits(X)
    y_pred = []
    feature_coef_stat = pd.DataFrame(data= None,index = X_col)
    feature_imp_stat = pd.DataFrame(data= None,index = X_col)
    i = 0
    for train_index, test_index in loo.split(X):
        X_train, X_test = X[train_index], X[test_index]
        y_train = np.array(y)[train_index]
        model.fit(X_train, y_train)
        
        if model_name in ['LASSO','Logistic','EN','SVM']:
            feature_coef = model.coef_
            feature_imp = np.abs(model.coef_)
            feature_imp /= feature_imp.sum()
            feature_coef_stat.loc[:,i]=feature_coef.ravel()
            feature_imp_stat.loc[:,i]=feature_imp.ravel()
        elif model_name in ['RF','XGB']:
            feature_coef = model.feature_importances_
            feature_imp = model.feature_importances_
            feature_coef_stat.loc[:,i]=feature_coef.ravel()
            feature_imp_stat.loc[:,i]=feature_imp.ravel()
        else:
            print('Feature weight is unavailable for [',model_name,'].')
        
        if pred_mth == 'predict':
            x_test_pred = model.predict(X_test)[0]
        elif pred_mth == 'predict_proba':
            x_test_pred = model.predict_proba(X_test)[0,1]
        y_pred.append(x_test_pred)
        i += 1
    y_pred = np.array(y_pred) # list to array
    fpr, tpr, threshold = roc_curve(np.array(y), y_pred) 
    roc_auc = auc_roc(np.array(y), y_pred)  

    return {'AUC':roc_auc,
            'FPR_TPR_TH':(fpr,tpr,threshold),
            'Ytrue_Ypred':(np.array(y),y_pred),
            'Coef_stat':(feature_coef_stat,feature_imp_stat)}


def kFold_CV(X,y,k,model,pred_mth='predict_proba'):
    import numpy as np
    from sklearn import model_selection 
    from sklearn.metrics import roc_curve

    if pred_mth == 'predict':
        y_pred = model_selection.cross_val_predict(model, X, y, cv=k, method='predict')
    elif pred_mth == 'predict_proba':
        y_pred = model_selection.cross_val_predict(model, X, y, cv=k, method='predict_proba')[:,1]
    y_pred = np.array(y_pred) # list to array
    fpr, tpr, threshold = roc_curve(np.array(y), y_pred) 
    roc_auc = auc_roc(np.array(y), y_pred)  

    return {'AUC':roc_auc,
            'FPR_TPR_TH':(fpr,tpr,threshold),
            'Ytrue_Ypred':(y,y_pred)}

def Split_Valid(X,y,model,test_size,pred_mth='predict_proba'):
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import roc_curve
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=42)
    clf = model.fit(X_train, y_train)
    if pred_mth == 'predict':
        y_pred = clf.predict(X_test)
        y_train_pred = clf.predict_proba(X_train)
    elif pred_mth == 'predict_proba':
        y_pred = clf.predict_proba(X_test)[:,1]
        y_train_pred = clf.predict_proba(X_train)[:,1]
    print('Train AUC:',auc_roc(np.array(y_train), y_train_pred))
    fpr, tpr, threshold = roc_curve(np.array(y_test), y_pred) 
    roc_auc = auc_roc(np.array(y_test), y_pred) 
    return {'AUC':roc_auc,
            'FPR_TPR_TH':(fpr,tpr,threshold),
            'Ytrue_Ypred':(y_test,y_pred)}

def Split_Valid_train_test(X_train, X_test, y_train, y_test,model,pred_mth='predict_proba'):
    #from sklearn.model_selection import train_test_split
    from sklearn.metrics import roc_curve
    #X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=test_size, random_state=42)
    clf = model.fit(X_train, y_train)
    if pred_mth == 'predict':
        y_pred = clf.predict(X_test)
        y_train_pred = clf.predict(X_train)
    elif pred_mth == 'predict_proba':
        y_pred = clf.predict_proba(X_test)[:,1]
        y_train_pred = clf.predict_proba(X_train)[:,1]
    print('Train AUC: %.3f' % auc_roc(np.array(y_train), y_train_pred))
    fpr, tpr, threshold = roc_curve(np.array(y_test), y_pred)
    roc_auc_trian = auc_roc(np.array(y_train), y_train_pred)
    roc_auc_test = auc_roc(np.array(y_test), y_pred)
    return {'AUC_Train':roc_auc_trian,'AUC_Test':roc_auc_test,
            'FPR_TPR_TH':(fpr,tpr,threshold),
            'Ytrue_Ypred':(y_test,y_pred)}
#%% model output
def save_feature_weight(X,y,model_name,model_obj,path,title):
    import os
    import numpy as np
    import pandas as pd
    if not os.path.exists(path):
        os.makedirs(path)  
    model_obj.fit(KNN_standardilze(X), y)
    if model_name in ['LASSO','Logistic','EN','SVM']:
        feature_coef = model_obj.coef_
        feature_imp = np.abs(model_obj.coef_)
        feature_imp /= feature_imp.sum()
        feature_importance_stat = pd.DataFrame({'Feature':X.columns,'Coefficient':feature_coef.ravel(),'Importance':feature_imp.ravel()})
        feature_importance_stat.to_csv(path+title+'.csv')
    elif model_name in ['RF','XGB']:
        feature_coef = model_obj.feature_importances_
        feature_imp = model_obj.feature_importances_
        feature_importance_stat = pd.DataFrame({'Feature':X.columns,'Coefficient':feature_coef.ravel(),'Importance':feature_imp.ravel()})
        feature_importance_stat.to_csv(path+title+'.csv')
    else:
        print('Feature weight is unavailable for [',model_name,'].',sep='')

def perf_table(y_actual, y_hat,path,output_csv):
    def perf_measure(y_actual, y_hat):
        import numpy as np
        TP = 0
        FP = 0
        TN = 0
        FN = 0
        for i in range(len(y_hat)): 
            if y_actual[i]==y_hat[i]==1:
                TP += 1
            if y_hat[i]==1 and y_actual[i]!=y_hat[i]:
                FP += 1
            if y_actual[i]==y_hat[i]==0:
               TN += 1
            if y_hat[i]==0 and y_actual[i]!=y_hat[i]:
               FN += 1
        if TP+FN != 0:
            TPR = TP/(TP+FN)# Sensitivity, hit rate, recall, or true positive rate 
        else:
            TPR = np.nan
        if TN+FP != 0:
            TNR = TN/(TN+FP) # Specificity or true negative rate
        else:
            TNR = np.nan
        if TP+FP != 0:
            PPV = TP/(TP+FP) # Precision or positive predictive value
        else:
            PPV = np.nan
        if TN+FN != 0:
            NPV = TN/(TN+FN) # Negative predictive value
        else:
            NPV = np.nan
        return[TP, FP, TN, FN,TPR,TNR,NPV,PPV]

    import pandas as pd
    theshold = np.flip((np.array(range(0,100,5))/100)[1:])
    perf_table = pd.DataFrame(data=None,index=theshold,columns=['TP','FP','TN','FN','Sensitivity','Specificity','NPV','PPV'])
    
    for th in theshold:
        new_y = []
        for y_p in y_hat:
            if y_p >= th:
                new_y += [1]
            else:
                new_y += [0]
        perf_table.loc[th,:]=perf_measure(y_actual, new_y)
    import os
    if not os.path.exists(path):
        os.makedirs(path)  
    perf_table.to_csv(path+output_csv)

#%% model parameter tuning

def model_tuning(X,y,model_name,model):
    param_grid_set = {'EN': {'alpha': [0.01, 0.1, 1, 10, 100],
                             'l1_ratio': [0.1, 0.3, 0.5, 0.7, 0.9]},
                      'LASSO': {'alpha': [0.01, 0.1, 1, 10]},
                      'Logistic': {'C': [0.01, 0.1, 1, 10], 
                                   'penalty': [None,'l1', 'l2']},
                      'SVM': {'C': [0.01, 0.1, 1, 10], 
                              #'kernel': ['linear', 'rbf', 'poly']
                              },
                      'RF': {'n_estimators': [10, 50, 100, 200],
                             'max_depth': [None, 5, 10, 20]},
                      'KNN': {'n_neighbors': [1, 3, 5, 7, 9], 
                              'p': [1, 2, 3]},
                      'XGB': {#'n_estimators': [50, 100, 200],
                              'max_depth': [3, 5, 7],
                              'learning_rate': [0.01, 0.1, 1],
                              'subsample': [0.5, 0.7, 0.9],
                              'colsample_bytree': [0.5, 0.7, 0.9]
                              },
                      'MLP': {'hidden_layer_sizes': [(50,50,50), (50,100,50), (100,)],
                              'activation': ['tanh', 'relu'],
                              'solver': ['sgd', 'adam'],
                              'alpha': [0.0001, 0.05],
                              'learning_rate': ['constant','adaptive'],
                              },
                      'HGBC': {'learning_rate': [0.05, 0.1, 0.2],
                               'max_depth': [3, 5, 7],
                               'max_leaf_nodes': [15, 31, 63],
                               }
                      }
    if model_name in param_grid_set.keys():
        from sklearn.model_selection import GridSearchCV
        clf = GridSearchCV(model, param_grid_set[model_name], ) #cv=10
        clf.fit(X, y)
        return (clf.best_estimator_,clf.best_score_,clf.best_params_)
    else:
        import numpy as np
        print('Parameter tuning of [',model_name,'] is not support.',sep='')
        return (model,np.nan,np.nan)

#%% multi models validation

def multi_models_validation(panel_set,summary_output):
    #panel_set = [('LN_RatioPanel_1',X_LN1,y_LN),('LN_RatioPanel_2',X_LN2,y_LN),('LN_RatioPanel_3',X_LN3,y_LN)]
    #input X: datatype: pandas.dataframe, (n samples, n features), no NA.
    #input y: datatype: pandas.dataframe, (n samples, 1) Control:0, Case:1
    import pandas as pd
    from sklearn import svm
    import sklearn.linear_model as lm
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.neural_network import MLPClassifier
    from xgboost import XGBClassifier
    from sklearn.ensemble import HistGradientBoostingClassifier
    #from scipy.stats import bootstrap
    from sklearn.model_selection import train_test_split

    import warnings
    warnings.filterwarnings('ignore')

    model1 = lm.Lasso(alpha=0.1)
    model2 = lm.LogisticRegression()
    model3 = lm.ElasticNet()
    model4 = svm.SVC(probability=True,kernel='linear')
    model5 = RandomForestClassifier(max_depth=7, random_state=0)
    model6 = KNeighborsClassifier(n_neighbors=5)
    model7 = XGBClassifier(random_state=42,eval_metric='logloss',use_label_encoder=False)
    model8 = MLPClassifier(random_state=1, max_iter=300)
    model9 = HistGradientBoostingClassifier()
    
    model_set = [('LASSO',model1),('Logistic',model2),('EN',model3),('SVM',model4),
                 ('RF',model5),('KNN',model6),('XGB',model7),('MLP',model8)]
    #model_set = [('HGBC',model9)]
    Valid_mth_set = ['LOO','5-Fold', '10-Fold',] #'TrainTest'

    ans = pd.DataFrame(data=None,columns=['Panel','Model','CV_method',
                                          'Threshold','TP','FN','FP','TN',
                                          'Sensitivity','Sensitivity_95CI-','Sensitivity_95CI+','Sensitivity_SE',
                                          'Specificity','Specificity_95CI-','Specificity_95CI+','Specificity_SE',
                                          'AUC','AUC_95CI-','AUC_95CI+','AUC_SE',
                                          'PPV','PPV_95CI-','PPV_95CI+','PPV_SE',
                                          'NPV','NPV_95CI-','NPV_95CI+','NPV_SE',
                                          'F1_Score','F1_95CI-','F1_95CI+','F1_SE',
                                          'Explained_Var_Score','Explained_Var_95CI-','Explained_Var_95CI+','Explained_Var_SE'
                                          ])

    for X_panel in panel_set:
        X = X_panel[1].values
        y = X_panel[2].values.ravel()
        
        '''
        train_sam, test_sam = train_test_split(X_panel[1].index,test_size=0.33,random_state=42)
        import os
        if not os.path.exists('TrainTestSplit/'):
            os.makedirs('TrainTestSplit/')  
        pd.DataFrame.from_dict({'TrainSet':train_sam,'Test_Set':test_sam},orient='index').T.to_csv('TrainTestSplit/'+X_panel[0]+'_TrainTestSampleID.csv')
        X_train = X_panel[1].loc[train_sam,:]
        y_train = np.array(X_panel[2].loc[train_sam,'Group'])
        #X_train = KNN_standardilze(X_train)
        X_test = X_panel[1].loc[test_sam,:]
        y_test = np.array(X_panel[2].loc[test_sam,'Group'])
        #X_test = KNN_standardilze(X_test)
        '''
        for model_panel in model_set:
            model,best_score,best_para = model_tuning(X,y,model_panel[0],model_panel[1])
            
            if model_panel[0] in ['LASSO','EN','XGB']:
                pred_mth = 'predict'
            else:
                pred_mth = 'predict_proba'
            for CV_method in Valid_mth_set:
                print(X_panel[0],model_panel[0],CV_method)
                if CV_method == 'LOO':
                    CV_rt = LOO_CV(X,y,model,pred_mth)
                    #CV_rt = LOO_CV_coef(X,y,model,pred_mth)
                elif CV_method == '5-Fold':
                    CV_rt = kFold_CV(X,y,5,model,pred_mth)
                elif CV_method == '10-Fold':
                    CV_rt = kFold_CV(X,y,10,model,pred_mth)
                #elif CV_method == 'TrainTest':
                #    CV_rt = Split_Valid2(X_train, X_test, y_train, y_test,model,pred_mth)
            
                y_true, y_pred = CV_rt['Ytrue_Ypred']
                import os
                if not os.path.exists('Ypred/'):
                    os.makedirs('Ypred/')  
                pd.DataFrame({'Sample':X_panel[1].index,
                              'Y_true':y_true,
                              'Y_pred':y_pred}).to_csv('Ypred/'+X_panel[0]+'_'+model_panel[0]+'_'+CV_method+'_Ypredict.csv')

                y_true = np.array(y_true)
                y_pred = np.array(y_pred)
                if roc_auc_score(y_true, y_pred)<0.5:
                    y_true = 1-y_true
                # Assuming you have two arrays: y_true and y_pred
                # y_true contains the true labels (0 or 1) and y_pred contains the predicted labels (0 or 1)
                # Calculate confusion matrix
                from sklearn.metrics import roc_curve
                fpr,tpr,thresholds = roc_curve(y_true,y_pred)
                best_thr = thresholds[np.argmax(np.abs(tpr-fpr))]
                if best_thr >= 1 or best_thr<= 0:
                    best_thr = 0.5
                y_bin = [1 if i>=best_thr else 0 for i in y_pred]
                tn, fp, fn, tp = confusion_matrix(y_true, y_bin).ravel()
                
                if tn == 0 or tp == 0:
                    print(X_panel[0]+' '+model_panel[0]+' '+CV_method,' FAILED.')
                    ans.loc[len(ans),:]=[X_panel[0],model_panel[0],CV_method,
                                         np.nan,tp,fn,fp,tn,
                                         np.nan,np.nan,np.nan,np.nan,
                                         np.nan,np.nan,np.nan,np.nan,
                                         np.nan,np.nan,np.nan,np.nan,
                                         np.nan,np.nan,np.nan,np.nan,
                                         np.nan,np.nan,np.nan,np.nan,
                                         np.nan,np.nan,np.nan,np.nan,
                                         np.nan,np.nan,np.nan,np.nan]
                else:
                 
                    save_AUC_fig(y_true,y_pred,'AUC_plot/',X_panel[0]+'_'+model_panel[0]+'_'+CV_method+'_ROC')
                    perf_table(y_true,y_pred,'PerformanceTable/',X_panel[0]+'_'+model_panel[0]+'_'+CV_method+'_ModelPerformance.csv')
                    save_feature_weight(X_panel[1],np.array(X_panel[2]),model_panel[0],model,'FeatureImportance/',X_panel[0]+'_'+model_panel[0]+'_'+CV_method+'FeatureImportance.csv')
                    dist_plot(y_true,y_pred,best_thr,'Dist_plot/',X_panel[0]+'_'+model_panel[0]+'_'+CV_method+'_FreqPlot')
                    
                    # Calculate sensitivity, specificity, PPV, NPV, F1 score, and AUC
                    sensitivity_value = sensitivity(y_true, y_pred)
                    specificity_value = specificity(y_true, y_pred)
                    ppv_value = ppv(y_true, y_pred)
                    npv_value = npv(y_true, y_pred)
                    f1_value = f1(y_true, y_pred)
                    auc_value = auc_roc(y_true, y_pred)
                    r2_value = explained_variance(y_true, y_pred)
                    # Calculate 95% confidence intervals using bootstrap
                    sensitivity_ci_lower, sensitivity_ci_upper, sensitivity_se = bootstrap_CI(y_true,y_pred, sensitivity)
                    specificity_ci_lower, specificity_ci_upper, specificity_se = bootstrap_CI(y_true,y_pred, specificity)
                    ppv_ci_lower, ppv_ci_upper, ppv_se = bootstrap_CI(y_true,y_pred, ppv)
                    npv_ci_lower, npv_ci_upper,npv_se = bootstrap_CI(y_true,y_pred, npv)
                    f1_ci_lower, f1_ci_upper,f1_se = bootstrap_CI(y_true,y_pred, f1)
                    auc_ci_lower, auc_ci_upper,auc_se = bootstrap_CI(y_true,y_pred, auc_roc)
                    r2_ci_lower, r2_ci_upper,r2_se = bootstrap_CI(y_true,y_pred, explained_variance)
                    
                    ans.loc[len(ans),:]=[X_panel[0],model_panel[0],CV_method,
                                         best_thr,tp,fn,fp,tn,
                                         sensitivity_value,sensitivity_ci_lower,sensitivity_ci_upper,sensitivity_se,
                                         specificity_value,specificity_ci_lower,specificity_ci_upper,specificity_se,
                                         auc_value,auc_ci_lower,auc_ci_upper,auc_se,
                                         ppv_value,ppv_ci_lower,ppv_ci_upper,ppv_se,
                                         npv_value,npv_ci_lower,npv_ci_upper,npv_se,
                                         f1_value,f1_ci_lower,f1_ci_upper,f1_se,
                                         r2_value,r2_ci_lower,r2_ci_upper,r2_se]
                    print('Validation AUC: %.3f' % auc_value)

    ans.to_csv(summary_output+'_model_performance.csv')
    print('File saved:',summary_output+'_model_performance.csv')
    return ans


    import os
    if not os.path.exists(path):
        os.makedirs(path)  
    import matplotlib.pyplot as plt
    import seaborn as sns
    import pandas as pd
    ans = pd.DataFrame({'Ytrue':y_true,'Ypred':y_pred})
    a = ans[ans['Ytrue']==0]['Ypred'].values.ravel()
    b = ans[ans['Ytrue']==1]['Ypred'].values.ravel()
    
    fig,ax = plt.subplots(figsize=(5,3))
    sns.set(style='dark')
    sns.set_style("dark", {"axes.facecolor": "#FFFFFF",
                           "axes.linewidth": 1,
                           "axes.edgecolor": 'black',
                           "xtick.color":'.15',
                           "ytick.color":'.15',
                           "xtick.major.size":5,
                           "ytick.major.size":5,})
    sns.distplot(a,
                   hist=True,
                   kde=True,#开启核密度曲线kernel density estimate (KDE)
                   bins=15,
                   kde_kws={'linestyle':'-','linewidth':'1','color':'#4682B4',#设置外框线属性                                               
                           },
                   color='#B0C4DE',
                   axlabel=r'Model Prediction',#设置x轴标题
                  )
    sns.distplot(b,
                   hist=True,
                   kde=True,#开启核密度曲线kernel density estimate (KDE)
                   bins=15,
                   kde_kws={'linestyle':'-','linewidth':'1','color':'#DC143C',#设置外框线属性                                               
                           },
                   color='#FFC0CB',
                   axlabel=r'Model Prediction',#设置x轴标题
                  )
    ax.vlines([threshold],0, ax.get_ylim()[1],linestyle='dashed',colors = 'Gray')
    plt.legend(labels=['Control','Case'], loc='best')
    plt.savefig(path+save_file_name+'.svg',format='svg',bbox_inches = 'tight'
#%%
import sys
import feature_selector as fs
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings('ignore')

#%%
sam1 = pd.read_csv('D:\\Product\\PA 20240818\\uni data\\APA.vs.BHA~.csv',index_col=0)
X1 = mp.KNN_standardilze(sam1.T).reset_index(drop=True)
y1 = pd.DataFrame([int(i.split('.')[0]) for i in sam1.columns.values],columns=['Group'])

panel_set = [('APA.vs.BHA',X1,y1)]
ans = mp.multi_models_validation(panel_set,'PA')
