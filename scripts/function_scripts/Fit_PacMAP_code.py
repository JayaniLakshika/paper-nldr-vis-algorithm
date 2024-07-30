def Fit_PaCMAP(nn_components, knn, init_n, MN_ratio_n, FP_ratio_n, path, path2):
  import numpy as np
  import pacmap
  from pandas import read_csv
  import pandas as pd
  import os
  


  d = read_csv(path)
 
  df = d.values


  # loading preprocessed coil_20 dataset
  # you can change it with any dataset that is in the ndarray format, with the shape (N, D)
  # where N is the number of samples and D is the dimension of each sample
  X = df.reshape(df.shape[0], -1)

  # initializing the pacmap instance
  # Setting n_neighbors to "None" leads to a default choice shown below in "parameter" section
  embedding = pacmap.PaCMAP(n_components = nn_components, n_neighbors = knn, MN_ratio = MN_ratio_n, FP_ratio = FP_ratio_n) 

  # fit the data (The index of transformed data corresponds to the index of the original data)
  X_transformed = embedding.fit_transform(X, init = init_n)
  
  r2 = nn_components
  li = ['PaCMAP']

  lst1 = [e for s in li for e in [s]*r2]
  lst2 = list(range(1, r2+1))
  lst2_n = list(map(str, lst2))
  
  test_list = [e+s for s in lst2_n for e in lst1]
  col_names = list(set(test_list))

  df_new = pd.DataFrame(X_transformed, columns = col_names)
  
  
  df_new.to_csv(path2,index=False)
