#https://github.com/deepcharles/ruptures/issues/242
import ruptures as rpt
import numpy as np
rpt.version.version
data_list = [0.0, 0.1, 1.2, 1.0]
N_data = len(data_list)
data_mat = np.array(data_list).reshape(N_data,1)
algo = rpt.Binseg(model="l2",min_size=1,jump=1).fit(data_mat)
computed_break_dict={}
for n_bkps in range(N_data):
    try:
        result = algo.predict(n_bkps=n_bkps)
    except Exception as e:
        result = e
    computed_break_dict[n_bkps] = result
computed_break_dict
expected_break_dict = {
    0:[4],
    1:[2,4],
    2:[2,3,4],
    3:[1,2,3,4]
    }




#https://github.com/deepcharles/ruptures/issues/244
import ruptures as rpt
import numpy as np
rpt.version.version
data_list = [0,0.3,0.2,0.1, 10,11,12,13]
N_data = len(data_list)
data_mat = np.array(data_list).reshape(N_data,1)
algo = rpt.Binseg(model="normal",jump=1).fit(data_mat)
computed_break_dict={n_bkps:algo.predict(n_bkps=n_bkps) for n_bkps in range(4)}
computed_break_dict
expected_break_dict = {
    0:[8],
    1:[4,8],
    2:[4,6,8],
    3:[2,4,6,8]
    }

#from article.Rnw
data_list=[1, -7, 8, 10, 2, 4]
N_data = len(data_list)
data_mat = np.array(data_list).reshape(N_data,1)
algo = rpt.Binseg(model="l2",min_size=1,jump=1).fit(data_mat)
computed_break_dict={n_bkps:algo.predict(n_bkps=n_bkps) for n_bkps in range(4)}
computed_break_dict
expected_break_dict = {
    0:[8],
    1:[4,8],
    2:[4,6,8],
    3:[2,4,6,8]
    }
