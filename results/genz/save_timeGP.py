import pandas as pd


gp_time = pd.DataFrame([[0]*6 for i in range(6)], columns=["1","2","3","5","10","20"])

for genz in range(1, 7):
	
	c = 0
	for dim in [1,2,3,5,10,20]:
			

		file_name="GPresults{}dim{}.csv".format(genz, dim)
		time_complex = pd.read_csv(file_name)
		time_taken = time_complex.runtimeGP[0]
		print(time_taken)
		gp_time.iloc[genz-1][c] = time_taken
		c += 1

gp_time.to_csv("GPTime.csv")

		


