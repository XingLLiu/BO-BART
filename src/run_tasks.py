import os
from subprocess import Popen
if __name__ == "__main__":

    print("Please indicate (space separated without break \
    \n[number of iterations]")
    args = int(input())
    num_iterations = args
    print(args)

    for genz in range(1, 7):

        for dim in [1, 2, 3, 5, 10, 20]:
 	    print("running", genz, dim)
	    directive = "Rscript GPRunTime.R {} {} {} &".format(dim, num_iterations, genz)

            os.system(directive)


