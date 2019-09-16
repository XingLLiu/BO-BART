import torch
import gpytorch as gp

class RBF_GP(gp.models.ExactGP):
    def __init__(self, train_x, train_y, likelihood):
        super(RBF_GP, self).__init__(train_x, train_y, likelihood)
        self.mean = gp.means.ConstantMean() 
        self.kernel = gp.kernels.RBFKernel()

    def forward(self, x):
        """forward pass of the GP

        """
        mean_x = self.mean(x)
        kernel_x = self.kernel(x)
        return gp.distribution.MultivariateNormal(mean_x, kernel_x)


