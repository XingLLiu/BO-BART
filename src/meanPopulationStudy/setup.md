# Mean Population experimental setup
We need to modify the kernel means and double integral of the kernel.

In `GPBQ.R`, we need to add another function called `computeGPBQEmpirical()` that takes in roughly the arguments of `computeGPBQ()` and addition ones such as the entire dataset `df`. When we compute empirical mean embedding, we compute the integral 

$$\int k(x,y) dy\approx \sum_{i=1}^n k(x, x_i)$$,

for $n$ being the size of the entire dataset (not just training). Similarly, for the double integral we have

$$\int\int k(x,y) dx dy\approx \sum_{i,j=1}^n k(x_i, x_j)$$.

Note the similarity to distribution regression. Having done this, the rest will follow from `computeGPBQ()` (make sure you don't sample from the uniform measure, but the candidate population set).