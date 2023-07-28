**INSTALL**

pip install git+https://github.com/cjun01/tetraploid_dosage.git

1. Introduction:
The calculate_probabilities function takes in two integers a and b, representing the number of reads of type 'A' and 'B' respectively. It then calculates the binomial probability for various dosage scenarios and prints out the results.

2. Setting up:
Ensure you have Python installed on your machine and import the necessary comb function from the math module.

3. Explanation of the Function:

  •	Heuristic Probabilities: These are pre-defined probabilities for different dosage scenarios. They are represented as:

  0:AAAA, 1:AAAB, 2:AABB, 3:ABBB, 4:BBBB
  
  •	Binomial Probability: An internal function binomial_probability calculates the binomial probability using the formula:
  ![formula](https://latex.codecogs.com/gif.latex?P%28X=k%29%20%3D%20%5Cbinom%7Bn%7D%7Bk%7D%20p%5Ek%20%281-p%29%5E%7Bn-k%7D)

  •	Special Cases: If there are no reads for 'A' or 'B', the function will handle these scenarios and assign probabilities accordingly.

  •	Calculation and Normalization: The function will calculate the probabilities for each scenario based on the heuristic and normalize the results to ensure they sum up to 1.

**Using the Function**

To use the function, simply call it with desired a and b values:

**from tetraploid_dosage import calculate_probabilities**


a, b = 1, 100
calculate_probabilities(a, b)
The function will then print out the probabilities for all dosage categories.

Expected Output:

The output will look like:

Dosage = 0: Probability = xx.xxxx
Dosage = 1: Probability = xx.xxxx
...
Where xx.xxxx represents the calculated probability.
