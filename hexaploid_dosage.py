from math import comb
def calculate_probabilities(a, b):
    # Probabilities from heuristic
    #Ref dosage 0:aaaaaa, 1:aaaaab, 2:aaaabb, 3:aaabbb, 4:aabbbb, 5:abbbbb, 6:bbbbbb
    e=0.001
    p_values = {
        0: 1-e,  
        1: 5/6+e,
        2: 4/6+e,
        3: 3/6+e,
        4: 2/6+e,
        5: 1/6+e,
        6: 0+e  
    }


    # Calculate binomial probability
    def binomial_probability(n, k, p):
        return comb(n, k) * (p ** k) * ((1 - p) ** (n - k))

    n = a + b

    # Check for the special cases
    if a == 0 and b == 0:
        print("No reads.")
        return None
    else:
        # Compute probabilities for each scenario
        probs = {}
        for dosage, p in p_values.items():


            probs[dosage] = binomial_probability(n, a, p)

            if probs[dosage] > 1:
                probs[dosage] = 0
            if probs[dosage] < 0:
                probs[dosage] = 0


        # Normalize the probabilities
        total_prob = sum(probs.values())
        probabilities = {k: v / total_prob for k, v in probs.items()}

    # Check if there are reads or not
    if not probabilities:
        # If no reads, stop the program
        exit()

    # Print the probabilities for all dosage categories
    for dosage, prob in probabilities.items():
        print(f"Dosage = {dosage}: Probability = {prob:.4f}")



a,b=1,50
calculate_probabilities(a,b)
