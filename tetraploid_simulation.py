import random
from math import comb
import matplotlib.pyplot as plt

###############################################################################
# ORIGINAL calculate_probabilities FUNCTION (UNMODIFIED)
###############################################################################
def calculate_probabilities(a, b):
    """
    Original function: unchanged. Prints out probabilities for each dosage.
    0:AAAA, 1:AAAB, 2:AABB, 3:ABBB, 4:BBBB
    """
    e = 0.001
    p_values = {
        0: 1 - e,       # AAAA
        1: 3/4 + e,     # AAAB
        2: 2/4 + e,     # AABB
        3: 1/4 + e,     # ABBB
        4: 0 + e        # BBBB
    }

    def binomial_probability(n, k, p):
        return comb(n, k) * (p ** k) * ((1 - p) ** (n - k))

    n = a + b

    if a == 0 and b == 0:
        print("No reads.")
        return

    probs = {}
    for dosage, p in p_values.items():
        val = binomial_probability(n, a, p)
        # Basic safety checks (though typically unnecessary):
        if val > 1:
            val = 0
        if val < 0:
            val = 0
        probs[dosage] = val

    total_prob = sum(probs.values())
    if total_prob == 0:
        print("No reads.")
        return

    probabilities = {k: v / total_prob for k, v in probs.items()}

    for dosage, prob in probabilities.items():
        print(f"Dosage = {dosage}: Probability = {prob:.4f}")

    return probabilities

###############################################################################
# SIMULATION FUNCTIONS
###############################################################################
def simulate_one_genotype(true_dosage, coverage, e=0.001, num_reps=100):
    """
    For a given true dosage (0..4) and coverage (e.g., 10x),
    generate 'num_reps' sets of (ref_reads, alt_reads), each time
    calling 'calculate_probabilities' to see what dosage is inferred.

    The simulated data *does* incorporate an error rate 'e':
      p_Ref = fraction_of_A_in_genotype*(1-e) + (1 - fraction_of_A_in_genotype)*e
    Example: AAAA has fraction_of_A_in_genotype = 1.0
             so p_Ref = 1.0*(1-e) + (0.0)*e = 1-e
    This ensures AAAA can still yield a small % of alternate reads.
    """
    # fraction_of_A_in_genotype: AAAA=1.0, AAAB=0.75, AABB=0.50, ABBB=0.25, BBBB=0.0
    fraction_map = {
        0: 1.0,   # AAAA
        1: 0.75,  # AAAB
        2: 0.50,  # AABB
        3: 0.25,  # ABBB
        4: 0.0    # BBBB
    }
    frac_A = fraction_map[true_dosage]

    # Effective proportion of reference reads, incorporating error:
    p_ref = frac_A * (1 - e) + (1 - frac_A) * e

    correct_calls = 0

    for _ in range(num_reps):
        # Simulate the number of reference reads out of 'coverage'
        # with probability p_ref
        ref_reads = sum(random.random() < p_ref for _ in range(coverage))
        alt_reads = coverage - ref_reads

        # Call the original function
        probs = calculate_probabilities(ref_reads, alt_reads)
        if probs is None:
            continue

        # Best call = dosage with highest posterior probability
        best_call = max(probs, key=probs.get)
        if best_call == true_dosage:
            correct_calls += 1

    return correct_calls / num_reps if num_reps > 0 else 0.0

def main():
    # Coverage from 1x to 100x
    coverages = range(1, 101)
    num_reps = 20  # Fewer reps to reduce console spam; increase as needed
    error_rate = 0.001

    # Store accuracy results for each genotype:
    #   accuracy[dosage] = [acc_at_cov_1, acc_at_cov_2, ..., acc_at_cov_100]
    accuracy = {g: [] for g in range(5)}

    # Run simulations
    for cov in coverages:
        print(f"\n=== Coverage: {cov}x ===")
        for true_dosage in range(5):
            frac_correct = simulate_one_genotype(
                true_dosage, coverage=cov, e=error_rate, num_reps=num_reps
            )
            accuracy[true_dosage].append(frac_correct)
            print(f"  True Dosage={true_dosage}, Accuracy={frac_correct:.3f}")

    # Plot the results
    plt.figure(figsize=(10, 6))
    dosage_labels = {0: "AAAA", 1: "AAAB", 2: "AABB", 3: "ABBB", 4: "BBBB"}
    for g in range(5):
        plt.plot(coverages, accuracy[g], label=dosage_labels[g])

    plt.title("Accuracy of Dosage Calls (1–100x) with Error in Simulation")
    plt.xlabel("Coverage")
    plt.ylabel("Fraction of Correct Calls")
    plt.ylim(0, 1.05)
    plt.legend()
    plt.grid(True)
    
    # Save the plot to PNG at 600 dpi
    plt.savefig("dosage_accuracy.png", dpi=600)
    plt.show()

if __name__ == "__main__":
    main()
