import matplotlib.pyplot as plt
from math import comb

###############################################################################
# ORIGINAL FUNCTION: calculate_probabilities (unchanged)
###############################################################################
def calculate_probabilities(a, b):
    # Probabilities from heuristic
    # 0:AAAA, 1:AAAB, 2:AABB, 3:ABBB, 4:BBBB
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

    # Print out the probabilities for all dosage categories
    for dosage, prob in probabilities.items():
        print(f"Dosage = {dosage}: Probability = {prob:.4f}")

    return probabilities

###############################################################################
# SIMULATION PARAMETERS AND READ GENERATION
###############################################################################
# Define ideal fractions for read generation without error:
# For genotype:
#   0 (AAAA): 1.0 (100% A, 0% B)
#   1 (AAAB): 0.75 (75% A, 25% B)
#   2 (AABB): 0.50 (50% A, 50% B)
#   3 (ABBB): 0.25 (25% A, 75% B)
# We omit genotype 4 (BBBB) because it is symmetric to AAAA.
ideal_fractions = {
    0: 1.0,   # AAAA (represents both AAAA and BBBB)
    1: 0.75,  # AAAB
    2: 0.50,  # AABB
    3: 0.25   # ABBB
}

# We'll simulate for these four genotype cases.
genotypes_to_simulate = [0, 1, 2, 3]

# Read depth (coverage) range: 1x to 100x
coverages = list(range(1, 101))

# Dictionary to store the inferred probability (from calculate_probabilities)
# for each genotype at each coverage.
# Format: probability_results[genotype] = [prob_at_depth1, prob_at_depth2, ..., prob_at_depth100]
probability_results = {geno: [] for geno in genotypes_to_simulate}

###############################################################################
# SIMULATION LOOP: For each read depth and for each genotype,
# generate ideal read counts and compute the inferred probability.
###############################################################################
for depth in coverages:
    for geno in genotypes_to_simulate:
        frac_A = ideal_fractions[geno]
        # Generate perfect read counts:
        # For example, for AAAB (geno 1) at a depth of 20, we aim for 75% A.
        ref_reads = int(round(depth * frac_A))
        alt_reads = depth - ref_reads

        # Call the original function (which applies the post-treatment error e=0.001)
        probs = calculate_probabilities(ref_reads, alt_reads)
        if probs is None:
            probability_results[geno].append(0)
        else:
            # Record the probability assigned to the true genotype (geno)
            probability_results[geno].append(probs.get(geno, 0))

###############################################################################
# DETERMINE THE MINIMUM COVERAGE WHERE EACH GENOTYPE HAS >= 99% INFERRED PROBABILITY
###############################################################################
# We want the minimum coverage at which, for all simulated genotypes, the inferred probability is at least 0.99.
threshold_coverage = None
for i, depth in enumerate(coverages):
    # Check if all four genotypes have probability >= 0.99 at this coverage
    if all(probability_results[geno][i] >= 0.99 for geno in genotypes_to_simulate):
        threshold_coverage = depth
        break

###############################################################################
# COMPOSITE PLOT: Plot the Inferred Probability vs. Coverage for Each Genotype
###############################################################################
plt.figure(figsize=(10, 6))
genotype_labels = {
    0: "AAAA/BBBB",  # representing AAAA and BBBB
    1: "AAAB",
    2: "AABB",
    3: "ABBB"
}

for geno in genotypes_to_simulate:
    plt.plot(coverages, probability_results[geno], marker='o', label=genotype_labels[geno])

# Mark the threshold coverage with a vertical line if found.
if threshold_coverage is not None:
    plt.axvline(x=threshold_coverage, color='red', linestyle='--',
                label=f"Threshold: {threshold_coverage}x (all ≥99%)")

plt.xlabel("Coverage (read depth)")
plt.ylabel("Inferred Probability of True Genotype")
plt.title("Inferred Genotype Probability vs. Coverage")
plt.ylim(0, 1.05)
plt.legend()
plt.grid(True)
plt.savefig(r"C:\Users\CFLZXC\OneDrive - Plant and Food Research\potato\potato_bioinformatics\results\ready\haplotype\dosage_accuracy.png", dpi=600)

plt.show()
