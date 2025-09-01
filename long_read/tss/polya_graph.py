import matplotlib.pyplot as plt

supported_polya_p120_file = "/supported_polya_p120.bed"
unsupported_polya_p120_file = "/unsupported_polya_p120.bed"
supported_polya_p0_file = "/supported_polya_p0.bed"
unsupported_polya_p0_file = "/unsupported_polya_p0.bed"

# row caount
def count_lines(file_path):
    with open(file_path, "r") as file:
        return sum(1 for line in file)


supported_p120_count = count_lines(supported_polya_p120_file)
unsupported_p120_count = count_lines(unsupported_polya_p120_file)
supported_p0_count = count_lines(supported_polya_p0_file)
unsupported_p0_count = count_lines(unsupported_polya_p0_file)


categories = ['Supported p120', 'Unsupported p120', 'Supported p0', 'Unsupported p0']
counts = [supported_p120_count, unsupported_p120_count, supported_p0_count, unsupported_p0_count]

plt.figure(figsize=(12, 6))
plt.bar(categories, counts, color=['green', 'red', 'green', 'red'])
plt.title("PolyA Supported vs Unsupported (Counts)")
plt.xlabel("PolyA Type (p120 vs p0)")
plt.ylabel("Count")
plt.savefig("polya_supported_vs_unsupported_counts.png")
plt.show()

total_p120 = supported_p120_count + unsupported_p120_count
total_p0 = supported_p0_count + unsupported_p0_count

supported_ratio_p120 = supported_p120_count / total_p120
unsupported_ratio_p120 = unsupported_p120_count / total_p120
supported_ratio_p0 = supported_p0_count / total_p0
unsupported_ratio_p0 = unsupported_p0_count / total_p0

ratios = [supported_ratio_p120, unsupported_ratio_p120, supported_ratio_p0, unsupported_ratio_p0]


plt.figure(figsize=(12, 6))
plt.bar(categories, ratios, color=['green', 'red', 'green', 'red'])
plt.title("PolyA Supported vs Unsupported (Proportions)")
plt.xlabel("PolyA Type (p120 vs p0)")
plt.ylabel("Proportion")
plt.savefig("polya_supported_vs_unsupported_proportions.png")
plt.show()
