import matplotlib.pyplot as plt

supported_tss_t120_file = "/supported_tss_t120.bed"
unsupported_tss_t120_file = "unsupported_tss_t120.bed"
supported_tss_t0_file = "/supported_tss_t0.bed"
unsupported_tss_t0_file = "/unsupported_tss_t0.bed"

# TSS 
def count_lines(file_path):
    with open(file_path, 'r') as file:
        return sum(1 for line in file)

# calculate
supported_tss_t120_count = count_lines(supported_tss_t120_file)
unsupported_tss_t120_count = count_lines(unsupported_tss_t120_file)
supported_tss_t0_count = count_lines(supported_tss_t0_file)
unsupported_tss_t0_count = count_lines(unsupported_tss_t0_file)

# 2. for graph
categories = ['Supported t120', 'Unsupported t120', 'Supported t0', 'Unsupported t0']
counts = [supported_tss_t120_count, unsupported_tss_t120_count, supported_tss_t0_count, unsupported_tss_t0_count]

# 3. grph with real count
plt.figure(figsize=(10,6))
plt.bar(categories, counts, color=['green', 'red', 'green', 'red'])

plt.title('TSS Supported vs Unsupported (Counts)')
plt.xlabel('TSS Type (t120 vs t0)')
plt.ylabel('Count')

plt.savefig('tss_supported_vs_unsupported_counts.png')  
plt.show()  

# 4. calculate ratios
total_t120 = supported_tss_t120_count + unsupported_tss_t120_count
total_t0 = supported_tss_t0_count + unsupported_tss_t0_count

supported_ratio_t120 = supported_tss_t120_count / total_t120
unsupported_ratio_t120 = unsupported_tss_t120_count / total_t120

supported_ratio_t0 = supported_tss_t0_count / total_t0
unsupported_ratio_t0 = unsupported_tss_t0_count / total_t0

ratios = [supported_ratio_t120, unsupported_ratio_t120, supported_ratio_t0, unsupported_ratio_t0]

# 5. graph for ratios
plt.figure(figsize=(10,6))
plt.bar(categories, ratios, color=['green', 'red', 'green', 'red'])

plt.title('TSS Supported vs Unsupported (Proportions)')
plt.xlabel('TSS Type (t120 vs t0)')
plt.ylabel('Proportion')

plt.savefig('tss_supported_vs_unsupported_proportions.png') 
plt.show() 
