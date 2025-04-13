import os
import re

# Step 1: All .rst file paths collect karo
rst_files = []
for root, dirs, files in os.walk("."):
    for file in files:
        if file.endswith(".rst"):
            rst_files.append(os.path.join(root, file))

# Step 2: label regex se labels dhundo
label_pattern = re.compile(r'^\s*\.\.\s+_([^:]+):')
labels = {}

for file_path in rst_files:
    with open(file_path, 'r', encoding='utf-8') as f:
        for i, line in enumerate(f, start=1):
            match = label_pattern.match(line)
            if match:
                label = match.group(1)
                labels.setdefault(label, []).append((file_path, i))

# ✅ Step 3: Duplicate check aur print
duplicates = {k: v for k, v in labels.items() if len(v) > 1}

if not duplicates:
    print("✅ कोई duplicate label नहीं मिला — सब कुछ सही है!")
else:
    for label, occurrences in duplicates.items():
        print(f"Label: {label}")
        for file, line_num in occurrences:
            print(f"  → {file} — Line {line_num}")
        print()
