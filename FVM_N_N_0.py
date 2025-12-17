import numpy as np
import matplotlib.pyplot as plt
import time

# ============================================================
# Generation of Hexagonal Lattice Spaces for Graphene
# ============================================================

t0 = time.time()

# ---------------------------
# Parameters
# ---------------------------
n = 133
p_bond = 0.91     # bond probability (≈9% defects)

# MATLAB uses 1-based indexing, so we pad the array
M = np.zeros((n + 1, n + 1), dtype=int)

# ---------------------------
# Random percolation lattice
# ---------------------------
for i in range(1, n + 1):
    for j in range(1, n + 1):
        if np.random.rand() <= p_bond:
            M[i, j] = 1

# ---------------------------
# Hexagonal lattice geometry constraints
# ---------------------------
M[1, 1] = 0
M[3:n:3, 1:n:2] = 0
M[4:n:3, 2:n:2] = 0

M[:, :5] = 0
M[:4, :] = 0
M[:, n-2:n+1] = 0
M[n-2:n+1, :] = 0

# ============================================================
# Initial cluster labeling (Hoshen–Kopelman–type)
# ============================================================

k = 1

for i in range(2, n + 1):
    for j in range(2, n + 1):

        if M[i, j] == 1:

            # Case A
            if (j % 2 == 1) and (i % 3 == 2):
                up = M[i - 1, j]
                left = M[i, j - 1]

                if up == 0 and left == 0:
                    M[i, j] = k
                    k += 1
                elif up > 0 and left == 0:
                    M[i, j] = up
                elif up == 0 and left > 0:
                    M[i, j] = left
                else:
                    M[i, j] = min(up, left)

            # Case B
            if (j % 2 == 0) and (i % 3 == 2):
                if M[i, j - 1] == 0:
                    M[i, j] = k
                    k += 1
                else:
                    M[i, j] = M[i, j - 1]

            # Case C
            if (j % 2 == 1) and (i % 3 == 1):
                if M[i - 1, j - 1] == 0:
                    M[i, j] = k
                    k += 1
                else:
                    M[i, j] = M[i - 1, j - 1]

            # Case D
            if (j % 2 == 0) and (i % 3 == 0):
                if M[i - 1, j] == 0:
                    M[i, j] = k
                    k += 1
                else:
                    M[i, j] = M[i - 1, j]

# ============================================================
# Cluster relaxation (iterative label minimization)
# ============================================================

RM = np.zeros_like(M)
RRM = np.zeros_like(M)
Mmin = np.zeros_like(M)

L = True
while L:
    RM[:, :] = M
    RRM[:, :] = M
    Mmin[:, :] = M

    for i in range(2, n):
        for j in range(2, n):

            if M[i, j] > 0:

                neighbors = []

                if (j % 2 == 1) and (i % 3 == 2):
                    neighbors = [(i-1,j), (i,j-1), (i,j+1)]

                elif (j % 2 == 0) and (i % 3 == 2):
                    neighbors = [(i,j-1), (i+1,j), (i,j+1)]

                elif (j % 2 == 1) and (i % 3 == 1):
                    neighbors = [(i-1,j-1), (i+1,j), (i-1,j+1)]

                elif (j % 2 == 0) and (i % 3 == 0):
                    neighbors = [(i+1,j-1), (i-1,j), (i+1,j+1)]

                values = [RM[a,b] for a,b in neighbors if RM[a,b] > 0]
                if values:
                    mmin = min(values + [RM[i,j]])
                    RRM[i,j] = mmin
                    for a,b in neighbors:
                        if M[a,b] > 0:
                            RRM[a,b] = mmin

    changes = 0
    for i in range(1, n + 1):
        for j in range(1, n + 1):
            if M[i, j] > RRM[i, j]:
                M[i, j] = RRM[i, j]
                changes += 1

    L = changes > 0

# ============================================================
# Extract largest connected cluster
# ============================================================

labels, counts = np.unique(M[M > 0], return_counts=True)
largest_label = labels[np.argmax(counts)]

M[M != largest_label] = 0
M[M == largest_label] = 1

# ============================================================
# Hexagonal coordinates
# ============================================================

A = np.sqrt(3) / 2
X, Y = np.meshgrid(np.arange(n + 1), np.arange(n + 1))
X = A * X

Y = Y + np.tile([0, 0.5], (n + 1, (n + 1)//2 + 1))[:, :n+1]

X[M == 0] = 0
Y[M == 0] = 0

# ============================================================
# Boundary identification (value = 2)
# ============================================================

for i in range(2, n):
    for j in range(2, n):
        if M[i, j] == 1:
            for a, b in [
                (i-1,j),(i+1,j),(i,j-1),(i,j+1),
                (i-1,j-1),(i+1,j+1),(i-1,j+1),(i+1,j-1)
            ]:
                if M[a, b] == 0:
                    M[a, b] = 2

# ============================================================
# Count active lattice sites
# ============================================================

LN = np.sum(M == 1)

# ============================================================
# Save data
# ============================================================

np.savez(
    "gra_per_10500_9_exact.npz",
    M=M,
    LN=LN,
    X=X,
    Y=Y
)

# ============================================================
# Plot
# ============================================================

plt.figure(figsize=(6, 6))
plt.plot(X, Y, 'k.', markersize=5)
plt.axis('equal')
plt.axis([-1, 150, -1, 150])
plt.show()

print("LN =", LN)
print("Elapsed time:", time.time() - t0)
