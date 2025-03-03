import seaborn as sns

# Run the solver and store values over time
A_history = [A0.copy()]
A = A0.copy()
for _ in range(100):  
    A = godunov_solver(A, dx, dt, 0.1)
    A_history.append(A.copy())

# Convert to numpy array for visualization
A_history = np.array(A_history)

# Create heatmap
plt.figure(figsize=(10, 6))
sns.heatmap(A_history, cmap="coolwarm", cbar=True)
plt.xlabel("River Distance")
plt.ylabel("Time Step")
plt.title("Flash Flood Evolution (Heatmap)")
plt.show()
