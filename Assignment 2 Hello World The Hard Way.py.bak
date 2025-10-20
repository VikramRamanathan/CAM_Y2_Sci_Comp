import numpy as np
import matplotlib.pyplot as plt

# Hello world, but in an unnecessarily difficult manner

user_input = input("Type in an input string, or type 0 to print Hello World!\n")

if user_input == "0":
    user_input = "Hello World!"
    
if len(user_input) < 3: # I cannot be bothered to handle these cases properly
    user_input = user_input.ljust(3)

y_values = [ord(i) for i in user_input]
print("\nASCII numbers of characters: " + str(y_values) + "\n")

# Encoding this string as points on a Chebyshev polynomial
degree = len(user_input)
x_values = np.linspace(-1, 1, degree)

eqn_mat = np.zeros((degree, degree))

eqn_mat[:, 0] = 1
eqn_mat[:, 1] = x_values
for i in range(2, degree): # Handy recurrence relation exists here
    eqn_mat[:, i] = 2 * x_values * eqn_mat[:, i-1] - eqn_mat[:, i-2]


# We now have an equation of the form eqn_mat * coefficients = y_values
inv = np.linalg.inv(eqn_mat)
coefficients = np.dot(inv, y_values)

print("Chebyshev polynomial basis coefficients: " + str(coefficients) + "\n")


# Chebyshev polynomial evaluation via recurrence relation
# This one took pen and paper and about 40 minutes to figure out, because math is hard
def evaluate_polynomial(coefficients, x):
    
    accumulator = 0
    
    coefficients = coefficients[::-1]
    highest_order = coefficients[0]
    second_highest_order = coefficients[1]
    
    
    for i in range(2, degree):
        reduction = highest_order
        highest_order = highest_order * (2 * x) + second_highest_order # reducing c1 by one step
        second_highest_order = coefficients[i] - reduction
    
    
    return highest_order * x + second_highest_order # final step is evaluating the base case from the construction


reconstructed_string = []
evaluation_points = []

# Recovering the message
for i in range(degree):
    x = i * 2 / (degree - 1) - 1
    evaluation_points.append(x)
    reconstructed_string.append(chr(round(evaluate_polynomial(coefficients, x))))
    
print("Evaluating polynomial at these normalised x-coordinates: " + str(evaluation_points) + "\n")

print("".join(reconstructed_string))

plotting_points = 100
x_for_plotting = np.linspace(-1, 1, plotting_points)
y_for_plotting = np.zeros((plotting_points))

for i in range(plotting_points):
    y_for_plotting[i] = evaluate_polynomial(coefficients, x_for_plotting[i])

plt.plot(x_for_plotting, y_for_plotting)
plt.title("The polynomial embedding the string")
plt.show()