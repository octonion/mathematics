import numpy as np
import math

def solve_for_p(a, b, max_solutions=10):
    """
    Finds the first positive integer (k,n) solutions for a given probability p=a/b.
    """
    if math.gcd(a, b) != 1:
        raise ValueError("a and b must be relatively prime.")

    c = b - a
    p = a / b
    print(f"--- Solving for p = {a}/{b} ---")

    # --- 1. Choose the Correct Seed ---
    # The threshold is (3-sqrt(5))/2
    threshold = (3 - 5**0.5) / 2
    
    if p > threshold:
        # "p large" seed from k=-2 case
        x0 = 3 * a - 2 * b
        y0 = a
        seed_type = "k=-2 (p > 0.382)"
    else:
        # "p small" seed from k=-1 case
        x0 = 3 * b - 2 * a
        y0 = -b
        seed_type = "k=-1 (p < 0.382)"

    x = x0
    y = y0
    k = (x - a - 8 * c) // (5 * c)
    n = k + (a + c * (k + 2) + y) // (2 * a)
    print("Starting k,n",k,n)

    print(f"Chosen Seed Type: {seed_type}")
    print(f"Starting Pell Solution (x₀, y₀): ({x0}, {y0})\n")
    
    # --- 2. Setup for Iteration ---
    # Pell multiplier matrix
    M = np.array([[9, 20], [4, 9]], dtype=object)
    
    # Initial vector for iteration
    v = np.array([x0, y0], dtype=object)
    
    positive_solutions = []
    m = 0
    
    # --- 3. Generate, Test, and Verify ---
    while len(positive_solutions) < max_solutions:
        if m > 0:
            v = M @ v # Generate next solution
        
        x, y = v[0], v[1]
        
        # Check divisibility for k
        if (x - a - 8 * c) % (5 * c) == 0:
            k = (x - a - 8 * c) // (5 * c)
            
            # Check divisibility for n (using the + sign for positive growth)
            if (a + c * (k + 2) + y) % (2 * a) == 0:
                n = k + (a + c * (k + 2) + y) // (2 * a)
                
                # Check if solution is positive and valid
                if k > 0 and n > 0 and n >= k + 2:
                    
                    # Verify the algebraic equation
                    lhs = c**2 * (k + 1) * (k + 2) + a * c * (k + 2) * (n - k)
                    rhs = a**2 * (n - k) * (n - k - 1)
                    
                    if lhs == rhs:
                        solution_data = {
                            "m": m,
                            "k": int(k),
                            "n": int(n),
                            "Verified": "Yes"
                        }
                        positive_solutions.append(solution_data)
                        print(f"Found Solution #{len(positive_solutions)} at m={m}: (k,n) = ({k:}, {n:})")
        m += 1
        if m > 1000: # Safeguard against extreme cases
            print("Iteration limit reached. No more solutions found.")
            break
            
    return positive_solutions

# --- Main Execution ---
if __name__ == "__main__":


    # Example: Run for p = 2/3
    solutions = solve_for_p(a=2, b=3)
    
    print("\n--- Final List of Solutions for p=2/3 ---")
    for sol in solutions:
        print(f"m={sol['m']}: (k, n) = ({sol['k']:}, {sol['n']:})")

    print()

    solutions = solve_for_p(a=2, b=7)
    
    print("\n--- Final List of Solutions for p=2/7 ---")
    for sol in solutions:
        print(f"m={sol['m']}: (k, n) = ({sol['k']:}, {sol['n']:})")

    print()
    
    solutions = solve_for_p(a=11, b=13)
    
    print("\n--- Final List of Solutions for p=11/13 ---")
    for sol in solutions:
        print(f"m={sol['m']}: (k, n) = ({sol['k']:}, {sol['n']:})")

    print()
    
    solutions = solve_for_p(a=5, b=41)
    
    print("\n--- Final List of Solutions for p=5/41 ---")
    for sol in solutions:
        print(f"m={sol['m']}: (k, n) = ({sol['k']:}, {sol['n']:})")

    print()
    
    solutions = solve_for_p(a=13, b=101)
    
    print("\n--- Final List of Solutions for p=13/101 ---")
    for sol in solutions:
        print(f"m={sol['m']}: (k, n) = ({sol['k']:}, {sol['n']:})")

    print()
    
    solutions = solve_for_p(a=91, b=101)
    
    print("\n--- Final List of Solutions for p=91/101 ---")
    for sol in solutions:
        print(f"m={sol['m']}: (k, n) = ({sol['k']:}, {sol['n']:})")

    print()

    solutions = solve_for_p(a=1, b=24)
    
    print("\n--- Final List of Solutions for p=1/24 ---")
    for sol in solutions:
        print(f"m={sol['m']}: (k, n) = ({sol['k']:}, {sol['n']:})")
