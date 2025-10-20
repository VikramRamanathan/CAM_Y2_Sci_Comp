import numpy as np
import scipy
from enum import Enum
import time
import re
import sys
import matplotlib.pyplot as plt

# The idea is to convolve pre-defined pdfs together

# This took like 8 hours, terrible decision

# With the optimisations I've done, you can find the pdf of 1 million 6-sided dice without too much difficulty. Don't use the Monte Carlo for that.

dice_dict = {}
MAX_FTT_SIZE = 100 # It can actually be much higher but I'm being safe

def convolve_n_pdfs(arrays):
    max_sequence_length = max([len(i) for i in arrays])
    mat = np.zeros((len(arrays), max_sequence_length))
    
    for i, array in enumerate(arrays):
        mat[i, :len(array)] = array
    
    # This breaks if I use rfft and it's due to some weird indexing behaviour I think
    transformed_convolution = np.fft.fft(mat, n = (sum([len(i) for i in arrays])) - len(arrays) + 1).prod(axis=0)
    convolution = np.fft.ifft(transformed_convolution).real
    
    # If it's one in a trillion who cares
    convolution[convolution < 10**-12] = 0
    
    size_before_trimming = convolution.shape[0]
    
    convolution = np.trim_zeros(convolution, trim="f")
    
    zero_offset = size_before_trimming - convolution.shape[0]
    
    convolution = np.trim_zeros(convolution, trim="b")
    
    return convolution, zero_offset
    
def chunk_list_up(arr, size):
    for i in range(0, len(arr), size):
        yield arr[i:i + size]
    

DICE_MODIFIER = Enum("DICE_MODIFIER", [("No_Modifier", 0), ("rr", 1), ("e", 2), ("kh", 3), ("kl", 4), ("kho", 5), ("klo", 6)])

# Creates a die characteristic array for an input die string
# This should really be an object but I'm lazy
# Format of a dice string is dX{rrN}{eM}[khLoI][klLoI]
# Square brackets are optional, curly brackets are any number of terms
# rrN means that die rolls of the number N will be rerolled once
# eM means die rolls of the number M will explode, i.e. count as an additional success
# khLoI and klLoI mean keep highest/lowest L dice results out of I dice total. This is the most complex feature.
def get_dice_characteristics(dice_string):
    
    def update_state(stored_op, num, keep_number):
        if stored_op == DICE_MODIFIER.No_Modifier:
            dice_dict["Dice sides"] = num
            
        elif stored_op == DICE_MODIFIER.rr:
            dice_dict["rr"].append(num)
            
        elif stored_op == DICE_MODIFIER.e:
            dice_dict["e"].append(num)
        
        elif stored_op == DICE_MODIFIER.kh:
            keep_number = num
            stored_op = DICE_MODIFIER.kho
        
        elif stored_op == DICE_MODIFIER.kl:
            keep_number = num
            stored_op = DICE_MODIFIER.klo

        elif stored_op == DICE_MODIFIER.kho:
            dice_dict["kho"] = (keep_number, num)

        elif stored_op == DICE_MODIFIER.klo:
            dice_dict["klo"] = (keep_number, num)
        
        if stored_op != DICE_MODIFIER.kh and stored_op != DICE_MODIFIER.kl:
            if dice_list[i] == "r":
                stored_op = DICE_MODIFIER.rr
            elif dice_list[i] == "e":
                stored_op = DICE_MODIFIER.e
            elif dice_list[i] == "k":
                if dice_list[i + 1] == "h":
                    stored_op = DICE_MODIFIER.kh
                elif dice_list[i + 1] == "l":
                    stored_op = DICE_MODIFIER.kl
        
        return dice_dict, stored_op, keep_number
    
    # If it's a die the first character is always d3
    dice_list = list(dice_string)[1:]
    
    dice_dict = {}
    dice_dict["rr"] = []
    dice_dict["e"] = []
    accumulator = []
    keep_number = None
    stored_op = DICE_MODIFIER.No_Modifier
    for i in range(len(dice_list)):
        if dice_list[i].isdigit():
            accumulator.append(dice_list[i])
        
        # Found an op
        else:
            
            if len(accumulator) == 0:
                continue
                
            else:
                num = int("".join(accumulator))
            
            dice_dict, stored_op, keep_number = update_state(stored_op, num, keep_number)
                
            accumulator = []
    
    num = int("".join(accumulator))
    dice_dict, _, _ = update_state(stored_op, num, keep_number)
    
    return dice_dict

def modify_dice_for_highest_and_lowest(prob_dist, keep_no, of_no, keep_high):
    
    def evaluate_partial_sum(order_n_add_one_sums, cdf, n, j):
        print("n:", n)
        print("j:", j)
        
        order_n_sums = []
        coeff = scipy.special.comb(n, j, exact=True)
        for i in range(len(order_n_add_one_sums)):
            order_n_sums.append(order_n_add_one_sums[i] + coeff * cdf[i] ** j * (1 - cdf[i]) ** (n-j))
        
        return order_n_sums
        
    def sum_highest_k_from_n_dice(prob_dist, cdf, keep_no, of_no, keep_high):
        dist_width = len(prob_dist)
        index_arr = np.ones(of_no).astype(int)
        out_pdf = np.zeros(keep_no * (dist_width - 1) + 1, dtype=np.float64)
        
        cdf = np.array(cdf, dtype=np.float64)
        prob_dist = np.array(prob_dist, dtype=np.float64)
        
        #print("Starting selection")
        start_time = time.time()
        #print(prob_dist)
        
        while index_arr[of_no - 1] < dist_width:
            prob_set = prob_dist[index_arr]
            
            # Should really hoist this check out of the loop but it doesn't matter
            
            if keep_high:
                highest_k_indices = np.argpartition(index_arr, -keep_no)[-keep_no:]
            
            else:
                highest_k_indices = np.argpartition(index_arr, keep_no)[:keep_no]
                
            value = index_arr[highest_k_indices].sum()
            probability = prob_set.prod()
            
            out_pdf[value] = out_pdf[value] + probability
            
            
            index_arr[0] = index_arr[0] + 1
            # Update the indexing array
            for i in range(of_no - 1):
                if index_arr[i] == dist_width:
                    index_arr[i] = 1
                    index_arr[i + 1] = index_arr[i + 1] + 1
        
        #print("Selection finished")
        #print("Time elapsed", time.time() - start_time)
        #print("Dist sums to", out_pdf.sum())
        
        # Excess precision unnecessary at the end
        return out_pdf.astype(np.float32)
    
    if keep_no >= of_no:
        pass # We are done
    
    cdf = [prob_dist[0]] # This should always be zero, but better to be safe
    for i in prob_dist[1:]:
        cdf.append(cdf[-1] + i)
    
    if keep_high and keep_no == 1:
        # This case is quite easy mathematically
        # CDF(max(X1, X2) < z) = CDF(X1 < z) * CDF(X2 < z)
        
        merged_cdf = [i**of_no for i in cdf]
        prob_dist = [0] + [merged_cdf[i+1] - merged_cdf[i] for i in range(len(merged_cdf)-1)]
        
    elif not keep_high and keep_no == 1:
        # CDF(min(X1, X2) > z) = CDF(X1 > z) * CDF(X2 > z)

        cdf_inverted = [1 - i for i in cdf]
        
        merged_cdf = [1 - i**of_no for i in cdf_inverted]
        prob_dist = [0] + [merged_cdf[i+1] - merged_cdf[i] for i in range(len(merged_cdf)-1)]
        
    elif keep_no == 0:
        return [0]
        
    elif keep_no > 1:
        if keep_no > 8:
            print("Please do not try and select multiple dice from more than 8, it won't end well. Falling back on no selection.")
            return prob_dist
        # I think doing this efficiently is an unsolved problem in combinatorics/probability
        # So here is a brute force solution
        
        prob_dist = sum_highest_k_from_n_dice(prob_dist, cdf, keep_no, of_no, keep_high)
        
        
        ## This is order statistics, using combinatorics I only marginally understand
        #dists_arr = [[i**of_no for i in cdf]]
        #for minus_j in range(1, keep_no):
        #    dists_arr.append(evaluate_partial_sum(dists_arr[minus_j - 1], cdf, of_no, of_no - minus_j))
        #
        #pdf_arr = []
        #for cdf_arr in dists_arr:
        #    pdf_arr.append([0] + [cdf_arr[i+1] - cdf_arr[i] for i in range(len(cdf_arr)-1)])
        #raise Exception("Crap I forgot that these aren't independant I can't just convolve them")
        #return pdf_arr
    
    return prob_dist

def generate_die(dice_string):
    dice_dict = get_dice_characteristics(dice_string)
    prob_dist = [0]
    prob_dist = prob_dist + [1/dice_dict["Dice sides"]] * dice_dict["Dice sides"]
    
    rerolls = dice_dict["rr"]
    if len(rerolls) != 0:
        prob_dist_copy = prob_dist[:]
        for i in rerolls:
            if i < len(prob_dist):
                prob_dist[i] = 0
        
        norm_factor = 1 - sum(prob_dist)
        for i in range(len(prob_dist)):
            prob_dist[i] = prob_dist[i] + norm_factor * prob_dist_copy[i]
    
    # Exploding dice don't play well with pdfs
    if len(dice_dict["e"]) != 0:
        print("My PDF-based methods don't play well with exploding dice, functionality currently unimplemented")
    
    if "kho" in dice_dict:
        prob_dist = modify_dice_for_highest_and_lowest(prob_dist, dice_dict["kho"][0], dice_dict["kho"][1], keep_high=True)
        
    elif "klo" in dice_dict:
        prob_dist = modify_dice_for_highest_and_lowest(prob_dist, dice_dict["klo"][0], dice_dict["klo"][1], keep_high=False)
        
    
    return prob_dist
    
# Returns the dice distribution, as well as the index of the "zero point" of the distribution
# If the input has no dice it will return [], None and it will be what you deserve
def parse_input_and_generate_dice(input_string, dice_dict, do_not_convolve = False):
    dice_regex = "\d+d"
    
    positive_dice = []
    negative_dice = []
    
    last_op = "+"
    input_arr = input_string.split(" ")
    for i in input_arr:
        
        if re.match(dice_regex, i):
            
            if last_op == None:
                print("Please correctly specify the operations between dice. Terminating.")
                sys.exit()
            
            if last_op == "+":
                positive_dice.append(i)
                
            else:
                negative_dice.append(i)
                
            last_op = None
        
        elif i == "+" or i == "-":
            last_op = i
        
    
    to_convolve_positive = []
    for i in positive_dice:
        num_dice = i.split("d")[0]
        dice_string = i[len(num_dice):]
        if dice_string not in dice_dict:
            dice_dict[dice_string] = generate_die(dice_string)
        
        to_convolve_positive = to_convolve_positive + [dice_dict[dice_string]] * int(num_dice)
    
    to_convolve_negative = []
    for i in negative_dice:
        num_dice = i.split("d")[0]
        dice_string = i[len(num_dice):]
        if dice_string not in dice_dict:
            dice_dict[dice_string] = generate_die(dice_string)
        
        to_convolve_negative = to_convolve_negative + [dice_dict[dice_string]] * int(num_dice)
    
    if do_not_convolve:
        return to_convolve_positive, to_convolve_negative
    
    # I run out of memory a lot if I don't do this
    max_fft_size_dynamic = MAX_FTT_SIZE
    max_fft_size_resizing_factor = 1 # Not using this anymore but hey might be useful at some point
    
    total_offset_positive = 0
    while len(to_convolve_positive) >= max_fft_size_dynamic:
        print("Distribution is massive, chunking iteratively")
        temp = []
        offsets = []
        for chunk in chunk_list_up(to_convolve_positive, max_fft_size_dynamic):
            convolved, zero_offset = convolve_n_pdfs(chunk)
            temp.append(convolved)
            offsets.append(zero_offset)
        to_convolve_positive = temp
        max_fft_size_dynamic = max(10, int(max_fft_size_dynamic * max_fft_size_resizing_factor))
        total_offset_positive = total_offset_positive + sum(offsets)
        print(total_offset_positive)
    
    
    max_fft_size_dynamic = MAX_FTT_SIZE
    total_offset_negative = 0
    while len(to_convolve_negative) >= max_fft_size_dynamic:
        print("Distribution is massive, chunking iteratively")
        temp = []
        offsets = []
        for chunk in chunk_list_up(to_convolve_negative, max_fft_size_dynamic):
            convolved, zero_offset = convolve_n_pdfs(chunk)
            temp.append(convolved)
            offsets.append(zero_offset)
        to_convolve_positive = temp
        max_fft_size_dynamic = max(10, int(max_fft_size_dynamic * max_fft_size_resizing_factor))
        total_offset_negative = total_offset_negative + sum(offsets)
    
    
    
    if len(to_convolve_positive) > 0:
        positive_convolution, positive_zero_offset = convolve_n_pdfs(to_convolve_positive)
        total_offset_positive = total_offset_positive + positive_zero_offset
        
        if len(to_convolve_negative) > 0:
            negative_convolution, negative_zero_offset = convolve_n_pdfs(to_convolve_negative)
            total_offset_negative = total_offset_negative + negative_zero_offset
            
            # Everywhere except here (ignoring clumping), the zero point is at the zero index.
            negative_convolution = negative_convolution[::-1]
            final_pdf, final_offset = convolve_n_pdfs([positive_convolution, negative_convolution])
            # The zero point moves from 0 to the max position in the negative
            
            zero_point = len(negative_convolution) - 1 - total_offset_positive + total_offset_negative - final_offset
            return final_pdf, zero_point
        
        return positive_convolution, -total_offset_positive
    
    if len(to_convolve_negative) > 0:
        negative_convolution, negative_zero_offset = convolve_n_pdfs(to_convolve_negative)
        total_offset_negative = total_offset_negative + negative_zero_offset
        
        return negative_convolution[::-1], len(negative_convolution) - 1 + total_offset_negative
    
    return [], None
    

def run_monte_carlo(positive_dice_pdfs, negative_dice_pdfs, num_trials):
    
    rng = np.random.default_rng()
    results = np.zeros(num_trials)
    for j in range(len(positive_dice_pdfs)):
        results[:] = results[:] + rng.choice(len(positive_dice_pdfs[j]), p = positive_dice_pdfs[j], size = num_trials)
        
    for j in range(len(negative_dice_pdfs)):
        results[:] = results[:] - rng.choice(len(negative_dice_pdfs[j]), p = negative_dice_pdfs[j], size = num_trials)
    
    return results


print("This program takes a configuration of dice and returns an exact probability distribution.")
print("It will then carry out a Monte Carlo simulation to estimate the probability distribution, and compare the two.")
print("Please enter a dice configuration you would like rolled.")
print("A few examples of valid dice configurations are listed below.\n")
print("1d6: Rolls a six-sided die")
print("2d6: Rolls two six-sided dice, summing the values")
print("2d6rr1rr2: Rolls two six-sided dice, rerolling values of one and two on each die")
print("1d6kh1o2: Rolls two six-sided dice, and keeps the highest result")
print("3d6kh1o2: Rolls three pairs of six-sided dice, and keeps the highest result from each pair")
print("1d6 + 1d3: Rolls a six-sided die and a three-sided die and adds the results")
print("1d6 - 1d3: Rolls a six-sided die and a three-sided die and subtracts the latter from the former")
print("10d13rr4kh2o3 - 10d5 - 2d3rr1: Exactly what it looks like")
print("\nInstruction codes are [rrX: Reroll values of X, khNoM: Keep Highest N rolls of M, klNoM: Keep Lowest N rolls of M]")
print("All operations (plus and minus) and dice must be space-separated")
print("\n")

input_string = input("Enter your dice configuration:\n")

# Input sanitisation is for losers
start_time = time.time()
positive_dice_pdfs, negative_dice_pdfs = parse_input_and_generate_dice(input_string, dice_dict, do_not_convolve=True)
pdf, zero_point = parse_input_and_generate_dice(input_string, dice_dict)
print(time.time() - start_time)
if len(pdf) == 0:
    print("Please enter at least one die")

else:
    print("The exact PDF is", pdf)
    print("(The zero of the distribution is at index", zero_point, "relative to the first entry in the PDF)")
    
    how_many_trials = "a"
    while isinstance(how_many_trials, str):
        how_many_trials = input("\nPlease enter how many Monte-Carlo trials to run, 0 to skip: ")
    
        try:
            how_many_trials = int(how_many_trials)
            
        except ValueError:
            print("That is not a number")
    
    if how_many_trials < 0:
        print("Invalid number entered, not running Monte Carlo.")
        
    elif how_many_trials == 0:
        print("Not running Monte Carlo")
    
    else:
        monte_carlo_results = run_monte_carlo(positive_dice_pdfs, negative_dice_pdfs, how_many_trials)
        plt.hist(monte_carlo_results, bins=np.arange(min(monte_carlo_results) - 0.5, max(monte_carlo_results) + 0.5), color="tab:orange", density=True, label="Monte Carlo")
        #plt.title("PDF of dice expression")
        #plt.xlabel("Value of dice expression")
        #plt.ylabel("Probability")
    
    pdf_max = max(pdf)
    x_values = []
    y_values = []
    start_recording = False
    last_above = 0
    
    # It is technically possible to make fiercely bimodal distributions, so I need to trunacate the distribution this way
    min_size_factor = 1000
    for i in range(len(pdf)):
        if pdf[i] >= pdf_max / min_size_factor:
            start_recording = True
            last_above = i
            
        if start_recording:
            x_values.append(i - zero_point)
            y_values.append(pdf[i])
    
    
    x_values = x_values[:len(x_values) - len(pdf) + last_above + 1]
    y_values = y_values[:len(y_values) - len(pdf) + last_above + 1]
    
    bar_width = 0.8
    if len(x_values) > 10:
        bar_width = 1
    plt.bar(x_values, y_values, width=bar_width, alpha = 0.7, label="Exact")
    plt.title("Exact PDF of dice expression")
    plt.xlabel("Value of dice expression")
    plt.ylabel("Probability")
    plt.legend()
    
    plt.show()

#print(get_dice_characteristics("d6rr1rr2e6e5kh1o2kl1o2"))
#print(generate_die("d6rr1kh1o2"))
#print(convolve_n_pdfs([d6, d6, d6]))