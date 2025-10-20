import sys
import itertools

# Primitive truth table generator
# With another hour it would probably have a sane order of operations and support brackets via recursion
# But you did say 2.5 hours


sentence = input("Enter a boolean expression of the form \n\"a OR b AND NOT c\" \nSentences evaluated left to right (with NOT having priority), \nplease do not use brackets or I will be sad. \n\nBoolean operators in upper case only. \n\n")

sentence = sentence.split(" ")

# You know I could just use eval() here, but that would be lazy

symbols = [i for i in sentence if i != "AND" and i != "NOT" and i != "OR"]
# Deduplicate symbols that appear multiple times
symbols = [symbols[i] for i in range(len(symbols)) if symbols[i] not in symbols[:i]]
print("\nRecognised symbols are: " + str(symbols))

# Primitive validation

if sentence[0] == "AND" or sentence[0] == "OR":
    print("Invalid boolean expression, exiting program.")
    sys.exit()

for i in range(len(sentence)-1):
    if sentence[i] in symbols and sentence[i+1] in symbols \
        or sentence[i] == "AND" and sentence[i+1] == "OR" \
        or sentence[i] == "OR" and sentence[i+1] == "AND" \
        or sentence[i] == "NOT" and sentence[i+1] == "AND" \
        or sentence[i] == "NOT" and sentence[i+1] == "OR" \
        or sentence[i] == "OR" and sentence[i+1] == "OR" \
        or sentence[i] == "AND" and sentence[i+1] == "AND":
        print("Invalid boolean expression, exiting program.")
        sys.exit()


# This is terrible

input_values = []
truth_values = []
for i in range(2**len(symbols)):
    values = tuple(bin(i)[2:].zfill(len(symbols))[::-1])
    values = [i == "1" for i in values] # Convert to bool
    input_values.append(values)
    
    working_sentence = sentence[:]
    for i in range(len(symbols)):
        symbol = symbols[i]
        working_sentence[:] = [values[i] if x == symbol else x for x in working_sentence]
        
    
    if working_sentence[0] == "NOT" and isinstance(working_sentence[1], int):
        working_sentence[1] = not working_sentence[1]
        del working_sentence[0]
    
    
    # In hindsight it's quite easy to support brackets, I just need to implement this using a recursion-like structure
    # But that would push this project over 2 hours
    while len(working_sentence) > 1:
        
        
        if working_sentence[1] == "AND":
            if isinstance(working_sentence[2], bool):
                working_sentence[0] = working_sentence[0] and working_sentence[2]
                del working_sentence[2]
                del working_sentence[1]
            else: # If it's not a bool then it must be a not
                working_sentence[0] = working_sentence[0] and not working_sentence[3]
                del working_sentence[3]
                del working_sentence[2]
                del working_sentence[1]
        
        if len(working_sentence) > 1 and working_sentence[1] == "OR":
            if isinstance(working_sentence[2], bool):
                working_sentence[0] = working_sentence[0] or working_sentence[2]
                del working_sentence[2]
                del working_sentence[1]
            else:
                working_sentence[0] = working_sentence[0] or not working_sentence[3]
                del working_sentence[3]
                del working_sentence[2]
                del working_sentence[1]
    
    truth_values.append(working_sentence[0])

print("\nOutput")
print("Symbols: " + str(symbols) + ", Truth value: T/F\n")
for i, j in zip(input_values, truth_values):
    print(str(i) + " yields " + str(j))
    

if True not in truth_values:
    print("The expression is not satisfiable")
else:
    print("The expression is satisfiable")