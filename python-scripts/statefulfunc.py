#!/anaconda3/bin/python

running_sum = 0

def print_sum(added_value):
    global running_sum

    running_sum += added_value
    print(running_sum)

