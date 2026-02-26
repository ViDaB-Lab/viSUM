import csv
import statistics
import sys

# script compares scores for realms and produces csv files: 1. a check csv file that contains all values and the results to make sure the script is working ok and 2. a csv file with results with transcript ID, length, highest score, realm prediction, and higher or lower than 1.25X median as column headers
# to run script: python deep6_compare_all_vs_all.py deep6output.csv nameofcheckcsv.csv nameofresultscsv.csv

arg1 = sys.argv[1]
arg2 = sys.argv[2]
arg3 = sys.argv[3]

checkcsv = open(arg2, 'a')
rescsv = open(arg3, 'a')

# Define a function to get the largest value, the header of the column with the largest value, and the median value of the values in columns 3, 5, 6, 7, and 8


def get_largest_value(row):
    # Get a list of the values in columns 3, 4, 5, 6, 7, and 8 by converting each value to a float
    values = [float(row[i]) for i in range(2, 8)]
    # Find the median value of the list of values using the statistics module
    median = statistics.median(values)
    # Find the largest value of the list of values
    largest = max(values)
    # Find the header of the largest value by getting the index of the largest value in the list of values and adding 2 (the first two columns are not being compared)
    largest_column = header[values.index(largest) + 2]
    # Return the largest value, header, and median
    return largest, largest_column, median


with open(arg1) as f:
    # Create a reader object to read the csv file
    reader = csv.reader(f)
    # Get the header row and store it in the header variable
    header = next(reader)
    # Loop through each row in the reader object
    for row in reader:
        # Print the values in each row with "Row" as the first column
        print(row[0] + "\t" + "\t".join(row[2:]))
        checkcsv.write(row[0] + "\t" + "\t".join(row[2:]))
        checkcsv.write("\n")
        # Get the largest value, header of the column with the largest value, and median value for the current row
        largest_value, largest_column, median=get_largest_value(row)
        # Check if the largest value is 1.25 times higher than the median
        if largest_value > 1.25 * median:
            # Print a message saying the largest value is 1.25 times higher than the median
            print(f"{row[0]},{row[1]},({largest_value}),{largest_column},higher")
            # Print a message saying the largest value is 1.25 times higher than the median in a new csv file with the name provided as argument 3
            rescsv.write(
                f"{row[0]},{row[1]},({largest_value}),{largest_column},higher\n")
            checkcsv.write(
                f"{row[0]},{row[1]},({largest_value}),{largest_column},higher\n")
        else:
            # Print a message saying the largest value is not 1.25 times higher than the median
            print(f"{row[0]},{row[1]},({largest_value}),{largest_column},lower")
            # Print a message saying the largest value is not 1.25 times higher than the median in a new csv file with the name provided as argument 3
            rescsv.write(
                f"{row[0]},{row[1]},({largest_value}),{largest_column},lower\n")
            checkcsv.write(
                f"{row[0]},{row[1]},({largest_value}),{largest_column},lower\n")

checkcsv.close()
rescsv.close()
