# test files

The first test file is an example of if there was PCR duplication. So our output would be one line

the second test file is where all of them are unique and they should all be printed out to the new SAM file

The thirds test file has a change in position. There is a soft clip which has to be accounted for which would make it the same PCR as the other strand
so this duplication should be accounted for and the second line should not be recorded

The fourth test file is reverse strand as well as has a soft clipping that needed to be accounted for. This will make it so that there should be 
one of the lines saved.

