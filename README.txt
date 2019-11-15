Ayca Begum Tascioglu
21600907

To run the code:

1) Please go to the directory of where source code(allalign.py) located.
 
2) Add sequence file that you want to use, in the same file with the source code.

3) After, in the terminal go to the directory where allalign.py and desired input files( i.e sequences.fasta and patterns.fa etc) are located.

4) type 'make'.

5) Program asks for the file which sequence1 and sequence2 are located in.
	Type the name of pattern file (should be located in same file with program code):

	>>sequences.fasta 

6) The program will run and waits for user input as modes 0,1,2,3,4 in a while loop.
	Type 1 for Global Alignment with naive gap scoring 
	Type 2 for Global Alignment with affine gap scoring
	Type 3 for Local Alignment with naive gap scoring 
	Type 4 for Local Alignment with affine gap scoring
	For Exit type 0

	>> 1

7) Program creates result.txt which is in same file with source code; saves alignment1 and alignment2 in it.
