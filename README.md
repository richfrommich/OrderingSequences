# OrderingSequences
Java program for sorting sequences so that the disorder is pushed to the bottom


In order to compile the program, simply use:
  javac GroupAlignmentsNoPal.java
  
To use the program, use:
  java GroupAlignmentsNoPal [input alignment file] [start] [end] \> [output alignment file]
  
where [input alignment file] is replaced by the name of the input alignment file, 
similarly for the output alignment file, and start and end are the numbers defining the
region used for ordering the sequences.

The ">" is, however, necessary. Otherwise the output alignment will appear on the screen
