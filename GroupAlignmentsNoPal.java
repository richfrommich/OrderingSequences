/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.TreeMap;


/**
 *
 * @author rgoldstein
 */
public class GroupAlignmentsNoPal {
    
    ArrayList<Sequence> sequenceList = new ArrayList<>();
    TreeMap<Double, SequenceCluster> sortedSequenceClusterMap = new TreeMap<>();
    int leftBound;
    int rightBound;
    int segLength;
    Random random = new Random();
    
    char[] dnaIntToChar = {'A', 'C', 'G', 'T', '-'};
    HashMap<Character, Integer> dnaCharToIntHash = new HashMap<>();

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        GroupAlignmentsNoPal groupSeq = new GroupAlignmentsNoPal();
        groupSeq.run(args);
    }
    
    
    void run(String[] args) {
        init(args);
        readFiles(args[0]);
        HashMap<String, SequenceCluster> seqClusterHash = groupIdentical();
//        hmmSortAndPrint(seqClusterHash);
        groupByGapsAndPrint(seqClusterHash);
    }
    
    void init(String[] args) {
        leftBound = Integer.parseInt(args[1]);
        rightBound = Integer.parseInt(args[2]);
        segLength = rightBound - leftBound;
        int index = 0;
        for (char dna: dnaIntToChar) {
            dnaCharToIntHash.put(dna, index++);
        }
    }

    HashMap<String, SequenceCluster> groupIdentical() {
        HashMap<String, SequenceCluster> seqClusterHash = new HashMap<>();      // HashMap of clusters labelled by sequence segment
        for (Sequence sequence : sequenceList) {
            if (seqClusterHash.containsKey(sequence.seqSegmentString)) {
                seqClusterHash.get(sequence.seqSegmentString).addSequence(sequence);
            } else {
                SequenceCluster newSeqCluster = new SequenceCluster(sequence.seqSegmentString);
                newSeqCluster.addSequence(sequence);
                seqClusterHash.put(newSeqCluster.seqClusterSegmentString, newSeqCluster);
            }
        }

        for (SequenceCluster seqCluster : seqClusterHash.values()){
//            System.out.println(seqCluster.seqClusterSegmentString);
            seqCluster.nSequences = seqCluster.seqList.size();
        }
        return seqClusterHash;
    }
        
        
    void groupByGapsAndPrint(HashMap<String, SequenceCluster> seqClusterHash) {
        HashMap<String, IndelGroup> indelGroupHash = new HashMap<>();
        for (SequenceCluster seqCluster : seqClusterHash.values()) {
            String cString = new String(seqCluster.seqClusterSegmentString);
            String indelString = cString.replaceAll("[A-Z]", "x").replaceAll("-", "_");
            seqCluster.indelString = indelString;
            if (indelGroupHash.containsKey(indelString)) {
                indelGroupHash.get(indelString).seqClusterList.add(seqCluster);
//                System.out.println(indelString + "\t" + seqCluster.seqClusterSegmentString);
            } else {
                IndelGroup newIndelGroup = new IndelGroup(indelString);
                newIndelGroup.seqClusterList.add(seqCluster);
                indelGroupHash.put(indelString,newIndelGroup);
            }   
        }
        TreeMap<Double, IndelGroup> sortIndelGroup = new TreeMap<>();
        for (IndelGroup indelGroup : indelGroupHash.values()) {
            indelGroup.computeNSequences();
//            System.out.println(indelGroup.indelString + "\t" + indelGroup.nSequences);
            sortIndelGroup.put((indelGroup.nSequences + 1.0E-10 * random.nextDouble()), indelGroup);
        }
        
        for (double dub : sortIndelGroup.descendingKeySet()) {
            IndelGroup indelGroup = sortIndelGroup.get(dub);
//            System.out.println(indelGroup.indelString + "\t" + indelGroup.nSequences);     
            HashMap<String, SequenceCluster> indelClassSeqClusterHash = new HashMap<>();
            for (SequenceCluster seqCluster : indelGroup.seqClusterList) {
                indelClassSeqClusterHash.put(seqCluster.seqClusterSegmentString, seqCluster);
//                System.out.println("\t" + seqCluster.seqClusterSegmentString);
            }
//            System.out.println("Size of indelClassSeqClusterHash " + indelClassSeqClusterHash.size() );
            hmmSortAndPrint(indelClassSeqClusterHash);
//            if (random.nextDouble() < 0.01) {
//                System.exit(1);
//            }
        }
        
    }
        
    void hmmSortAndPrint(HashMap<String, SequenceCluster> seqClusterHash) {    
        double[][] hmm = new double[segLength][6];
        for (SequenceCluster seqCluster : seqClusterHash.values()){
            for (int i = 0; i < segLength; i++) {
                hmm[i][seqCluster.seqClusterSegmentInts[i]] += seqCluster.nSequences;
            }
        }

        for (int i = 0; i < segLength; i++) {   
            double tot = hmm[i][0] + hmm[i][1] + hmm[i][2] + hmm[i][3] + hmm[i][4] + hmm[i][5];
            for (int  iBase = 0; iBase < 6; iBase++) {
                hmm[i][iBase] = Math.log((hmm[i][iBase] + 1.0E-5)/ tot);
            }
//            System.out.println(i + "\t" + Arrays.toString(hmm[i]));
        }   
        sortedSequenceClusterMap = new TreeMap<>();
        for (SequenceCluster seqCluster : seqClusterHash.values()){
            double prob = 1.0E-8 * random.nextDouble();
            for (int i = 0; i < segLength; i++) { 
                prob += hmm[i][seqCluster.seqClusterSegmentInts[i]];
            }
            seqCluster.matchToHmm = prob;
            sortedSequenceClusterMap.put(prob, seqCluster);
        }
        for (double dub  : sortedSequenceClusterMap.descendingKeySet()) {
//            System.out.format("xxx\t%.4f\t%s\t%s\n", dub, new String(sortedSequenceClusterMap.get(dub).seqClusterSegmentChars),
//                    sortedSequenceClusterMap.get(dub).indelString);
            SequenceCluster seqCluster = sortedSequenceClusterMap.get(dub);
            seqCluster.printSequences();
        }     
    }
    
    
    void readFiles(String fileName) {
//        System.out.println("Starting input");
        try {
            FileReader file = new FileReader(fileName);
            BufferedReader buff = new BufferedReader(file);
            boolean eof = false;

            String line;            
            String currentName = "";        // Holder for name
            String readSequence = "";       // Holder for sequence

            while (!eof) {
                line = buff.readLine();
                if (line == null) {         // No more sequences
                    eof = true;
                    // Store last read sequence
                    Sequence currentSeq = new Sequence(currentName, readSequence);
                    sequenceList.add(currentSeq);
                } else {
                    if (line.startsWith(">")){                                          // End of previous sequence, start of next                                                    // End of current sequence - process and prepare for next
                            if (currentName.length() > 0) {
                                Sequence currentSeq = new Sequence(currentName, readSequence);
                                sequenceList.add(currentSeq);
                            }
                            currentName = line.replaceFirst(">","");
                            readSequence = "";
                    } else {    // Current line is extension of sequence from previous line
                        readSequence = readSequence + line;
                    } 
                }
            }
            buff.close();
            file.close();
        } catch (IOException ioe) {
            System.out.println("Error -- " + ioe.toString());
            System.exit(1);
        }

//        System.out.println("Finished input");
//        System.exit(1);
    }

    
    class IndelGroup {
        int nSequences;
        String indelString;
        ArrayList<SequenceCluster> seqClusterList = new ArrayList<>();
        
        
        IndelGroup(String indelString) {
            this.indelString = indelString;
        }
        
        void computeNSequences() {
            nSequences = 0;
            for (SequenceCluster seqCluster : seqClusterList) {
                nSequences += seqCluster.nSequences;
            }
        }
        
    }
    
    
    
    
    
    class SequenceCluster {
        
        ArrayList<Sequence> seqList = new ArrayList<>();
        String seqClusterSegmentString;
        String indelString;
        char[] seqClusterSegmentChars;
        int[] seqClusterSegmentInts;
        double matchToHmm = 0;
        int nSequences;
        
        SequenceCluster(String seqClusterSegmentString) {
//            System.out.println("New Sequence Cluster " + seqClusterSegmentString);
            this.seqClusterSegmentString = seqClusterSegmentString;
            this.seqClusterSegmentChars = seqClusterSegmentString.toCharArray();
            this.seqClusterSegmentInts = new int[seqClusterSegmentChars.length];
            for (int i = 0; i < seqClusterSegmentChars.length; i++) {
                seqClusterSegmentInts[i] = 5;
                if (dnaCharToIntHash.containsKey(seqClusterSegmentChars[i])) {
//                    System.out.println(i + "\t" +  seqClusterSegmentChars[i]);
                    seqClusterSegmentInts[i] = dnaCharToIntHash.get((seqClusterSegmentChars[i]));
                }
            }
//            System.out.println("\t" + Arrays.toString(seqClusterSegmentChars));
        }
        
        void addSequence(Sequence newSeq) {
            seqList.add(newSeq);
        }
        
        void printSequences() {
            for (Sequence seq : seqList) {
                System.out.println("> " + seq.seqName);
                System.out.println(seq.seq);
            }
        }
        
    }
    

    class Sequence{
        String seqName;
        String seq;
        String seqSegmentString;
        char[] seqSegmentChars;
        
        Sequence(String seqName, String seq){
            this.seqName = seqName;
            this.seq = seq;
            this.seqSegmentString = seq.substring(leftBound, rightBound);
            this.seqSegmentChars = this.seqSegmentString.toCharArray();
//            System.out.println(this.seq);
        }
    }
    
    
    
    
}
