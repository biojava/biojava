/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence;

import java.util.ArrayList;
import java.util.Collection;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;

/**
 *
 * @author Scooter
 */
public class MultipleSequenceAlignment<C extends AbstractSequence,D extends Compound> {
    int alignedSequenceLength = -1;
    ArrayList<C> sequences = new ArrayList<C>();

    public void addAlignedSequence(C sequence){
        if(alignedSequenceLength == -1){
            alignedSequenceLength = sequence.getLength();
        }
        if(sequence.getLength() != alignedSequenceLength){
            throw new SequenceLengthError(sequence.getAccession() + " length = " + sequence.getLength()  +  " not equal to msa length = " + alignedSequenceLength);
        }
        sequences.add(sequence);
    }

    public boolean removeAlignedSequence(C sequence){
        return sequences.remove(sequence);
    }

    public Collection<C> getSequences(){
        return sequences;
    }

    public int getAlignedSequenceLength(){
        return alignedSequenceLength;
    }

    public Collection<D> getCompoundsAtColumn(int column){
        if(column > alignedSequenceLength){
            throw new SequenceLengthError("Column=" + column + " is greater than the sequence length=" + alignedSequenceLength);
        }

        ArrayList<D> columnList = new ArrayList<D>();

        for(C sequence : sequences){
            D compound = (D)sequence.getCompoundAt(column);
            columnList.add(compound);
        }

        return columnList;
    }


    public static void main(String[] args) {
        MultipleSequenceAlignment<ProteinSequence,AminoAcidCompound> msaProteins = new MultipleSequenceAlignment<ProteinSequence,AminoAcidCompound>();
        msaProteins.addAlignedSequence(new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJX"));
        msaProteins.addAlignedSequence(new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJX"));
        msaProteins.addAlignedSequence(new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJX"));
        msaProteins.addAlignedSequence(new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJX"));
        msaProteins.addAlignedSequence(new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJX"));
        msaProteins.addAlignedSequence(new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJX"));
        msaProteins.addAlignedSequence(new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJX"));
        msaProteins.addAlignedSequence(new ProteinSequence("ARNDCEQGHILKMFPSTWYVBZJX"));

        Collection<AminoAcidCompound> columnAA = msaProteins.getCompoundsAtColumn(3);

        System.out.println(columnAA);

        MultipleSequenceAlignment<DNASequence,NucleotideCompound> msaDNA = new MultipleSequenceAlignment<DNASequence,NucleotideCompound>();
        msaDNA.addAlignedSequence(new DNASequence("ATCGATCGATCGATCG"));
        msaDNA.addAlignedSequence(new DNASequence("ATCGATCGATCGATCG"));
        msaDNA.addAlignedSequence(new DNASequence("ATCGATCGATCGATCG"));
        msaDNA.addAlignedSequence(new DNASequence("ATCGATCGATCGATCG"));
        msaDNA.addAlignedSequence(new DNASequence("ATCGATCGATCGATCG"));
        msaDNA.addAlignedSequence(new DNASequence("ATCGATCGATCGATCG"));
        msaDNA.addAlignedSequence(new DNASequence("ATCGATCGATCGATCG"));
        msaDNA.addAlignedSequence(new DNASequence("ATCGATCGATCGATCG"));

        Collection<NucleotideCompound> columnDNA = msaDNA.getCompoundsAtColumn(3);
        System.out.println(columnDNA);
    }
}
