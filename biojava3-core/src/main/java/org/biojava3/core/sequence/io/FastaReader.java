/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.biojava3.core.sequence.io;

import org.biojava3.core.sequence.io.template.HeaderParserInterface;
import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collection;
import org.biojava3.core.sequence.ProteinSequence;
import org.biojava3.core.sequence.compound.AminoAcidCompoundSet;
import org.biojava3.core.sequence.io.template.SequenceCreatorInterface;
import org.biojava3.core.sequence.template.AbstractSequence;
import org.biojava3.core.sequence.template.Compound;

/**
 *
 * @author Scooter Willis <willishf at gmail dot com>
 */
public class FastaReader<S extends AbstractSequence> {

    SequenceCreatorInterface sequenceCreator;
    HeaderParserInterface headerParser;
    BufferedReader br;

    public FastaReader(BufferedReader br, HeaderParserInterface headerParser, SequenceCreatorInterface sequenceCreator) {
        this.headerParser = headerParser;
        this.br = br;
        this.sequenceCreator = sequenceCreator;
    }

    public Collection<S> process() throws Exception {
        ArrayList<S> sequences = new ArrayList<S>();

        String line = "";
        String header = "";
        StringBuffer sb = new StringBuffer();
        int maxSequenceLength = -1;
        do {
            line = line.trim(); // nice to have but probably not needed
            if (line.length() != 0) {
                if (line.startsWith(">")) {
                    if (sb.length() > 0) {
                        S sequence = (S)sequenceCreator.getSequence(sb.toString());
                        headerParser.parseHeader(header, sequence);
                        sequences.add(sequence);
                        if (maxSequenceLength < sb.length()) {
                            maxSequenceLength = sb.length();
                        }
                        sb = new StringBuffer(maxSequenceLength);
                    }
                    header = line.substring(1);
                } else if (line.startsWith(";")) {
                } else {
                    sb.append(line);
                }
            }
            line = br.readLine();
        } while (line != null);
        return sequences;
    }

    public static void main(String[] args) {
        try{
            FileReader fileReader = new FileReader("/Users/Scooter/mutualinformation/project/nuclear_receptor/PF00104_small.fasta");
            BufferedReader bufferedReader = new BufferedReader(fileReader);

            FastaReader<ProteinSequence> fastaReader = new FastaReader<ProteinSequence>(bufferedReader, new GenericFastaHeaderParser(), new ProteinSequenceCreator(AminoAcidCompoundSet.getAminoAcidCompoundSet()));
            Collection<ProteinSequence> proteinSequences = fastaReader.process();
            bufferedReader.close();
            fileReader.close();
            
            System.out.println(proteinSequences);
        }catch(Exception e){
            e.printStackTrace();
        }
    }

}
