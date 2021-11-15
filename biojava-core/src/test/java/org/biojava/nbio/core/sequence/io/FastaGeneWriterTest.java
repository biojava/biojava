package org.biojava.nbio.core.sequence.io;

import static org.junit.jupiter.api.Assertions.assertEquals;

import java.io.ByteArrayOutputStream;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

import org.biojava.nbio.core.sequence.AccessionID;
import org.biojava.nbio.core.sequence.ChromosomeSequence;
import org.biojava.nbio.core.sequence.GeneSequence;
import org.biojava.nbio.core.sequence.Strand;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound;
import org.junit.jupiter.api.Test;

class FastaGeneWriterTest {

    @Test
    void basicGeneWriterTest() throws Exception {

        List<GeneSequence> sequences = new ArrayList<GeneSequence>();
        ChromosomeSequence seq1 = new ChromosomeSequence(
                "ATATATATATATATATATATATATATATATATACGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCATATATATATATATATATATATACGCGCGCGCGCGCGCGCATATATATATATATATATATATATATATATATACGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCATATATATATATATATATATATACGCGCGCGCGCGCGCGC");
        GeneSequence gene1 = seq1.addGene(new AccessionID("gene1"), 1, 20, Strand.POSITIVE);

        gene1.addExon(new AccessionID("t1_1_10"), 1, 10);
        gene1.addExon(new AccessionID("t1_12_15"), 12, 15);
        GeneSequence gene2 = seq1.addGene(new AccessionID("gene2"), 1, 20, Strand.NEGATIVE);

        gene2.addExon(new AccessionID("t2_1_10"), 1, 10);
        gene2.addExon(new AccessionID("t2_12_15"), 12, 15);
        sequences.add(gene1);
        sequences.add(gene2);

        ByteArrayOutputStream os = new ByteArrayOutputStream();
        FastaGeneWriter fastaWriter = new FastaGeneWriter(os, sequences,
                new GenericFastaHeaderFormat<GeneSequence, NucleotideCompound>(), true);
        fastaWriter.process();

        String output = new String(os.toByteArray(), "UTF-8");
        String [] lines =  output.split("\n");
        assertEquals(4,lines.length);
        assertEquals(">gene1", lines[0]);
        assertEquals("ATATATATATaTATAtatat", lines[1]);
        assertEquals(">gene2", lines[2]);
        assertEquals("tatatATATaTATATATATA", lines[3]);
    }
}
