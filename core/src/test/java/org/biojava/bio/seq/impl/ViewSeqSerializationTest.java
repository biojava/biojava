/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */

package org.biojava.bio.seq.impl;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Iterator;

import junit.framework.TestCase;

import org.biojava.bio.Annotation;
import org.biojava.bio.seq.CircularView;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Feature;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.RangeLocation;
import org.biojava.bio.symbol.SymbolList;

/**
 * Tests for Serialization
 *
 * @author Mark Schreiber
 * @since 1.3
 */

public class ViewSeqSerializationTest extends TestCase
{
    protected Sequence seq;
    protected Sequence seq2;


    public ViewSeqSerializationTest(String name) {
        super(name);
    }

    protected void setUp() throws Exception {
        seq = new SimpleSequence(DNATools.createDNA("aacgtaggttccatgc"),
                                       "fragment1",
                                       "fragment1",
                                       Annotation.EMPTY_ANNOTATION);

        Feature.Template sft = new Feature.Template();
        sft.type = "test";
        sft.source = "test";
        sft.annotation = Annotation.EMPTY_ANNOTATION;
        sft.location = new RangeLocation(1, 3);
        seq.createFeature(sft);

        sft.location = new RangeLocation(10, 12);
        seq.createFeature(sft);
        seq = new CircularView(seq);
        sft.location = new RangeLocation(5,8);
        seq.createFeature(sft);

        ByteArrayOutputStream os = new ByteArrayOutputStream();
        ObjectOutputStream oos = new ObjectOutputStream(os);
        oos.writeObject(seq);
        oos.flush();
        oos.close();
        ObjectInputStream ois = new ObjectInputStream(
                new ByteArrayInputStream(os.toByteArray()));
        seq2 = (Sequence)ois.readObject();
    }

    private boolean compareSymbolList(SymbolList sl1, SymbolList sl2) {
        if (sl1.length() != sl2.length()) {
            return false;
        }

        Iterator si1 = sl1.iterator();
        Iterator si2 = sl2.iterator();
        while (si1.hasNext()) {
            if (! (si1.next() == si2.next())) {
                return false;
            }
        }

        return true;
    }

    public void testSymbols()
        throws Exception
    {
        assertTrue(compareSymbolList(seq,seq2));
    }

    public void testAlphabet()throws Exception{
        assertTrue(seq.getAlphabet()==seq2.getAlphabet());
    }

    public void testFeatures() throws Exception{
        assertTrue(seq.countFeatures() == seq2.countFeatures());
    }
}
