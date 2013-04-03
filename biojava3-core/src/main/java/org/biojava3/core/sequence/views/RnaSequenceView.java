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
 * Created on 10-18-2010
 *
 * @author Andy Yates
 */
package org.biojava3.core.sequence.views;

import java.util.HashMap;
import java.util.Map;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.compound.RNACompoundSet;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.template.ProxySequenceReader;
import org.biojava3.core.sequence.template.SequenceProxyView;
import org.biojava3.core.sequence.template.Sequence;
import org.biojava3.core.sequence.template.SequenceMixin;

/**
 * Attempts to do on the fly translation of RNA by not requesting the compounds
 * until asked.
 *
 * @author ayates
 */
public class RnaSequenceView extends SequenceProxyView<NucleotideCompound> implements ProxySequenceReader<NucleotideCompound> {

    private CompoundSet<NucleotideCompound> rnaCompounds;
    private Map<NucleotideCompound, NucleotideCompound> dnaToRna = null;
    private Map<NucleotideCompound, NucleotideCompound> rnaToDna = null;

    public RnaSequenceView(Sequence<NucleotideCompound> sourceDna) {
        this(sourceDna, RNACompoundSet.getRNACompoundSet());
    }

    public RnaSequenceView(Sequence<NucleotideCompound> sourceDna,
            CompoundSet<NucleotideCompound> rnaCompounds) {
        super(sourceDna);
        this.rnaCompounds = rnaCompounds;
    }

    @Override
    public String getSequenceAsString() {
        return SequenceMixin.toString(this);
    }

    @Override
    public NucleotideCompound getCompoundAt(int position) {
        NucleotideCompound dna = getViewedSequence().getCompoundAt(position);
        return getDnaToRna().get(dna);
    }

    @Override
    public int getIndexOf(NucleotideCompound compound) {
        return getViewedSequence().getIndexOf(getRnaToDna().get(compound));
    }

    @Override
    public int getLastIndexOf(NucleotideCompound compound) {
        return getViewedSequence().getLastIndexOf(getRnaToDna().get(compound));
    }

    public Map<NucleotideCompound, NucleotideCompound> getRnaToDna() {
        if(rnaToDna == null) {
            buildTranslators();
        }
        return rnaToDna;
    }

    public Map<NucleotideCompound, NucleotideCompound> getDnaToRna() {
        if(dnaToRna == null) {
            buildTranslators();
        }
        return dnaToRna;
    }

    protected void buildTranslators() {
        Map<NucleotideCompound, NucleotideCompound> localDnaToRna =
                new HashMap<NucleotideCompound, NucleotideCompound>();
        Map<NucleotideCompound, NucleotideCompound> localRnaToDna =
                new HashMap<NucleotideCompound, NucleotideCompound>();

        NucleotideCompound thymine =
                getViewedSequence().getCompoundSet().getCompoundForString("T");
        NucleotideCompound lowerThymine =
                getViewedSequence().getCompoundSet().getCompoundForString("t");

        for (NucleotideCompound dnaBase : getViewedSequence().getCompoundSet().getAllCompounds()) {
            NucleotideCompound equivalent;
            if (dnaBase.equals(thymine)) {
                equivalent = rnaCompounds.getCompoundForString("U");
            }
            else if (dnaBase.equals(lowerThymine)) {
                equivalent = rnaCompounds.getCompoundForString("u");
            }
            else {
                equivalent = rnaCompounds.getCompoundForString(
                    dnaBase.toString());
            }
            localDnaToRna.put(dnaBase, equivalent);
            localRnaToDna.put(equivalent, dnaBase);
        }
        this.dnaToRna = localDnaToRna;
        this.rnaToDna = localRnaToDna;
    }

    @Override
    public void setCompoundSet(CompoundSet<NucleotideCompound> compoundSet) {
        throw new UnsupportedOperationException("Unsupported operation; create a new viewed sequence");
    }

    @Override
    public void setContents(String sequence) {
        throw new UnsupportedOperationException("Unsupported operation; create a new viewed sequence");
    }
}
