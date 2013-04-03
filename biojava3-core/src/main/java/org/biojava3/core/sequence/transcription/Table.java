package org.biojava3.core.sequence.transcription;

import java.util.Arrays;
import java.util.List;

import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.Compound;
import org.biojava3.core.sequence.template.CompoundSet;

/**
 * Provides a way of separating us from the specific {@link IUPACTable} even
 * though this is the only implementing class for the interface.
 *
 * @author ayates
 */
public interface Table {

  List<Codon> getCodons(CompoundSet<NucleotideCompound> nucelotides,
      CompoundSet<AminoAcidCompound> aminoAcids);

  CompoundSet<Codon> getCodonCompoundSet(
      final CompoundSet<NucleotideCompound> rnaCompounds,
      final CompoundSet<AminoAcidCompound> aminoAcidCompounds);

  /**
   * Returns true if the given compound could have been a start amino acid;
   * this does not assert if the codon that actually coded for the amino
   * acid was a start codon. This is as accurate a call as we can make with an
   * {@link AminoAcidCompound}.
   */
  boolean isStart(AminoAcidCompound compound);

  /**
   * Instance of a Codon which is 3 {@link NucleotideCompound}s, its
   * corresponding {@link AminoAcidCompound} and if it is a start or stop codon.
   * The object implements hashCode & equals but according to the nucleotide
   * compounds & not to the designation of it being a start, stop & amino
   * acid compound
   *
   * @author ayates
   *
   */
  public static class Codon implements Compound {

    private final NucleotideCompound one;
    private final NucleotideCompound two;
    private final NucleotideCompound three;
    private final boolean            start;
    private final boolean            stop;
    private final AminoAcidCompound  aminoAcid;
    private final String             stringified;

    public Codon(NucleotideCompound one, NucleotideCompound two,
        NucleotideCompound three, AminoAcidCompound aminoAcid, boolean start,
        boolean stop) {
      this.one = one;
      this.two = two;
      this.three = three;
      this.start = start;
      this.stop = stop;
      this.aminoAcid = aminoAcid;
      stringified = one.getBase().toUpperCase()
          + two.getBase().toUpperCase()
          + three.getBase().toUpperCase();
    }

    public Codon(NucleotideCompound one, NucleotideCompound two,
        NucleotideCompound three) {
      this(one,two,three,null,false,false);
    }

    public NucleotideCompound getOne() {
      return one;
    }

    public NucleotideCompound getTwo() {
      return two;
    }

    public NucleotideCompound getThree() {
      return three;
    }

    public boolean isStart() {
      return start;
    }

    public boolean isStop() {
      return stop;
    }

    public AminoAcidCompound getAminoAcid() {
      return aminoAcid;
    }

    public List<NucleotideCompound> getAsList() {
      return Arrays.asList(one,two,three);
    }

    public boolean equalsNucelotides(NucleotideCompound... compounds) {
      return compounds[0].equalsIgnoreCase(one)
          && compounds[1].equalsIgnoreCase(two)
          && compounds[2].equalsIgnoreCase(three);
    }

    public boolean equals(Object o) {
      if(o == null) {
        return false;
      }
      if(o instanceof Codon) {
        Codon them = (Codon)o;
        return equalsNucelotides(them.getOne(), them.getTwo(), them.getThree());
      }
      return false;
    }

    @Override
    public int hashCode() {
      int result = getOne().getShortName().hashCode();
      result = 37 * result + getTwo().getShortName().hashCode();
      result = 37 * result + getThree().getShortName().hashCode();
      return result;
    }

    @Override
    public String toString() {
      return stringified;
    }


    public boolean equalsIgnoreCase(Compound compound) {
      return toString().equalsIgnoreCase(compound.toString());
    }


    public String getDescription() {
      throw new UnsupportedOperationException("Not supported");
    }


    public String getLongName() {
      throw new UnsupportedOperationException("Not supported");
    }


    public Float getMolecularWeight() {
      throw new UnsupportedOperationException("Not supported");
    }


    public String getShortName() {
      return stringified;
    }


    public void setDescription(String description) {
      throw new UnsupportedOperationException("Not supported");
    }


    public void setLongName(String longName) {
      throw new UnsupportedOperationException("Not supported");
    }


    public void setMolecularWeight(Float molecularWeight) {
      throw new UnsupportedOperationException("Not supported");
    }


    public void setShortName(String shortName) {
      throw new UnsupportedOperationException("Not supported");
    }
  }

}
