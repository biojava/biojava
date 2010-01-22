package org.biojava3.core.sequence;

import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.CompoundSet;

public class TranslationTable {


  public static class Codon {

    private final NucleotideCompound baseOne;
    private final NucleotideCompound baseTwo;
    private final NucleotideCompound baseThree;
    private boolean isStartCodon;
    private final AminoAcidCompound aminoAcid;

    public Codon(NucleotideCompound baseOne, NucleotideCompound baseTwo,
        NucleotideCompound baseThree, boolean isStartCodon, AminoAcidCompound aminoAcid) {
      super();
      this.baseOne = baseOne;
      this.baseTwo = baseTwo;
      this.baseThree = baseThree;
      this.isStartCodon = isStartCodon;
      this.aminoAcid = aminoAcid;
    }

    public NucleotideCompound getBaseOne() {
      return baseOne;
    }

    public NucleotideCompound getBaseTwo() {
      return baseTwo;
    }

    public NucleotideCompound getBaseThree() {
      return baseThree;
    }

    public boolean isStartCodon() {
      return isStartCodon;
    }

    public AminoAcidCompound getAminoAcid() {
      return aminoAcid;
    }
  }

  /**
   * Available translations
   *
   * <ul>
   * <li>1 - UNIVERSAL</li>
   * <li>2 - VERTEBRATE_MITOCHONDRIAL</li>
   * <li>3 - YEAST_MITOCHONDRIAL</li>
   * <li>4 - MOLD_MITOCHONDRIAL</li>
   * <li>5 - INVERTEBRATE_MITOCHONDRIAL</li>
   * <li>6 - CILIATE_NUCLEAR</li>
   * <li>9 - ECHINODERM_MITOCHONDRIAL</li>
   * <li>10 - EUPLOTID_NUCLEAR</li>
   * <li>11 - BACTERIAL</li>
   * <li>12 - ALTERNATIVE_YEAST_NUCLEAR</li>
   * <li>13 - ASCIDIAN_MITOCHONDRIAL</li>
   * <li>14 - FLATWORM_MITOCHONDRIAL</li>
   * <li>15 - BLEPHARISMA_MACRONUCLEAR</li>
   * <li>16 - 2CHLOROPHYCEAN_MITOCHONDRIAL</li>
   * <li>21 - TREMATODE_MITOCHONDRIAL</li>
   * <li>23 - SCENEDESMUS_MITOCHONDRIAL</li>
   * </ul>
   *
   * Taken from <a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c">NCBI</a>
   *
   * Takes in an ID, name, amino acid string and the locations of amino acids
   * which acts as start codons in the translation table. You can give the
   * 3 codon position strings that correspond to the amino acid string or
   * if you are using the default IUPAC codes you can use the hardcoded ones
   * which are consistent amongst all
   * <a href="http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c">
   * codon tables</a>.
   *
   * @author ayates
   *
   */
  public static class IUPACParser {

    private final int id;
    private final String name;
    private final String aminoAcidString;
    private final String startCodons;
    private final String baseOne;
    private final String baseTwo;
    private final String baseThree;

    public IUPACParser(String name, int id, String aminoAcidString,
        String startCodons, String baseOne, String baseTwo, String baseThree) {
      this.aminoAcidString = aminoAcidString;
      this.startCodons = startCodons;
      this.name = name;
      this.id = id;
      this.baseOne = baseOne;
      this.baseTwo = baseTwo;
      this.baseThree = baseThree;
    }

    public IUPACParser(String name, int id, String aminoAcidString,
        String startCodons) {
      this(
          name, id, aminoAcidString, startCodons,
          "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
          "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
          "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
      );
    }

    public int getId() {
      return id;
    }

    public String getName() {
      return name;
    }

    public List<Codon> getCodons(CompoundSet<NucleotideCompound> nucelotides, CompoundSet<AminoAcidCompound> aminoAcids) {

      List<String> aminoAcidStrings  = aminoAcids();
      List<String> startCodonStrings = startCodons();
      List<String> codonStrings      = codonStrings();

      List<Codon> codons = new ArrayList<Codon>(aminoAcidStrings.size());

      for(int i=0; i<aminoAcidStrings.size(); i++) {

        String codonString = codonStrings.get(i);
        NucleotideCompound one   =
          nucelotides.getCompoundForString(codonString.substring(0,0));
        NucleotideCompound two   =
          nucelotides.getCompoundForString(codonString.substring(1,2));
        NucleotideCompound three =
          nucelotides.getCompoundForString(codonString.substring(2,2));
        boolean isStart =
          (startCodonStrings.get(i) == "*") ? true : false;
        AminoAcidCompound aminoAcid =
          aminoAcids.getCompoundForString(aminoAcidStrings.get(i));

        codons.add(new Codon(one, two, three, isStart, aminoAcid));
      }

      return codons;
    }

    private List<String> codonStrings() {
      List<String> codons = new ArrayList<String>();
      for(int i=0; i<baseOne.length(); i++) {
        String codon =
            baseOne.substring(i,i) +
            baseTwo.substring(i,i) +
            baseThree.substring(i,i);
        codons.add(codon);
      }
      return codons;
    }

    private List<String> aminoAcids() {
      return toList(aminoAcidString, "");
    }

    private List<String> startCodons() {
      return toList(startCodons, "");
    }

    private List<String> toList(String string, String token) {
      List<String> strings = new ArrayList<String>(string.length());
      StringTokenizer tokenizer = new StringTokenizer(string, token);
      while(tokenizer.hasMoreElements()) {
        strings.add(tokenizer.nextToken());
      }
      return strings;
    }

  }

}
