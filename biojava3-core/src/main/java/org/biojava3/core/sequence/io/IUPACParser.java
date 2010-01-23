package org.biojava3.core.sequence.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.StringTokenizer;

import org.biojava3.core.exceptions.ParserException;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.transcription.TranslationTable.Codon;

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
 * @author Andy Yates
 */
public class IUPACParser {

  private static final String IUPAC_LOCATION = "/org/biojava3/core/sequence/iupac.txt";

  private InputStream is;
  private boolean cleanup;
  private List<Table> tables;

  public IUPACParser() {
    is = getClass().getClassLoader().getResourceAsStream(IUPAC_LOCATION);
    cleanup = true;
  }

  public IUPACParser(InputStream is) {
    this.is = is;
    cleanup = false;
  }

  public List<Table> getTables() {
   if(tables == null) {
     tables = parseTables();
   }
   return tables;
  }

//  public void populateTranslationTable(TranslationTable table,
//      CompoundSet<NucleotideCompound> nucleotideCompounds,
//      CompoundSet<AminoAcidCompound> aminoAcids) {
//    for(Table t: getTables()) {
//
//    }
//  }

  private List<Table> parseTables() {
    List<Table> tables = new ArrayList<Table>();
    try {
      BufferedReader br = new BufferedReader(new InputStreamReader(is));
      String line;
      Integer id = null;
      String name, aa, starts, baseone, basetwo, basethree;
      name = aa = starts = baseone = basetwo = basethree = null;
      while( (line = br.readLine()) != null) {
        if(line.equalsIgnoreCase("//")) {
          tables.add(new Table(name, id, aa, starts, baseone, basetwo, basethree));
          name = aa = starts = baseone = basetwo = basethree = null;
          id = null;
        }
        else {
          String[] keyValue = line.split("\\s+=\\s+");
          if(keyValue[0].equals("AAs")) {
            aa = keyValue[1];
          }
          else if(keyValue[0].equals("Starts")) {
            starts = keyValue[1];
          }
          else if(keyValue[0].equals("Base1")) {
            baseone = keyValue[1];
          }
          else if(keyValue[0].equals("Base2")) {
            basetwo = keyValue[1];
          }
          else if(keyValue[0].equals("Base3")) {
            basethree = keyValue[1];
          }
          else {
            name = keyValue[0];
            id = Integer.parseInt(keyValue[1]);
          }
        }
      }
    }
    catch(IOException e) {
      throw new ParserException("Problem whilst parsing IUPAC", e);
    }
    finally {
      if(cleanup) {
        try {
          is.close();
        } catch (IOException e) {
          //TODO Swallowed exception need to deal with
        }
      }
    }

    return tables;
  }

  /**
   * Holds the concept of a codon table from the IUPAC format
   *
   * @author Andy Yates
   */
  public static class Table {

    private final Integer id;
    private final String name;
    private final String aminoAcidString;
    private final String startCodons;
    private final String baseOne;
    private final String baseTwo;
    private final String baseThree;

    public Table(String name, int id, String aminoAcidString,
        String startCodons, String baseOne, String baseTwo, String baseThree) {
      this.aminoAcidString = aminoAcidString;
      this.startCodons = startCodons;
      this.name = name;
      this.id = id;
      this.baseOne = baseOne;
      this.baseTwo = baseTwo;
      this.baseThree = baseThree;
    }

    public Table(String name, Integer id, String aminoAcidString,
        String startCodons) {
      this(
          name, id, aminoAcidString, startCodons,
          "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
          "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
          "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
      );
    }

    public Integer getId() {
      return id;
    }

    public String getName() {
      return name;
    }

    public List<Codon> getCodons(CompoundSet<NucleotideCompound> nucelotides, CompoundSet<AminoAcidCompound> aminoAcids) {

      List<String> aminoAcidStrings   = aminoAcids();
      List<String> startCodonStrings  = startCodons();
      List<List<String>> codonStrings = codonStrings();

      List<Codon> codons = new ArrayList<Codon>(aminoAcidStrings.size());

      for(int i=0; i<aminoAcidStrings.size(); i++) {

        List<String> codonString = codonStrings.get(i);
        NucleotideCompound one   =
          nucelotides.getCompoundForString(codonString.get(0));
        NucleotideCompound two   =
          nucelotides.getCompoundForString(codonString.get(1));
        NucleotideCompound three =
          nucelotides.getCompoundForString(codonString.get(2));
        boolean isStart = (startCodonStrings.get(i) == "*");
        AminoAcidCompound aminoAcid =
          aminoAcids.getCompoundForString(aminoAcidStrings.get(i));

        codons.add(new Codon(one, two, three, isStart, aminoAcid));
      }

      return codons;
    }

    private List<List<String>> codonStrings() {
      List<List<String>> codons = new ArrayList<List<String>>();
      for(int i=0; i<baseOne.length(); i++) {
        List<String> codon = Arrays.asList(
            baseOne.substring(i,i),
            baseTwo.substring(i,i),
            baseThree.substring(i,i));
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