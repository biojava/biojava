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
 * Created on 01-21-2010
 */
package org.biojava3.core.sequence.io;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava3.core.exceptions.ParserException;
import org.biojava3.core.sequence.compound.AminoAcidCompound;
import org.biojava3.core.sequence.compound.NucleotideCompound;
import org.biojava3.core.sequence.io.util.ClasspathResource;
import org.biojava3.core.sequence.io.util.IOUtils;
import org.biojava3.core.sequence.template.AbstractCompoundSet;
import org.biojava3.core.sequence.template.CompoundSet;
import org.biojava3.core.sequence.transcription.Table;


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
 * Taken from <a
 * href="http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c"
 * >NCBI</a> with slight modification and put into the classpath resource.
 *
 * Takes in an ID, name, amino acid string and the locations of amino acids
 * which acts as start codons in the translation table. You can give the 3 codon
 * position strings that correspond to the amino acid string or if you are using
 * the default IUPAC codes you can use the hardcoded ones which are consistent
 * amongst all <a
 * href="http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c"> codon
 * tables</a>.
 *
 * The generated {@link IUPACTable} objects do not parse the data further until
 * requested so if you do not use a translation table your only penalty is the
 * loading of the IUPAC data from the classpath.
 *
 * @author Andy Yates
 */
public class IUPACParser {

  private static class IOD {
    public static final IUPACParser INSTANCE = new IUPACParser();
  }

  public static IUPACParser getInstance() {
    return IOD.INSTANCE;
  }

  public static final String      IUPAC_LOCATION = "org/biojava3/core/sequence/iupac.txt";

  private InputStream              is;
  private List<IUPACTable>         tables;
  private Map<String, IUPACTable>  nameLookup;
  private Map<Integer, IUPACTable> idLookup;

  /**
   * Default version and uses the classpath based IUPAC table
   */
  public IUPACParser() {
    //use the preCache version to make sure we don't keep a IO handle open
    is = new ClasspathResource(IUPAC_LOCATION, true).getInputStream();
  }

  /**
   * Allows you to specify a different IUPAC table.
   */
  public IUPACParser(InputStream is) {
    this.is = is;
  }

  /**
   * Returns a list of all available IUPAC tables
   */
  public List<IUPACTable> getTables() {
    if (tables == null) {
      tables = parseTables();
    }
    return tables;
  }

  /**
   * Returns a table by its name
   */
  public IUPACTable getTable(String name) {
    populateLookups();
    return nameLookup.get(name);
  }

  /**
   * Returns a table by its identifier i.e. 1 means universal codon tables
   */
  public IUPACTable getTable(Integer id) {
    populateLookups();
    return idLookup.get(id);
  }

  private void populateLookups() {
    if(nameLookup == null) {
      nameLookup = new HashMap<String, IUPACTable>();
      idLookup = new HashMap<Integer, IUPACTable>();
      for(IUPACTable t: getTables()) {
        nameLookup.put(t.getName(), t);
        idLookup.put(t.getId(), t);
      }
    }
  }

  private List<IUPACTable> parseTables() {
    List<IUPACTable> localTables = new ArrayList<IUPACTable>();
    List<String> lines = IOUtils.getList(is);
    Integer id = null;
    String name, aa, starts, baseone, basetwo, basethree;
    name = aa = starts = baseone = basetwo = basethree = null;
    for (String line : lines) {
      if (line.equalsIgnoreCase("//")) {
        localTables.add(new IUPACTable(name, id, aa, starts, baseone, basetwo,
            basethree));
        name = aa = starts = baseone = basetwo = basethree = null;
        id = null;
      }
      else {
        String[] keyValue = line.split("\\s*=\\s*");
        if (keyValue[0].equals("AAs")) {
          aa = keyValue[1];
        }
        else if (keyValue[0].equals("Starts")) {
          starts = keyValue[1];
        }
        else if (keyValue[0].equals("Base1")) {
          baseone = keyValue[1];
        }
        else if (keyValue[0].equals("Base2")) {
          basetwo = keyValue[1];
        }
        else if (keyValue[0].equals("Base3")) {
          basethree = keyValue[1];
        }
        else {
          name = keyValue[0];
          id = Integer.parseInt(keyValue[1]);
        }
      }
    }

    return localTables;
  }

  /**
   * Holds the concept of a codon table from the IUPAC format
   *
   * @author Andy Yates
   */
  public static class IUPACTable implements Table {

    private final Integer      id;
    private final String       name;
    private final String       aminoAcidString;
    private final String       startCodons;
    private final String       baseOne;
    private final String       baseTwo;
    private final String       baseThree;

    private final List<Codon>  codons    = new ArrayList<Codon>();
    private CompoundSet<Codon> compounds = null;

    public IUPACTable(String name, int id, String aminoAcidString,
        String startCodons, String baseOne, String baseTwo, String baseThree) {
      this.aminoAcidString = aminoAcidString;
      this.startCodons = startCodons;
      this.name = name;
      this.id = id;
      this.baseOne = baseOne;
      this.baseTwo = baseTwo;
      this.baseThree = baseThree;
    }

    /**
     * Constructor which uses the basic IUPAC codon table format. Useful
     * if you need to specify your own IUPAC table with minimal
     * definitions from your side.
     */
    public IUPACTable(String name, Integer id, String aminoAcidString,
        String startCodons) {
      this(name, id, aminoAcidString, startCodons,
          "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG",
          "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG",
          "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG");
    }

    public Integer getId() {
      return id;
    }

    public String getName() {
      return name;
    }

    /**
     * Returns true if the given compound was a start codon in this
     * codon table. This will report true if the compound could ever have
     * been a start codon.
     *
     * @throws IllegalStateException Thrown if
     * {@link #getCodons(CompoundSet, CompoundSet)} was not called first.
     */
        @Override
    public boolean isStart(AminoAcidCompound compound) throws IllegalStateException {
      if(this.codons.isEmpty()) {
        throw new IllegalStateException("Codons are empty; please request getCodons() fist before asking this");
      }
      for(Codon codon: codons) {
        //Only check if the codon was a start codon and then ask if the compound was encoded by it
        if(codon.isStart()) {
          if(codon.getAminoAcid().equalsIgnoreCase(compound)) {
            return true;
          }
        }
      }
      return false;
    }

    /**
     * Returns a list of codons where the source and target compounds
     * are the same as those given by the parameters.
     *
     * @param nucleotides The nucleotide set to use when building BioJava 
     * representations of codons
     * @param aminoAcids The target amino acid compounds objects
     */
        @Override
    public List<Codon> getCodons(CompoundSet<NucleotideCompound> nucelotides,
        CompoundSet<AminoAcidCompound> aminoAcids) {

      if (this.codons.isEmpty()) {
        List<String> aminoAcidStrings = aminoAcids();
        List<String> startCodonStrings = startCodons();
        List<List<String>> codonStrings = codonStrings();

        for (int i = 0; i < aminoAcidStrings.size(); i++) {

          List<String> codonString    = codonStrings.get(i);
          NucleotideCompound one      = getCompound(codonString, 0, nucelotides);
          NucleotideCompound two      = getCompound(codonString, 1, nucelotides);
          NucleotideCompound three    = getCompound(codonString, 2, nucelotides);
          boolean start               = ("M".equals(startCodonStrings.get(i)));
          boolean stop                = ("*".equals(aminoAcidStrings.get(i)));
          AminoAcidCompound aminoAcid = aminoAcids
              .getCompoundForString(aminoAcidStrings.get(i));
          codons.add(new Codon(new CaseInsensitiveTriplet(one, two, three), aminoAcid, start, stop));
        }
      }

      return codons;
    }

    private NucleotideCompound getCompound(List<String> compounds,
        int position, CompoundSet<NucleotideCompound> nucelotides) {
      String compound = compounds.get(position);
      NucleotideCompound returnCompound = nucelotides
          .getCompoundForString(compound);
      if (returnCompound == null) {
        if ("T".equalsIgnoreCase(compound)) {
            returnCompound = nucelotides.getCompoundForString("U");
        }
        else {
          throw new ParserException("Cannot find a compound for string "
              + compound);
        }
      }
      return returnCompound;
    }

    /**
     * Returns the compound set of codons
     */
    public CompoundSet<Codon> getCodonCompoundSet(
        final CompoundSet<NucleotideCompound> rnaCompounds,
        final CompoundSet<AminoAcidCompound> aminoAcidCompounds) {
      if (compounds == null) {
        compounds = new AbstractCompoundSet<Codon>() {
          {
            for (Codon c : getCodons(rnaCompounds, aminoAcidCompounds)) {
              addCompound(c);
            }
          }
        };
      }
      return compounds;
    }

    private List<List<String>> codonStrings() {
      List<List<String>> codons = new ArrayList<List<String>>();
      for (int i = 0; i < baseOne.length(); i++) {
        List<String> codon = Arrays.asList(Character
            .toString(baseOne.charAt(i)),
            Character.toString(baseTwo.charAt(i)), Character.toString(baseThree
                .charAt(i)));
        codons.add(codon);
      }
      return codons;
    }

    private List<String> aminoAcids() {
      return split(aminoAcidString);
    }

    private List<String> startCodons() {
      return split(startCodons);
    }

    private List<String> split(String string) {
      List<String> split = new ArrayList<String>();
      for (int i = 0; i < string.length(); i++) {
        split.add(Character.toString(string.charAt(i)));
      }
      return split;
    }
  }
}