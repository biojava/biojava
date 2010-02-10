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

package org.biojava.bio.symbol;

/**
 * <p>Encapsulates the mapping from a source to a destination
 * alphabet.</p>
 *
 * <p>A TranslationTable is in effect a map or function with the
 * source domain being getSourceAlphabet() and the target domain being
 * getTargetAlphabet().  The method translate() maps a single symbol
 * from source to target.</p>
 *
 * <p>It is presumed that there will be some explicit declaration of
 * the mapping for attomic symbols, and that the mapping for all other
 * symbols will be infered from these.</p>
 *
 * <p>If you wish to translate every symbol in a symbol list then use
 * TranslatedSymbolList to automate the job. If you want to translate
 * windowed regions then first construct a WindowedSymbolList from the
 * original sequence and then build a TranslatedSymbolList from this
 * windowed view.</p>
 *
 * <p>There are a range of named tables available. Their names are
 * specified by the public static final fields in this interface. Each
 * table can be retrieved by calling the static method
 * RNATools.getGeneticCode(name) which will return the appropriate
 * TanslationTable instance.</p>
 *
 * @author Matthew Pocock
 * @author Keith James
 */
public interface TranslationTable {
    /**
     * Translation table name for the universal genetic code.
     */
    public static final String UNIVERSAL = "UNIVERSAL";

    /**
     * Translation table name for the vertebrate mitochondrial genetic
     * code.
     */
    public static final String VERT_MITO = "VERTEBRATE_MITOCHONDRIAL";

    /**
     * Translation table name for the yeast mitochondrial genetic
     * code.
     */
    public static final String YEAST_MITO = "YEAST_MITOCHONDRIAL";

    /**
     * Translation table name for the mold mitochondrial genetic code.
     */
    public static final String MOLD_MITO = "MOLD_MITOCHONDRIAL";

    /**
     * Translation table name for the invertebrate mitochondrial
     * genetic code.
     */
    public static final String INVERT_MITO = "INVERTEBRATE_MITOCHONDRIAL";

    /**
     * Translation table name for the ciliate nuclear genetic code.
     */
    public static final String CILIATE_NUC = "CILIATE_NUCLEAR";

    /**
     * Translation table name for the echinoderm mitochondrial genetic
     * code.
     */
    public static final String ECHIN_MITO = "ECHINODERM_MITOCHONDRIAL";

    /**
     * Translation table name for the euplotid nuclear genetic code.
     */
    public static final String EUPL_NUC = "EUPLOTID_NUCLEAR";
    
    /**
     * Translation table name for the bacterial and plant plastid genetic code.
     */
    public static final String BACTERIAL = "BACTERIAL";

    /**
     * Translation table name for the alternative yeast nuclear
     * genetic code.
     */
    public static final String ALT_YEAST_NUC = "ALTERNATIVE_YEAST_NUCLEAR";

    /**
     * Translation table name for the ascidian mitochondrial genetic
     * code.
     */
    public static final String ASCID_MITO = "ASCIDIAN_MITOCHONDRIAL";

    /**
     * Translation table name for the flatworm mitochondrial genetic
     * code.
     */
    public static final String FWORM_MITO = "FLATWORM_MITOCHONDRIAL";

    /**
     * Translation table name for the blepharisma macronuclear genetic
     * code.
     */
    public static final String BLEPH_MNUC = "BLEPHARISMA_MACRONUCLEAR";
    /**
     * Translation table name for the chlorophycean mitochondrial genetic
     * code.
     */
    public static final String CHLORO_MITO = "CHLOROPHYCEAN_MITOCHONDRIAL";
    /**
     * Translation table name for the trematode mitochondrial genetic
     * code.
     */
    public static final String TREMA_MITO = "TREMATODE_MITOCHONDRIAL";
    /**
     * Translation table name for the scenedesmus obliquus mitochondrial genetic
     * code.
     */
    public static final String SCENE_MITO = "SCENEDESMUS_MITOCHONDRIAL";

  /**
   * The alphabet of Symbols that can be translated.
   *
   * @return the source Alphabet
   */
  public Alphabet getSourceAlphabet();
  
  /**
   * The alphabet of Symbols that will be produced.
   *
   * @return the target Alphabet
   */
  public Alphabet getTargetAlphabet();
  
  /**
   * Translate a single symbol from source alphabet to the target alphabet.
   *
   * @param sym the Symbol to translate (member of source alphabet)
   * @return the translated version of sym (member of target alphabet)
   * @throws IllegalSymbolException if sym is not a member of the source
   *         alphabet
   */
  public Symbol translate(Symbol sym) throws IllegalSymbolException;
}
