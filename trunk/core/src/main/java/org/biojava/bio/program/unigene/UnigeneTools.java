package org.biojava.bio.program.unigene;

import java.io.IOException;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.AnnotationType;
import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.CardinalityConstraint;
import org.biojava.bio.PropertyConstraint;
import org.biojava.bio.program.tagvalue.ChangeTable;
import org.biojava.bio.program.tagvalue.LineSplitParser;
import org.biojava.bio.program.tagvalue.ParserListener;
import org.biojava.bio.program.tagvalue.RegexParser;
import org.biojava.bio.program.tagvalue.RegexSplitter;
import org.biojava.bio.program.tagvalue.SimpleTagValueWrapper;
import org.biojava.bio.program.tagvalue.TagDelegator;
import org.biojava.bio.program.tagvalue.TagValueContext;
import org.biojava.bio.program.tagvalue.TagValueListener;
import org.biojava.bio.program.tagvalue.TagValueParser;
import org.biojava.bio.program.tagvalue.ValueChanger;
import org.biojava.utils.ParserException;

/**
 * <p>Usefull tools for working with Unigene.</p>
 *
 * <p>This class is the main port-of-call for users of the Unigene package. It
 * provides the core APIs for finding a Unigene database as well as registering
 * your own Unigene drivers. Additionaly, it contains methods to return parsers
 * for each of the main Unigene flat-file types. If you wish to bypass the
 * biojava object model entirely, you can choose to use these parsers instead.
 * </p>
 *
 * <h2>Example use</h2>
 *
 * <p>Creating a Unigene instance from your local Unigene directory (assuming
 * that you have read/write privileges to the directory)</p>
 *
 * <pre>
 * UnigeneDB unigene = UnigeneTools.createUnigene(
 *   new URL("file:///usr/local/biodata/unigene") );
 * </pre>
 *
 * <p>Fetch a unigene cluster</p>
 *
 * <pre>
 * UnigeneDB unigene = UnigeneTools.loadUnigene(
 *   new URL("file:///usr/local/biodata/unigene") );
 * UnigeneCluster cluster = unigenge.getCluster("Aga001");
 * System.out.println("Title: " + cluster.getTitle());
 * </pre>
 *
 * <p>Parse a data file yourself</p>
 *
 * <pre>
 * BufferedReader br = new BufferedReader(new FileReader(unigeneFile));
 * Parser = new Parser();
 * TagValueListener echo = new Echo();
 * ParserListener pl = UnigeneTools.buildDataParser(echo);
 *
 * while(parser.read(br, pl.getParser(), pl.getListener())) {
 *   // read an entry
 * }
 * </pre>
 *
 * @author Matthew Pocock
 */
public class UnigeneTools {
  /**
   * <p>
   * Annotation schema for all UnigeneCluster instances. This states what
   * propperties can be expected to be associated with a cluster and how many
   * values they may have.
   * </p>
   */
  public static final AnnotationType UNIGENE_ANNOTATION;
  
  /**
   * <p>
   * Annotation schema for all Unigene libraries. This states what propperties
   * can be expected to be associated with a library and how many values they
   * may have.
   * </p>
   */
  public static final AnnotationType LIBRARY_ANNOTATION;

  private static final List factories;
  private static final Map shortName2SpeciesName;
  
  static {
    factories = new ArrayList();
    registerFactory(new FlatFileUnigeneFactory());

    shortName2SpeciesName = new HashMap();
    
    shortName2SpeciesName.put("Aga", "Anophelese gambiae");
    shortName2SpeciesName.put("Hs", "Homo sapiens");
    shortName2SpeciesName.put("Aga", "Anopheles gambiae");
    shortName2SpeciesName.put("Bt", "Bos taurus");
    shortName2SpeciesName.put("Dm", "Drosophila melanogaster");
    shortName2SpeciesName.put("Dr", "Danio rario");
    shortName2SpeciesName.put("Mm", "Mus musculus");
    shortName2SpeciesName.put("Rn", "Rattus norvegicus");
    shortName2SpeciesName.put("Xl", "Xenopus laevis");
    shortName2SpeciesName.put("At", "Arabidopsis thaliana");
    shortName2SpeciesName.put("Gma", "Glycine max");
    shortName2SpeciesName.put("Hv", "Hordeum vulgare");
    shortName2SpeciesName.put("Les", "Lycopersicon esculentum");
    shortName2SpeciesName.put("Mtr", "Medicago truncatula");
    shortName2SpeciesName.put("Os", "Oryza sativa");
    shortName2SpeciesName.put("Ta", "Triticum aestivum");
    shortName2SpeciesName.put("Zm", "Zea mays");
    
    // start to build this annotation type for .data files & UnigeneCluster
    // annotation bundles
    PropertyConstraint pc_string = new PropertyConstraint.ByClass(String.class);
    PropertyConstraint pc_int = new PropertyConstraint.ByClass(Integer.class);
    
    AnnotationType.Impl at_sts = new AnnotationType.Impl();
    at_sts.setConstraints("NAME",   pc_string, CardinalityConstraint.ONE);
    at_sts.setConstraints("ACC",    pc_string, CardinalityConstraint.ZERO_OR_ONE);
    at_sts.setConstraints("DSEG",   pc_string, CardinalityConstraint.ZERO_OR_ONE);
    at_sts.setConstraints("UNISTS", pc_string, CardinalityConstraint.ONE);
    PropertyConstraint pc_sts = new PropertyConstraint.ByAnnotationType(at_sts);
    
    AnnotationType.Impl at_txmap = new AnnotationType.Impl();
    at_txmap.setConstraints("MARKER", pc_string, CardinalityConstraint.ONE);
    at_txmap.setConstraints("RHPANEL", pc_string, CardinalityConstraint.ONE);
    PropertyConstraint pc_txmap = new PropertyConstraint.ByAnnotationType(at_txmap);
    
    AnnotationType.Impl at_protsim = new AnnotationType.Impl();
    at_protsim.setConstraints("ORG", pc_string, CardinalityConstraint.ONE);
    at_protsim.setConstraints("PROTGI", pc_string, CardinalityConstraint.ONE);
    at_protsim.setConstraints("PROTID", pc_string, CardinalityConstraint.ONE);
    at_protsim.setConstraints("PCT", pc_string, CardinalityConstraint.ONE);
    at_protsim.setConstraints("ALN", pc_int, CardinalityConstraint.ONE);
    PropertyConstraint pc_prosim = new PropertyConstraint.ByAnnotationType(at_protsim);
    
    AnnotationType.Impl at_sequence = new AnnotationType.Impl();
    at_sequence.setConstraints("ACC", pc_string, CardinalityConstraint.ONE);
    at_sequence.setConstraints("NID", pc_string, CardinalityConstraint.ONE);
    at_sequence.setConstraints("PID", pc_string, CardinalityConstraint.ZERO_OR_ONE);
    at_sequence.setConstraints("CLONE", pc_string, CardinalityConstraint.ZERO_OR_ONE);
    at_sequence.setConstraints("END", pc_string, CardinalityConstraint.ZERO_OR_ONE);
    at_sequence.setConstraints("LID", pc_string, CardinalityConstraint.ZERO_OR_ONE);
    at_sequence.setConstraints("MGC", pc_string, CardinalityConstraint.ZERO_OR_ONE);
    PropertyConstraint pc_sequence = new PropertyConstraint.ByAnnotationType(at_sequence);
    
    AnnotationType.Impl unigene = new AnnotationType.Impl();
    unigene.setConstraints("ID", pc_string, CardinalityConstraint.ONE);
    unigene.setConstraints("TITLE", pc_string, CardinalityConstraint.ONE);
    unigene.setConstraints("GENE", pc_string, CardinalityConstraint.ONE);
    unigene.setConstraints("CYTOBAND", pc_string, CardinalityConstraint.ONE);
    unigene.setConstraints("EXPRESS", pc_string, CardinalityConstraint.ONE);
    unigene.setConstraints(
      "GNM_TERMINUS",
      new PropertyConstraint.Enumeration(new Object[] { "T", "I", "S" } ),
      CardinalityConstraint.ONE);
    unigene.setConstraints("LOCUSLINK", pc_string, CardinalityConstraint.ONE);
    unigene.setConstraints("CHROMOSOME", pc_string, CardinalityConstraint.ONE);
    unigene.setConstraints("STS", pc_sts, CardinalityConstraint.ANY);
    unigene.setConstraints("TXMAP", pc_txmap, CardinalityConstraint.ANY);
    unigene.setConstraints("PROSIM", pc_prosim, CardinalityConstraint.ANY);
    unigene.setConstraints("SCOUNT", pc_int, CardinalityConstraint.ONE);
    unigene.setConstraints("SEQUENCE", pc_sequence, CardinalityConstraint.ANY);
    
    UNIGENE_ANNOTATION = unigene;
    
    AnnotationType.Impl library = new AnnotationType.Impl();
    library.setConstraints("ID", pc_string, CardinalityConstraint.ONE);
    library.setConstraints("TITLE", pc_string, CardinalityConstraint.ONE);
    library.setConstraints("TISSUE", pc_string, CardinalityConstraint.ONE);
    library.setConstraints("VECTOR", pc_string, CardinalityConstraint.ONE);
    
    LIBRARY_ANNOTATION = library;
  }
  
  /**
   * Converts short species names (like Hs) to long species names (like Homo
   * Sapiens).
   *
   * @param name  the short name
   * @return the long name
   */
  public static String getSpeciesForShortName(String name) {
    return (String) shortName2SpeciesName.get(name);
  }
  
  /**
   * Generate a tag-value parser for unigene data files that will pass all
   * parsing events on to your listener.
   *
   * @param listener the TagValueListener to pass events onto
   * @return a ParserListener that is ready to consume unigene data documents
   */
  public static ParserListener buildDataParser(TagValueListener listener)
  throws ParserException {
    try {
      LineSplitParser entryParser = (LineSplitParser) LineSplitParser.GENBANK.clone();
      entryParser.setTrimValue(true);
      entryParser.setEndOfRecord("//");
      
      ChangeTable changeT = new ChangeTable();
      changeT.setSplitter(
        "EXPRESS",
        new RegexSplitter(Pattern.compile("([^;]+)"), 1)
      );
      changeT.setChanger("ALN", ChangeTable.STRING_TO_INT);
      changeT.setChanger("SCOUNT", ChangeTable.STRING_TO_INT);
      ValueChanger changer = new ValueChanger(listener, changeT);
      
      SplitAndProp splitAndProp = new SplitAndProp(
        listener,
        Pattern.compile("(\\S+?)=([^;\\s]*)")
      );
      TagDelegator entryListener = new TagDelegator(changer);
      entryListener.setListener("STS", splitAndProp);
      entryListener.setListener("PROTSIM", splitAndProp);
      entryListener.setListener("SEQUENCE", splitAndProp);
      entryListener.setListener("TXMAP", new HandleMapInterval(listener));
      
      return new ParserListener(entryParser, entryListener);
    } catch (CloneNotSupportedException cnse) {
      throw new BioError(cnse);
    }
  }
  
  /**
   * Generate a tag-value parser for the library info unigene files.
   *
   * @param listener the TagValueListener to pass events onto
   * @return a ParserListener that is ready to consume unigene lib.info files
   */
  public static ParserListener buildLibInfoParser(TagValueListener listener)
  throws IOException, ParserException{
    RegexParser parser = new RegexParser();
    parser.setContinueOnEmptyTag(false);
    parser.setEndOfRecord(TagValueParser.EMPTY_LINE_EOR);
    parser.setMergeSameTag(false);
    parser.setPattern(Pattern.compile("([^=]+)=(.*)"));
    parser.setTagGroup(1);
    parser.setValueGroup(2);
    
    return new ParserListener(parser, listener);
  }  
  
  private static class SplitAndProp
  extends SimpleTagValueWrapper {
    private Pattern splitPattern;
    
    public SplitAndProp(TagValueListener delegate, Pattern splitPattern) {
      super(delegate);
      this.splitPattern = splitPattern;
    }
    
    public void value(TagValueContext tvc, Object value)
    throws ParserException {
      TagValueListener delegate = super.getDelegate();
      
      delegate.startRecord();
      
      String sv = (String) value;
      Matcher m = splitPattern.matcher(sv);
      while(m.find()) {
        String k = m.group(1);
        String v = m.group(2);
        
        delegate.startTag(k);
        delegate.value(tvc, v);
        delegate.endTag();
      }
      
      delegate.endRecord();
    }
  }
  
  private static class HandleMapInterval
  extends SimpleTagValueWrapper {
    private Pattern pattern;
    public HandleMapInterval(TagValueListener tvl) {
      super(tvl);
      pattern = Pattern.compile("([^-]+-[^;]+);\\s+\\w+=([^;]+);\\s+\\w+=(\\S+)");
    }
    
    public void value(TagValueContext tvc, Object value)
    throws ParserException {
      TagValueListener delegate = super.getDelegate();
      
      delegate.startRecord();
      
      String sv = (String) value;
      Matcher m = pattern.matcher(sv);
      if(!m.find()) {
        throw new ParserException("Could not parse line: " + sv);
      }
      
      delegate.startTag("INTERVAL");
      delegate.value(tvc, m.group(1));
      delegate.endTag();
      
      delegate.startTag("MARKER");
      delegate.value(tvc, m.group(2));
      delegate.endTag();

      delegate.startTag("RHPANEL");
      delegate.value(tvc, m.group(3));
      delegate.endTag();
      
      delegate.endRecord();
    }
  }

  /**
   * <p>Register a UnigeneFactory.</p>
   *
   * <p>This method is for developers who have written their own UnigeneFactory
   * implementations. By default, jdbc and file URLs are handled by built-in
   * factories.</p>
   *
   * <p>When you register a factory, it will be used for all URLs that is can
   * accept. If a factory is registered afterwards that can accept the same URL,
   * the first factory registered will be used.</p>
   *
   * @param factory  the UnigeneFactory to register
   */
  public static void registerFactory(UnigeneFactory factory) {
    factories.add(factory);
  }

  /**
   * <p>Register a UnigeneFactory.</p>
   *
   * <p>This method is for developers who wish to unregister a factory.</p>
   *
   * @param factory  the UnigeneFactory to unregister
   */
  public static void unregisterFactory(UnigeneFactory factory) {
    factories.remove(factory);
  }

  /**
   * Load a UnigeneDB instance referred to by a URL.
   *
   * @param dbURL the URL location the database
   * @return a UnigeneDB instance
   * @throws BioException if there was no UnigeneFactory able to process that
   *         URL or if there was some error connecting
   */
  public static UnigeneDB loadUnigene(URL dbURL)
  throws BioException {
    return findFactory(dbURL).loadUnigene(dbURL);
  }

  /**
   * Create a new UnigeneDB instance referred to by a URL.
   *
   * @param dbURL the URL location the database
   * @return a UnigeneDB instance
   * @throws BioException if there was no UnigeneFactory able to process that
   *         URL or if there was some error creating it
   */
  public static UnigeneDB createUnigene(URL dbURL)
  throws BioException {
    return findFactory(dbURL).createUnigene(dbURL);
  }

  /**
   * Find the UnigeneFactory that can accept a URL.
   *
   * <p><em>This method is for developers only.</em> The normal way to interact
   * with factories is to call UnigeneTools.loadUnigene() and
   * UnigeneTools.createUnigene()</p>
   *
   * @param dbURL  the URL to find a factory for
   * @return the UnigeneFactory that accepts that URL
   * @throws BioException if there is no factory for that type of URL
   */
  public static UnigeneFactory findFactory(URL dbURL)
  throws BioException {
    for(Iterator i = factories.iterator(); i.hasNext(); ) {
      UnigeneFactory factory = (UnigeneFactory) i.next();
      if(factory.canAccept(dbURL)) {
        return factory;
      }
    }

    throw new BioException("No factory for unigene url: " + dbURL);
  }
}
