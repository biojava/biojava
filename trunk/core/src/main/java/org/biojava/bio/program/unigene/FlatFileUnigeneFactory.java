package org.biojava.bio.program.unigene;

import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.IOException;
import java.net.URL;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.BioException;
import org.biojava.bio.program.indexdb.BioStore;
import org.biojava.bio.program.indexdb.BioStoreFactory;
import org.biojava.bio.program.indexdb.IndexStore;
import org.biojava.bio.program.tagvalue.Indexer;
import org.biojava.bio.program.tagvalue.Parser;
import org.biojava.bio.program.tagvalue.ParserListener;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.io.FastaFormat;
import org.biojava.bio.seq.io.SeqIOAdapter;
import org.biojava.bio.seq.io.SequenceBuilder;
import org.biojava.bio.seq.io.SequenceBuilderFactory;
import org.biojava.bio.seq.io.StreamReader;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.utils.CommitFailure;
import org.biojava.utils.ParserException;
import org.biojava.utils.io.CountedBufferedReader;
import org.biojava.utils.io.RAF;

/**
 * <p>A UnigeneFactory that will use flat-file indexing of the unigene ascii-art
 * files.</p>
 *
 * <p><em>This class is for developers and power-users.</em> Usually you will
 * not use this class directly, but rather use UnigeneTools.loadDatabase() with
 * a file URL.</p>
 *
 * <p>This will create all the index files necisary to look up records in a timely
 * manner. It requires read/write access to the unigene directory. No files
 * will be deleted during this opperation. The indexing strategy used is
 * compattible with the OBDA flat-file indexing spec and uses the package
 * org.biojava.bio.program.indexdb and parsers that are compattible with the
 * tag-value API.</p>
 *
 * @author Matthew Pocock
 */
public class FlatFileUnigeneFactory
implements UnigeneFactory {
  private static final String DATA_INDEX = "data.index";
  private static final String LIB_INFO_INDEX = "libInfo.index";
  private static final String UNIQUE_INDEX = "unique.index";
  private static final String ALL_INDEX = "all.index";

  /**
   * Accepts all URLs that are of the file protocol.
   */
  public boolean canAccept(URL unigeneLoc) {
    return unigeneLoc.getProtocol().equals("file");
  }

  public UnigeneDB loadUnigene(URL unigeneLoc)
  throws BioException {
    if(!unigeneLoc.getProtocol().equals("file")) {
      throw new BioException(
        "Can't create unigene from non-file URL: " +
        unigeneLoc
      );
    }

    File unigeneDir = new File(unigeneLoc.getPath());
    if(!unigeneDir.exists()) {
      throw new BioException("Could not locate directory: " + unigeneDir);
    }
    if(!unigeneDir.isDirectory()) {
      throw new BioException("Expecting a directory at: " + unigeneDir);
    }


    // load a pre-made unigene file set
    try {
      return new FlatFileUnigeneDB(
        new BioStore(new File(unigeneDir, DATA_INDEX), true),
        new BioStore(new File(unigeneDir, LIB_INFO_INDEX), true),
        new BioStore(new File(unigeneDir, UNIQUE_INDEX), true),
        new BioStore(new File(unigeneDir, ALL_INDEX), true)
      );
    } catch (IOException ioe) {
      throw new BioException("Could not instantiate flat file unigene db",ioe);
    }
  }

  public UnigeneDB createUnigene(URL unigeneLoc)
  throws BioException {
    if(!unigeneLoc.getProtocol().equals("file")) {
      throw new BioException(
        "Can't create unigene from non-file URL: " +
        unigeneLoc
      );
    }

    File unigeneDir = new File(unigeneLoc.getPath());
    if(!unigeneDir.exists()) {
      throw new BioException("Could not locate directory: " + unigeneDir);
    }
    if(!unigeneDir.isDirectory()) {
      throw new BioException("Expecting a directory at: " + unigeneDir);
    }

    try {
      indexAll(unigeneDir);
      indexUnique(unigeneDir);
      indexData(unigeneDir);
      indexLibInfo(unigeneDir);
    } catch (IOException ioe) {
      throw new BioException("Failed to index data",ioe);
    }

    return loadUnigene(unigeneLoc);
  }

  private void indexData(File unigeneDir)
  throws BioException, IOException {
    // create index file for all *.data files
    File dataIndexFile = new File(unigeneDir, DATA_INDEX);
    BioStoreFactory dataBSF = new BioStoreFactory();
    dataBSF.setPrimaryKey("ID");
    dataBSF.addKey("ID", 10);
    dataBSF.setStoreLocation(dataIndexFile);
    BioStore dataStore = dataBSF.createBioStore();
    File[] dataFiles = unigeneDir.listFiles(new FileFilter() {
      public boolean accept(File pathName) {
        return pathName.getName().endsWith(".data");
      }
    });
    for(int i = 0; i < dataFiles.length; i++) {
      File f = dataFiles[i];
      try {
        Indexer indexer = new Indexer(f, dataStore);
        indexer.setPrimaryKeyName("ID");
        Parser parser = new Parser();
        ParserListener pl = UnigeneTools.buildDataParser(indexer);
        while(parser.read(
          indexer.getReader(),
          pl.getParser(),
          pl.getListener()
        )) { ; }
      } catch (ParserException pe) {
        throw new BioException("Failed to parse " + f, pe);
      }
    }
    try {
      dataStore.commit();
    } catch (CommitFailure ne) {
      throw new BioException(ne);
    }
  }

  private void indexLibInfo(File unigeneDir)
  throws BioException, IOException {
    // create index for all *.lib.info files
    File liIndexFile = new File(unigeneDir, LIB_INFO_INDEX);
    BioStoreFactory liBSF = new BioStoreFactory();
    liBSF.setPrimaryKey("ID");
    liBSF.addKey("ID", 7);
    liBSF.setStoreLocation(liIndexFile);
    BioStore liStore = liBSF.createBioStore();
    File[] liFiles = unigeneDir.listFiles(new FileFilter() {
      public boolean accept(File pathName) {
        return pathName.getName().endsWith(".lib.info");
      }
    });
    for(int i = 0; i < liFiles.length; i++) {
      File f = liFiles[i];
      try {
        Indexer indexer = new Indexer(f, liStore);
        indexer.setPrimaryKeyName("ID");
        Parser parser = new Parser();
        ParserListener pl = UnigeneTools.buildLibInfoParser(indexer);
        while(parser.read(
            indexer.getReader(),
            pl.getParser(),
            pl.getListener()
        )) { ; }
      } catch (ParserException pe) {
        throw new BioException("Failed to parse " + f, pe);
      }
    }
    try {
      liStore.commit();
    } catch (CommitFailure ne) {
      throw new BioException(ne);
    }
  }

  private void indexUnique(File unigeneDir)
  throws BioException, IOException {
    File uniqueIndex = new File(unigeneDir, UNIQUE_INDEX);
    BioStoreFactory uniqueBSF = new BioStoreFactory();
    uniqueBSF.setStoreLocation(uniqueIndex);
    uniqueBSF.setPrimaryKey("ID");
    uniqueBSF.addKey("ID", 10);
    BioStore uniqueStore = uniqueBSF.createBioStore();
    File[] uniqueFiles = unigeneDir.listFiles(new FileFilter() {
      public boolean accept(File pathName) {
        return pathName.getName().endsWith(".seq.uniq");
      }
    });
    for(int i = 0; i < uniqueFiles.length; i++) {
      File f = uniqueFiles[i];
      RAF raf = new RAF(f, "r");
      FastaIndexer indexer = new FastaIndexer(
        raf,
        uniqueStore,
        Pattern.compile("#(\\S+)"),
        1
      );
      FastaFormat format = new FastaFormat();
      SymbolTokenization tok = DNATools.getDNA().getTokenization("token");
      StreamReader sreader = new StreamReader(
        indexer.getReader(),
        format,
        tok,
        indexer
      );
      while(sreader.hasNext()) {
        sreader.nextSequence();
      }
    }
    try {
      uniqueStore.commit();
    } catch (CommitFailure ne) {
      throw new BioException(ne);
    }
  }

  private void indexAll(File unigeneDir)
  throws BioException, IOException {
    File allIndex = new File(unigeneDir, ALL_INDEX);
    BioStoreFactory allBSF = new BioStoreFactory();
    allBSF.setStoreLocation(allIndex);
    allBSF.setPrimaryKey("ID");
    allBSF.addKey("ID", 10);
    BioStore allStore = allBSF.createBioStore();
    File[] allFiles = unigeneDir.listFiles(new FileFilter() {
      public boolean accept(File pathName) {
        return pathName.getName().endsWith(".seq.all");
      }
    });
    Pattern pattern = Pattern.compile("/gb=(\\S+)");
    for(int i = 0; i < allFiles.length; i++) {
      File f = allFiles[i];
      RAF raf = new RAF(f, "r");
      CountedBufferedReader reader = new CountedBufferedReader(new FileReader(f));

      long offset = -1;
      String id = null;
      for(String line = reader.readLine(); line != null; line = reader.readLine()) {
        if(line.startsWith("#")) {
          long nof = reader.getFilePointer();
          if(id != null) {
            allStore.writeRecord(raf, offset, (int) (nof - offset), id, Collections.EMPTY_MAP);
          }
          Matcher matcher = pattern.matcher(line);
          matcher.find();
          id = matcher.group(1);
          offset = nof;
        }
      }
    }
    try {
      allStore.commit();
    } catch (CommitFailure cf) {
      throw new BioException(cf);
    }
  }

  private static class FastaIndexer implements SequenceBuilderFactory {
    private final Map map = new HashMap();
    private final RAF raf;
    private final IndexStore store;
    private final CountedBufferedReader reader;
    private final Pattern idPattern;
    private final int idGroup;

    public FastaIndexer(RAF raf, IndexStore store, Pattern idPattern, int idGroup)
    throws IOException {
      this.raf = raf;
      this.store = store;
      this.idPattern = idPattern;
      this.idGroup = idGroup;
      reader = new CountedBufferedReader(
        new FileReader(
          raf.getFile()
        )
      );
    }

    public CountedBufferedReader getReader() {
      return reader;
    }

    public SequenceBuilder makeSequenceBuilder() {
      return new SeqIOIndexer();
    }

    class SeqIOIndexer extends SeqIOAdapter implements SequenceBuilder {
      long offset = 0L;
      String id;

      public void startSequence() {
        id = null;
        offset = reader.getFilePointer();
      }

      public void addSequenceProperty(Object key, Object value) {
        if(key.equals(FastaFormat.PROPERTY_DESCRIPTIONLINE)) {
          String line = (String) value;
          Matcher m = idPattern.matcher(line);
          m.find();
          id = m.group(idGroup);
        }
      }

      public void endSequence() {
        long nof = reader.getFilePointer();
        store.writeRecord(raf, offset, (int) (nof - offset), id, map);
        offset = nof;
      }

      public Sequence makeSequence() {
        return null;
      }
    }
  }
}
