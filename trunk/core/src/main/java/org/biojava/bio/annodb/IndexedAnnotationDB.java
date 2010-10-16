package org.biojava.bio.annodb;

import java.beans.XMLDecoder;
import java.beans.XMLEncoder;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.lang.reflect.Method;
import java.lang.reflect.Modifier;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.Annotation;
import org.biojava.bio.AnnotationTools;
import org.biojava.bio.AnnotationType;
import org.biojava.bio.BioException;
import org.biojava.bio.program.indexdb.BioStore;
import org.biojava.bio.program.indexdb.BioStoreFactory;
import org.biojava.bio.program.indexdb.Record;
import org.biojava.bio.program.tagvalue.AnnotationBuilder;
import org.biojava.bio.program.tagvalue.Index2Model;
import org.biojava.bio.program.tagvalue.Indexer2;
import org.biojava.bio.program.tagvalue.Parser;
import org.biojava.bio.program.tagvalue.ParserListener;
import org.biojava.bio.program.tagvalue.TagValueListener;
import org.biojava.bio.seq.io.filterxml.XMLAnnotationTypeHandler;
import org.biojava.bio.seq.io.filterxml.XMLAnnotationTypeWriter;
import org.biojava.utils.AssertionFailure;
import org.biojava.utils.CommitFailure;
import org.biojava.utils.ParserException;
import org.biojava.utils.io.RandomAccessReader;
import org.biojava.utils.stax.SAX2StAXAdaptor;
import org.biojava.utils.xml.PrettyXMLWriter;
import org.biojava.utils.xml.XMLWriter;
import org.xml.sax.SAXException;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.XMLReaderFactory;

/**
 * <p>A database of Annotation instances backed by an indexed file set.</p>
 *
 * @author Matthew Pocock
 * @since 1.3
 */
public class IndexedAnnotationDB
implements AnnotationDB {
  private final BioStore store;
  private final AnnotationType schema;
  private final ParserListenerFactory plFactory;
  private final ParserListener parserListener;
  private final AnnotationBuilder annBuilder;
  private final Parser recordParser;

  /**
   * Create a new IndexedAnnotationDB.
   *
   * @param dbName
   * @param storeLoc
   * @param model
   * @param toIndex
   * @param maxKeyLen
   * @param schema
   * @param plFactory
   * @throws BioException
   * @throws CommitFailure
   * @throws IOException
   * @throws ParserException
   */
  public IndexedAnnotationDB(
    String dbName,
    File storeLoc,
    Index2Model model,
    List toIndex,
    int maxKeyLen,
    AnnotationType schema,
    ParserListenerFactory plFactory
  ) throws BioException, CommitFailure, IOException, ParserException {
    // state
    BioStoreFactory bsf = new BioStoreFactory();
    bsf.setStoreName(dbName);
    bsf.setPrimaryKey(model.getPrimaryKeyName());
    bsf.setStoreLocation(storeLoc);
    
    for(Iterator i = model.getKeys().iterator(); i.hasNext(); ) {
      String key = (String) i.next();
      bsf.addKey(key, maxKeyLen);
    }
    
    this.store = bsf.createBioStore();
    this.schema = schema;
    this.plFactory = plFactory;
    this.annBuilder = new AnnotationBuilder(schema);
    this.parserListener = plFactory.getParserListener(annBuilder);
    this.recordParser = new Parser();
    
    // persistance
    File factoryFile = new File(store.getLocation(), "ParserListenerFactory.xml");
    XMLEncoder xmlEnc = new XMLEncoder(
      new BufferedOutputStream(
        new FileOutputStream(
          factoryFile
        )
      )
    );
    xmlEnc.writeObject(plFactory);
    xmlEnc.close();
    
    File schemaFile = new File(store.getLocation(), "schema.xml");
    PrintWriter schemaPW = new PrintWriter(
      new FileWriter(
        schemaFile
      )
    );
    XMLWriter schemaWriter = new PrettyXMLWriter(schemaPW);
    XMLAnnotationTypeWriter schemaTW = new XMLAnnotationTypeWriter();
    schemaTW.writeAnnotationType(schema, schemaWriter);
    schemaPW.flush();
    schemaPW.close();
    
    for(Iterator fi = toIndex.iterator(); fi.hasNext(); ) {
      File file = (File) fi.next();
      
      Indexer2 ndx = new Indexer2(file, store, model);
      ParserListener pl = plFactory.getParserListener(ndx);
      Parser parser = new Parser();
      while(parser.read(ndx.getReader(), pl.getParser(), pl.getListener())) {
        ;
      }
    }
    
    store.commit();
  }

  /**
   * Initialise the db from a store.
   *
   * @param store         the BioStore to initalise from
   * @throws IOException  if there was an IO fault accessing the store
   * @throws SAXException if the XML configuration file is corrupted
   */
  public IndexedAnnotationDB(BioStore store) throws IOException, SAXException {
    this.store = store;
    
    File factoryFile = new File(store.getLocation(), "ParserListenerFactory.xml");
    XMLDecoder xmlDec = new XMLDecoder(
      new BufferedInputStream(
        new FileInputStream(
          factoryFile
        )
      )
    );
    this.plFactory = (ParserListenerFactory) xmlDec.readObject();
    xmlDec.close();
    
    XMLReader parser = XMLReaderFactory.createXMLReader();
    XMLAnnotationTypeHandler annTypeH = new XMLAnnotationTypeHandler();
    parser.setContentHandler(
      new SAX2StAXAdaptor(
        annTypeH
      )
    );
    this.schema = annTypeH.getAnnotationType();
    
    this.annBuilder = new AnnotationBuilder(schema);
    this.parserListener = plFactory.getParserListener(annBuilder);
    this.recordParser = new Parser();
  }
  
  public String getName() {
    return store.getName();
  }

  public AnnotationType getSchema() {
    return schema;
  }
  
  public Iterator iterator() {
    return new Iterator() {
      Iterator rli = store.getRecordList().iterator();
      
      public boolean hasNext() {
        return rli.hasNext();
      }
      
      public Object next() {
        try {
          return process((Record) rli.next());
        } catch (Exception e) {
          throw new RuntimeException(e);
        }
      }
      
      public void remove() {
        throw new UnsupportedOperationException();
      }
    };
  }
  
  public int size() {
    return store.getRecordList().size();
  }
  
  public AnnotationDB filter(AnnotationType at) {
    AnnotationType schema = AnnotationTools.intersection(at, this.schema);
    
    if(schema != AnnotationType.NONE) {
      return new LazyFilteredAnnotationDB("", this, schema);
    } else {
      return AnnotationDB.EMPTY;
    }
  }
  
  public AnnotationDB search(AnnotationType at) {
    return new LazySearchedAnnotationDB("", this, at);
  }

  /**
   * Get the ParserListenerFactory used by this IndexedAnnotationDB.
   *
   * @return the ParserListenerFactory
   */
  public ParserListenerFactory getParserListenerFactory() {
    return plFactory;
  }
  
  private Annotation process(Record rec)
  throws IOException, ParserException {
    RandomAccessReader rar = new RandomAccessReader(rec.getFile());
    rar.seek(rec.getOffset());
    BufferedReader reader = new BufferedReader(rar);
    recordParser.read(reader, parserListener.getParser(), parserListener.getListener());
    return annBuilder.getLast();
  }

  /**
   * A factory for retrieving parsers and listeners.
   *
   * @author Matthew Pocock
   * @since 1.3
   */
  public static interface ParserListenerFactory
  extends Serializable {
    /**
     * Get the ParserListener for a TagValueListener.
     *
     * @param listener the TagValueListener to process
     * @return the ParserListener for this
     */
    public ParserListener getParserListener(TagValueListener listener);
  }

  /**
   * An implementation of ParserListenerFactory that uses a static method.
   *
   * @author Matthew Pocock
   * @since 1.3
   */
  public static class StaticMethodRPFactory
  implements ParserListenerFactory {
    private final  Method method;

    /**
     * Create a new StaticMethodRPFactory for a method.
     *
     * @param method  a Method to use
     * @throws IllegalArgumentException  if the Method is not statically scoped,
     *    or does not return a ParserListener or take a single argument of type
     *    TagValueListener
     */
    public StaticMethodRPFactory(Method method)
    throws IllegalArgumentException {
      if( (method.getModifiers() & Modifier.STATIC) != Modifier.STATIC ) {
        throw new IllegalArgumentException("Method must be static");
      }
      
      if(method.getReturnType() != ParserListener.class) {
        throw new IllegalArgumentException("Method must return a ParserListener instance");
      }
      
      if(
        method.getParameterTypes().length != 1 ||
        method.getParameterTypes()[0] != TagValueListener.class
      ) {
        throw new IllegalArgumentException("Method must accept a single TagValueListener as it's sole parameter");
      }
      
      this.method = method;
    }

    /**
     * Get the Method used.
     *
     * @return  the Method used.
     */
    public Method getMethod() {
      return method;
    }
    
    public ParserListener getParserListener(TagValueListener tvl) {
      try {
        return (ParserListener) method.invoke(null, new Object[] { tvl });
      } catch (Exception e) {
        throw new AssertionFailure("Could not invoke underlying method.", e);
      }
    }
  }
}

