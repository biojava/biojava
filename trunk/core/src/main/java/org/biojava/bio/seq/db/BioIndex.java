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
package org.biojava.bio.seq.db;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.io.RandomAccessFile;
import java.util.AbstractList;
import java.util.AbstractSet;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.seq.io.SequenceBuilderFactory;
import org.biojava.bio.seq.io.SequenceFormat;
import org.biojava.bio.seq.io.SymbolTokenization;

/**
 * The original object for indexing sequence files.
 *
 * <p>This class may not be thread-safe.</p>
 *
 * @author Matthew Pocock
 * @author Thomas Down
 */
public class BioIndex implements IndexStore {
  private static Comparator STRING_CASE_SENSITIVE_ORDER = new Comparator() {
    public int compare(Object a, Object b) {
      return ((Comparable) a).compareTo(b);
    }
  };

  private File indexDirectory;

  private int fileCount;
  private File[] fileIDToFile;

  private FileAsList indxList;

  private Set idSet = new ListAsSet();

  private String name;
  private SequenceFormat format;
  private SequenceBuilderFactory sbFactory;
  private SymbolTokenization symbolTokenization;

  {
    fileCount = 0;
    fileIDToFile = new File[4];
  }

  public BioIndex(
    File indexDirectory,
    String namespace,
    int idLength
  ) throws IOException, BioException {
    if(indexDirectory.exists()) {
      throw new BioException(
        "Can't create new index as directory already exists: " +
        indexDirectory
      );
    }

    // create directory
    indexDirectory.mkdirs();

    // create BIOINDEX.dat
    {
      File bioindex = new File(indexDirectory, "BIOINDEX.dat");
      bioindex.createNewFile();
      PrintWriter pw = new PrintWriter(new FileWriter(bioindex));
      pw.println("index\tflat/1");
      pw.close();
    }

    // create fileids.dat
    PrintWriter fileidsWriter;
    {
      File fileids = new File(indexDirectory, "fileids.dat");
      fileids.createNewFile();
      fileidsWriter = new PrintWriter(
        new FileWriter(
          fileids
        )
      );
    }

    // create config.dat
    PrintWriter configWriter;
    {
      File config = new File(indexDirectory, "config.dat");
      config.createNewFile();
      configWriter = new PrintWriter(new FileWriter(config));
      configWriter.println("namespace\t" + namespace);
    }

    // create index file
    {
      String uniqueName = "key_" + namespace + ".key";
      File unique = new File(indexDirectory, uniqueName);
      unique.createNewFile();

      int recordLen =
        idLength +                                   // id
        1 +                                          // tab
        4 +                                          // 9999 files
        1 +                                          // tab
        String.valueOf(Long.MAX_VALUE).length() +    // space for any long
        1 +                                          // tab
        String.valueOf(Integer.MAX_VALUE).length() + // space for any int
        "\n".length()                                // new line (os dependant)
        ;

      indxList = new IndexFileAsList(
        new RandomAccessFile(unique, "rw"),
        recordLen
      );

      fileidsWriter.println(uniqueName + "\t" + recordLen);
    }

    // other field initialization to get things going
    fileCount = 0;
    fileIDToFile = new File[4];

    configWriter.close();
    fileidsWriter.close();
  }

  /**
   * Load an existing index file.
   *
   * If indexDirectory does not exist, or is not a bioindex stoore, this will
   * barf.
   */
  public BioIndex(
    File indexDirectory
  ) throws IOException, BioException {
    this.indexDirectory = indexDirectory;

    if(!indexDirectory.exists()) {
      throw new BioException(
        "Tried to load non-existant index: " +
        indexDirectory
      );
    }

    // read in the global config
    {
      System.out.println("Global");
      Map config = new HashMap();
      BufferedReader fi = new BufferedReader(
        new FileReader(
          new File(indexDirectory, "config.dat")
        )
      );
      for(String line = fi.readLine(); line != null; line = fi.readLine()) {
        int tab = line.indexOf("\t");
        config.put(line.substring(0, tab), line.substring(tab + 1));
      }
      String namespace = (String) config.get("namespace");
      RandomAccessFile indxFile = new RandomAccessFile("key_" + namespace + ".key", "rw");
      int recLen = guessRecLen(indxFile);
      indxList = new IndexFileAsList(indxFile, recLen);
    }

    // set up file set
    {
      System.out.println("Files");
      fileCount = 0;
      fileIDToFile = new File[4];

      BufferedReader fi = new BufferedReader(
        new FileReader(
          new File(indexDirectory, "fileids.dat")
        )
      );
      for(String line = fi.readLine(); line != null; line = fi.readLine()) {
        StringTokenizer sTok = new StringTokenizer("\t");
        int id = Integer.parseInt(sTok.nextToken());
        File file = new File(sTok.nextToken());
        long fileLength = Long.parseLong(sTok.nextToken());

        if(file.length() != fileLength) {
          throw new BioException("File length changed: " + file + " "
          + file.length() + " vs " + fileLength);
        }

        fileIDToFile[id] = file;
      }
    }
  }

  private File getFileForID(int fileId) {
    return fileIDToFile[fileId];
  }

  private int getIDForFile(File file) {
    // scan list
    for(int i = 0; i < fileCount; i++) {
      if(file.equals(fileIDToFile[i])) {
        return i;
      }
    }

    // extend fileIDToFile array
    if(fileCount >= fileIDToFile.length) {
      File[] tmp = new File[fileIDToFile.length + 4]; // 4 is magic number
      System.arraycopy(fileIDToFile, 0, tmp, 0, fileCount);
      fileIDToFile = tmp;
    }

    // add the unseen file to the list
    fileIDToFile[fileCount] = file;
    return fileCount++;
  }

  public String getName() {
    return this.name;
  }

  public int guessRecLen(RandomAccessFile file)
  throws IOException {
    file.seek(0l);
    int b = 0;
    while(b != '\n' && b != '\r') {
      b  = file.read();
    }

    int offset = (int) file.getFilePointer();

    if(b == '\n') {          // \n
      return offset + 1;
    } else {
      b = file.read();
      if(b == '\n') {        // \r\n
        return offset + 2;
      } else {               // \r
        return offset + 1;
      }
    }
  }

  public Index fetch(String id)
  throws IllegalIDException, BioException {
    int indx = Collections.binarySearch(
      indxList,
      id,
      indxList.getComparator()
    );

    if(indx < 0) {
      throw new IllegalIDException("Can't find sequence for " + id);
    }

    return (Index) indxList.get(indx);
  }

  public void store(Index indx) {
    indxList.add(indx);
  }

  public void commit()
  throws BioException {
    indxList.commit();
    try {
      // write files
      {
        PrintStream fo = new PrintStream(
          new FileOutputStream(
            new File(indexDirectory, "fileids.dat")
          )
        );
        for(int i = 0; i < fileCount; i++) {
          fo.print(i);
          fo.print('\t');
          fo.print(fileIDToFile[i]);
          fo.print('\t');
          fo.print(fileIDToFile[i].length());
          fo.println();
        }
        fo.close();
      }
    } catch (Exception e) {
      rollback();
      throw new BioException("Unable to commit. Rolled back to be safe",e);
    }
  }

  public void rollback() {
    indxList.rollback();
  }

  public Set getIDs() {
    return idSet;
  }

  public Set getFiles() {
    return new HashSet(Arrays.asList(fileIDToFile));
  }

  public SequenceFormat getFormat() {
    return format;
  }

  public SequenceBuilderFactory getSBFactory() {
    return sbFactory;
  }

  public SymbolTokenization getSymbolParser() {
    return symbolTokenization;
  }

  private interface Commitable {
    public void commit()
    throws BioException;

    public void rollback();
  }

  // records stored as:
  // seqID(\w+) \t fileID(\w+) \t start(\d+) \t length(\d+) ' ' * \n
  private abstract class FileAsList
  extends AbstractList
    implements /* RandomAccess, */ Commitable {
    private RandomAccessFile mappedFile;
    private int commitedRecords;
    private int lastIndx;
    private Object lastRec;
    private byte[] buffer;

    public FileAsList(RandomAccessFile mappedFile, int recordLength) {
      this.mappedFile = mappedFile;
      buffer = new byte[recordLength];
    }

    public Object get(int indx) {
      if(indx < 0 || indx >= size()) {
        throw new IndexOutOfBoundsException();
      }

      if(indx == lastIndx) {
        return lastRec;
      }

      long offset = indx * buffer.length;
      try {
        mappedFile.seek(offset);
        mappedFile.readFully(buffer);
      } catch (IOException ioe) {
        throw new BioError("Failed to seek for record",ioe);
      }

      lastRec = parseRecord(buffer);
      lastIndx = indx;
      return lastRec;
    }

    public int size() {
      try {
        return (int) (mappedFile.length() / (long) buffer.length);
      } catch (IOException ioe) {
        throw new BioError("Can't read file length",ioe);
      }
    }

    public boolean add(Object o) {
      generateRecord(buffer, o);

      try {
        mappedFile.seek(mappedFile.length());
        mappedFile.write(buffer);
      } catch (IOException ioe) {
        throw new BioError("Failed to write index",ioe);
      }

      return true;
    }

    public void commit() {
      Collections.sort(indxList, indxList.getComparator());
      commitedRecords = indxList.size();
    }

    public void rollback() {
      try {
        mappedFile.setLength((long) commitedRecords * (long) buffer.length);
      } catch (Throwable t) {
        throw new BioError(
          "Could not roll back. " +
          "The index store will be in an inconsistent state " +
          "and should be discarded. File: " + mappedFile, t
        );
      }
    }

    protected abstract Object parseRecord(byte[] buffer);
    protected abstract void generateRecord(byte[] buffer, Object item);
    protected abstract Comparator getComparator();
  }

  private class IndexFileAsList extends FileAsList {
    private Comparator INDEX_COMPARATOR = new Comparator() {
      public int compare(Object a, Object b) {
        String as;
        String bs;

        if(a instanceof Index) {
          as = ((Index) a).getID();
        } else {
          as = (String) a;
        }

        if(b instanceof Index) {
          bs = ((Index) b).getID();
        } else {
          bs = (String) b;
        }

        return STRING_CASE_SENSITIVE_ORDER.compare(as, bs);
      }
    };

    public IndexFileAsList(RandomAccessFile file, int recordLength) {
      super(file, recordLength);
    }

    protected Object parseRecord(byte[] buffer) {
      int lastI = 0;
      int newI = 0;
      while(buffer[newI] != '\t') {
        newI++;
      }
      String id = new String(buffer, lastI, newI);

      while(buffer[newI] != '\t') {
        newI++;
      }
      File file = getFileForID(Integer.parseInt(new String(buffer, lastI, newI).trim()));

      while(buffer[newI] != '\t') {
        newI++;
      }
      long start = Long.parseLong(new String(buffer, lastI, newI));

      int length = Integer.parseInt(
        new String(buffer, newI + 1, buffer.length)
      );

      return new SimpleIndex(file, start, length, id);
    }

    protected void generateRecord(byte[] buffer, Object item) {
      Index indx = (Index) item;

      String id = indx.getID();
      int fileID = getIDForFile(indx.getFile());
      String start = String.valueOf(indx.getStart());
      String length = String.valueOf(indx.getLength());

      int i = 0;
      byte[] str;

      str = id.getBytes();
      for(int j = 0; j < str.length; j++) {
        buffer[i++] = str[j];
      }

      buffer[i++] = '\t';

      str = String.valueOf(fileID).getBytes();
      for(int j = 0; j < str.length; j++) {
        buffer[i++] = str[j];
      }

      buffer[i++] = '\t';

      str = start.getBytes();
      for(int j = 0; j < str.length; j++) {
        buffer[i++] = str[j];
      }

      buffer[i++] = '\t';

      str = length.getBytes();
      for(int j = 0; j < str.length; j++) {
        buffer[i++] = str[j];
      }

      while(i < buffer.length - 1) {
        buffer[i++] = ' ';
      }

      buffer[i] = '\n';
    }

    public Comparator getComparator() {
      return INDEX_COMPARATOR;
    }
  }

  private class ListAsSet
  extends AbstractSet {
    public Iterator iterator() {
      return indxList.iterator();
    }

    public int size() {
      return indxList.size();
    }
  }
}
