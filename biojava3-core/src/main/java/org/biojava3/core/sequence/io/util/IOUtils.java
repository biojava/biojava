package org.biojava3.core.sequence.io.util;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.zip.GZIPInputStream;

import org.biojava3.core.exceptions.ParserException;

public class IOUtils {

  private static final int BUFFER = 4096;

  /**
   * Closes any Object which implements the interface @ link Closeable} and
   * sending any error to the logger but not forcing any explicit catching of
   * stream errors.
   *
   * @param c The stream to close
   */
  public static void close(Closeable c) {
    try {
      if (c != null) {
        c.close();
      }
    } catch (IOException e) {
      Logger log = Logger.getLogger(IOUtils.class.getName());
      log.log(Level.WARNING, "Cannot close down the given Closeable object", e);
    }
  }

  /**
   * Moves the bytes from input to output using a 4KB byte array.
   *
   * @param input Input stream of bytes
   * @param output Output stream of bytes
   * @throws IOException If anything occurs in the case of the reads and writes
   */
  public static void copy(InputStream input, OutputStream output)
      throws IOException {
    byte[] buffer = new byte[BUFFER];
    int n = 0;
    while (-1 != (n = input.read(buffer))) {
      output.write(buffer, 0, n);
    }
  }

  /**
   * Takes in a reader and a processor, reads every line from the given
   * file and then invokes the processor. What you do with the lines is
   * dependent on your processor.
   *
   * The code will automatically close the given BufferedReader.
   *
   * @param br The reader to process
   * @param processor The processor to invoke on all lines
   * @throws ParserException Can throw this if we cannot parse the given reader
   */
  public static void processReader(BufferedReader br, ReaderProcessor processor) throws ParserException {
    String line;
    try {
      while( (line = br.readLine()) != null ) {
        processor.process(line);
      }
    }
    catch(IOException e) {
      throw new ParserException("Could not read from the given BufferedReader");
    }
    finally {
      close(br);
    }
  }

  /**
   * Returns the contents of a buffered reader as a list of strings
   *
   * @param br BufferedReader to read from; <strong>will be closed</strong>
   * @return List of Strings
   * @throws ParserException Can throw this if we cannot parse the given reader
   */
  public static List<String> getList(BufferedReader br) throws ParserException {
    final List<String> list = new ArrayList<String>();
    processReader(br, new ReaderProcessor() {
      public void process(String line) {
        list.add(line);
      }
    });
    return list;
  }

  /**
   * Delegates to {@link #getList(BufferedReader)} by wrapping the InputStream
   * in a valid reader. No encoding is mentioned so if you need anything
   * more advanced then use the other version of this method.
   *
   * @param is InputStream which is a text file
   * @return List of Strings representing the lines of the files
   * @throws ParserException Can throw this if the file is not a file or we
   * cannot parse it
   */
  public static List<String> getList(InputStream is) throws ParserException {
    return getList(new BufferedReader(new InputStreamReader(is)));
  }

  /**
   * Delegates to {@link #getList(InputStream)} by wrapping the File
   * in a valid stream. No encoding is mentioned so if you need anything
   * more advanced then use the other version of this method. Since this
   * uses {@link #openFile(File)} this code can support GZipped and plain
   * files.
   *
   * @param file File which is a text file
   * @return List of Strings representing the lines of the files
   * @throws ParserException Can throw this if the file is not a file or we
   * cannot parse it
   */
  public static List<String> getList(File file) throws ParserException {
    return getList(openFile(file));
  }

  /**
   * For a filename this code will check the extension of the file for a
   * .gz extension. If it finds one then the InputStream given back
   * is a {@link GZIPInputStream}. Otherwise we return a normal
   * {@link FileInputStream}.
   *
   * @param file File which may or may not be GZipped
   * @return The final stream
   * @throws ParserException Can throw this if the file is not a file or we
   * cannot open it for processing
   */
  public static InputStream openFile(File file) throws ParserException {
    final InputStream is;
    if(!file.isFile()) {
      throw new ParserException("The file "+file+" is not a file.");
    }
    String name = file.getName();
    try {
      if(name.endsWith(".gz")) {
        is = new GZIPInputStream(new FileInputStream(file));
      }
      else {
        is = new FileInputStream(file);
      }
    }
    catch(IOException e) {
      throw new ParserException("Cannot open "+file+" for processing", e);
    }
    return is;
  }

  /**
   * Closure interface used when working with
   * {@link IOUtils#processReader(String)}. Each time a line is encountered
   * the object that implements this interface will be invoked.
   *
   * @author ayates
   */
  public static interface ReaderProcessor {
    void process(String line) throws IOException;
  }

}
