package org.biojava3.core.sequence.io.util;

import java.io.BufferedReader;
import java.io.Closeable;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.List;

import org.biojava3.core.exceptions.ParserException;

public class IOUtils {

  private static final int BUFFER = 4096;

  /**
   * Closes any Object which implements the interface @ link Closeable} and
   * sending any error to the logger but not forcing any explicit catching of
   * stream errors.
   *
   * @param c
   *          The stream to close
   */
  public static void close(Closeable c) {
    try {
      if (c != null) {
        c.close();
      }
    } catch (IOException e) {
      // TODO Get logger
    }
  }

  public static void copy(InputStream input, OutputStream output)
      throws IOException {
    byte[] buffer = new byte[BUFFER];
    int n = 0;
    while (-1 != (n = input.read(buffer))) {
      output.write(buffer, 0, n);
    }
  }

  public static void processReader(BufferedReader br, ReaderProcessor processor) {
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
   */
  public static List<String> getList(BufferedReader br) {
    final List<String> list = new ArrayList<String>();
    processReader(br, new ReaderProcessor() {
      public void process(String line) {
        list.add(line);
      }
    });
    return list;
  }

  public static List<String> getList(InputStream is) {
    return getList(new BufferedReader(new InputStreamReader(is)));
  }

  public static interface ReaderProcessor {
    void process(String line);
  }

}
