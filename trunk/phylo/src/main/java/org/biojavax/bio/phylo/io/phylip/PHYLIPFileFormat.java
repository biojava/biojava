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
package org.biojavax.bio.phylo.io.phylip;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.Reader;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.seq.io.ParseException;

/**
 * Reads PHYLIP interleaved alignment files and fires events at a PHYLIPFileListener object.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public class PHYLIPFileFormat {
  
  private static int MAX_NAME_LENGTH = 10;
  
//Prevent instances.
  private PHYLIPFileFormat() {
  }
  
  public static void parseFile(final PHYLIPFileListener listener, final File inputFile) throws IOException, ParseException {
    final FileReader fr = new FileReader(inputFile);
    try {
      PHYLIPFileFormat.parseReader(listener, fr);
    } finally {
      fr.close();
    }
  }
  
  public static void parseInputStream(final PHYLIPFileListener listener,
      final InputStream inputStream) throws IOException, ParseException {
    PHYLIPFileFormat.parseReader(listener,
        new InputStreamReader(inputStream));
  }
  
  public static void parseReader(final PHYLIPFileListener listener,
      final Reader inputReader) throws IOException, ParseException {
    PHYLIPFileFormat.parse(listener, inputReader instanceof BufferedReader ? (BufferedReader) inputReader
            : new BufferedReader(inputReader));
  }
  
  public static void parse(final PHYLIPFileListener listener, BufferedReader reader) throws IOException, ParseException {
   listener.startFile();
   List sequenceNames = new ArrayList();
   String headerline = reader.readLine();
   Pattern pattern = Pattern.compile("\\s*(\\d+)\\s+(\\d+)\\s*");
   Matcher matcher = pattern.matcher(headerline);
   if (!matcher.matches()) {
     throw new ParseException("Invalid header line.");
   }
   int sequenceCount = Integer.parseInt(matcher.group(1));
   listener.setSequenceCount(sequenceCount);
   int sitesCount = Integer.parseInt(matcher.group(2));
   listener.setSitesCount(sitesCount);
   int currentSequenceIndex = 0;
   boolean collectedAllNames = false;
   String line = reader.readLine();
   while (line != null) {
     if (line.length() == 0) continue;
     if (!collectedAllNames) {
       String name = line.substring(0, MAX_NAME_LENGTH).trim();
       sequenceNames.add(name);
       line = line.substring(MAX_NAME_LENGTH).replaceAll("\\s", "");
     }
     listener.setCurrentSequenceName((String)sequenceNames.get(currentSequenceIndex));
     listener.receiveSequence(line);
     if (sequenceNames.size() == sequenceCount) collectedAllNames = true;
     currentSequenceIndex++;
     if (!(currentSequenceIndex < sequenceCount)) currentSequenceIndex = 0;
     line = reader.readLine();
   }
   listener.endFile();
  }
  
  /**
   * Writes the given Alignment in PHYLIP format to a file.
   * 
   * @param file
   *            the file to write to.
   * @param alignment
   *            the Alignment object to write.
   * @throws IOException
   *             if there is a problem during writing.
   */
  public static void writeFile(final File file, final Alignment alignment)
      throws IOException {
    final FileWriter fw = new FileWriter(file);
    try {
      PHYLIPFileFormat.writeWriter(fw, alignment);
    } finally {
      fw.close();
    }
  }

  /**
   * Writes the given Alignment in PHYLIP format to a stream.
   * 
   * @param os
   *            the stream to write to.
   * @param alignment
   *            the Alignment object to write.
   * @throws IOException
   *             if there is a problem during writing.
   */
  public static void writeStream(final OutputStream os,
      final Alignment alignment) throws IOException {
    final OutputStreamWriter ow = new OutputStreamWriter(os);
    PHYLIPFileFormat.writeWriter(ow, alignment);
  }

  /**
   * Writes the given Alignment in PHYLIP format to a writer.
   * 
   * @param writer
   *            the writer to write to.
   * @param alignment
   *            the Alignment object to write.
   * @throws IOException
   *             if there is a problem during writing.
   */
  public static void writeWriter(final Writer writer,
      final Alignment alignment) throws IOException {
    String lineSep = System.getProperty("line.separator");
    writer.write("" + (alignment.getLabels().size()));
    writer.write("   ");
    writer.write("" + (alignment.length()) + lineSep);
    for (Iterator i = alignment.getLabels().iterator(); i.hasNext();) {
      String label = (String)i.next();
      String sequence = alignment.symbolListForLabel(label).seqString();
      writer.write(PHYLIPFileFormat.formatSequenceLabel(label));
      writer.write(sequence);
      writer.write(lineSep);
    } 
    writer.flush();
  }
  
  private static String formatSequenceLabel(String label) {
    if (label.length() > MAX_NAME_LENGTH) {
      return label.substring(0, MAX_NAME_LENGTH);
    } else if (label.length() < MAX_NAME_LENGTH) {
      StringBuffer buffer = new StringBuffer(label);
      while (buffer.length() < MAX_NAME_LENGTH) {
        buffer.append(" ");
      }
      return buffer.toString();
    }
    else {
      return label;
    }
  }
}
