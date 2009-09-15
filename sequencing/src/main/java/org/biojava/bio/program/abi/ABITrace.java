 
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
 
package org.biojava.bio.program.abi;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.BufferedInputStream;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

import org.biojava.bio.BioError;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.AtomicSymbol;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SymbolList;


/**
 * Title: ABITrace<br><br>
 * ABITrace is a class for managing ABI file information,
 *  it is capable of opening an ABI file and storing
 *  the most important fields, which can be recalled as simple java types. It can also return
 *  an image corresponding to the trace.
 *  It has three constructors with input types <code>File, URL, and byte[]</code>.<br><br>
 *  ABI files contain two sets of basecall and sequence data, one that was originally
 *  created programatically and the other, which is an editable copy. This version of this object
 *  only references the original unedited data.<br>
 *
 * Copyright (c) 2001
 * @author David H. Klatte, Ph.D.
 * @author Matthew Pocock
 * @version 0.5alpha
 */
public class ABITrace
{

  //the next three lines are the important persistent data
  private String sequence;
  private int A[], G[], C[], T[], Basecalls[];
  private int TraceLength, SeqLength;

  //This is the actual file data.
  private byte[] TraceData;

  private int maximum = 0;

  //the next four declaration lines comprise the file index information
  private int MacJunk=0; //sometimes when macintosh files are
                         //FTPed in binary form, they have 128 bytes
                         //of crap pre-pended to them. This constant
                         //allows ABITrace to handle that in a way that
                         //is invisible to the user.
  private static int AbsIndexBase=26; //The file location of the Index pointer
  private int IndexBase,  PLOC;

  //the next declaration is for the actual file pointers
  private  int DATA9, DATA10, DATA11, DATA12, PBAS2, FWO;

/**
 * The File constructor opens a local ABI file and parses the content.
 * @param ABIFile is a <code>java.io.File</code> on the local file system.
 * @throws IOException if there is a problem reading the file.
 * @throws IllegalArgumentException if the file is not a valid ABI file.
 */
  public ABITrace( File ABIFile ) throws IOException
  {
    byte[] bytes = null;
    ByteArrayOutputStream baos = new ByteArrayOutputStream();
    FileInputStream fis = new FileInputStream(ABIFile);
    BufferedInputStream bis = new BufferedInputStream(fis);
    int b;
    while ((b = bis.read()) >= 0)
    {
      baos.write(b);
    }
    bis.close(); fis.close(); baos.close();
    bytes = baos.toByteArray();
    initData(bytes);
  }

/**
 * The URL constructor opens an ABI file from any URL.
 * @param ABIFile is a <code>java.net.URL</code> for an ABI trace file.
 * @throws IOException if there is a problem reading from the URL.
 * @throws IllegalArgumentException if the URL does not contain a valid ABI file.
 */
  public ABITrace( URL ABIFile ) throws IOException
  {
    byte[] bytes = null;
    ByteArrayOutputStream baos = new ByteArrayOutputStream();
    InputStream is = ABIFile.openStream();
    BufferedInputStream bis = new BufferedInputStream(is);
    int b;
    while ((b = bis.read()) >= 0)
    {
      baos.write(b);
    }
    bis.close(); is.close(); baos.close();
    bytes = baos.toByteArray();
    initData(bytes);
  }

/**
 * The <code>byte[]</code> constructor parses an ABI file represented as a byte array.
 * @throws IllegalArgumentException if the data does not represent a valid ABI file.
 */
  public ABITrace(byte[] ABIFileData)
  {
    initData(ABIFileData);
  }

/**
 * Returns the length of the sequence (number of bases) in this trace.
 */
  public int getSequenceLength() { return SeqLength; }

/**
 * Returns the length of the trace (number of x-coordinate points in the graph).
 */
  public int getTraceLength() { return TraceLength; }

/**
 * Returns an <code>int[]</code> array that represents the basecalls - each int in the
 * array corresponds to an x-coordinate point in the graph that is a peak (a base location).
 */
  public int[] getBasecalls() { return Basecalls; }

/**
 * Returns the original programatically determined (unedited) sequence as a <code>SymbolList</code>.
 */
  public SymbolList getSequence() throws BioError
  { 
    try {
      return DNATools.createDNA(sequence); 
    }
    catch (IllegalSymbolException ise) {
      // this should be impossible!
      throw new BioError(ise);
    }
  }

/**
 * Returns one of the four traces - all of the y-coordinate values,
 * each of which correspond to a single x-coordinate relative to the
 * position in the array, so that if element 4 in the array is 972, then
 * x is 4 and y is 972 for that point.
 *
 * @param base  the DNA AttomicSymbol to retrieve the trace values for
 * @return an array of ints giving the entire trace for that base
 * @throws IllegalSymbolException if the base is not valid
 */
  public int[] getTrace (AtomicSymbol base) throws IllegalSymbolException
  {
    if (base == DNATools.a()) {
      return A;
    } else if (base == DNATools.c()) {
      return C;
    } else if (base == DNATools.g()) {
      return G;
    } else if (base == DNATools.t()) {
      return T;
    } else {
      DNATools.getDNA().validate(base);
      throw new IllegalSymbolException("Don't know symbol: " + base);
    }
  }

/**
 * Returns a BufferedImage that represents the entire trace. The height can be set precisely in
 * pixels, the width in pixels is determined by the scaling factor times the number
 * of points in the trace (<code>getTraceLength()</code>). The entire trace is represented
 * in the returned image.
 *
 * @param imageHeight is the desired height of the image in pixels.
 * @param widthScale indiates how many horizontal pixels to use to represent a single x-coordinate (try 2).
 */
  public BufferedImage getImage(int imageHeight, int widthScale)
  {
    BufferedImage out = new BufferedImage(TraceLength * widthScale, imageHeight, BufferedImage.TYPE_BYTE_INDEXED);
    Graphics2D g = out.createGraphics();
    Color acolor = Color.green.darker();
    Color ccolor = Color.blue;
    Color gcolor = Color.black;
    Color tcolor = Color.red;
    Color ncolor = Color.pink;
    double scale = calculateScale(imageHeight);
    int[] bc = Basecalls;
    char[] seq = sequence.toCharArray();
    g.setBackground(Color.white);
    g.clearRect(0, 0, TraceLength * widthScale, imageHeight);
    int here = 0;
    int basenum = 0;
    for (int q = 1; q <= 5; q++)
    {
      for (int x = 0; x <= TraceLength - 2; x++)
      {
        if (q==1)
        {
          g.setColor(acolor);
          g.drawLine(2*x, transmute(A[x], imageHeight, scale),
                     2*(x + 1), transmute(A[x+1], imageHeight, scale));
        }
        if (q==2)
        {
          g.setColor(ccolor);
          g.drawLine(2*x, transmute(C[x], imageHeight, scale),
                     2*(x + 1), transmute(C[x+1], imageHeight, scale));
        }
        if (q==3)
        {
          g.setColor(tcolor);
          g.drawLine(2*x, transmute(T[x], imageHeight, scale),
                     2*(x + 1), transmute(T[x+1], imageHeight, scale));
        }
        if (q==4)
        {
          g.setColor(gcolor);
          g.drawLine(2*x, transmute(G[x], imageHeight, scale),
                     2*(x + 1), transmute(G[x+1], imageHeight, scale));
        }
        if (q==5)
        {
          if ((here > bc.length-1) || (basenum > seq.length-1)) break;
          if (bc[here] == x)
          {
            g.drawLine(2*x, transmute(-2, imageHeight, 1.0),
                       2*x, transmute(-7, imageHeight, 1.0));
            if ((basenum+1)%10 == 0) //if the basecount is divisible by ten
                            //add a number
            {
              g.drawLine(2*x, transmute(-20, imageHeight, 1.0),
                         2*x, transmute(-25, imageHeight, 1.0));
              g.drawString(Integer.toString(basenum+1),
                           2*x-3, transmute(-36, imageHeight, 1.0));
            }
            switch (seq[basenum])
            {
              case 'A': case 'a': g.setColor(acolor); break;
              case 'C': case 'c': g.setColor(ccolor); break;
              case 'G': case 'g': g.setColor(gcolor); break;
              case 'T': case 't': g.setColor(tcolor); break;
              default: g.setColor(ncolor);
            }
            g.drawChars(seq, basenum, 1,
                    2*x-3, transmute(-18, imageHeight, 1.0));
            g.setColor(Color.black);
            here++; basenum++;
          }
        }
      }
    }
    return out;
  }

/**
 * Initialize all of the data fields for this object.
 * @throws IllegalArgumentException which will propagate to all of the constructors.
 */
  private void initData(byte[] fileData)
  {
    TraceData = fileData;
    if (isABI())
    {
      setIndex();
      setBasecalls();
      setSeq();
      setTraces();
    }
    else throw new IllegalArgumentException("Not a valid ABI file.");
  }

/**
 * A utility method which fills array b with data from the trace starting at traceDataOffset.
 */
  private void getSubArray(byte[] b, int traceDataOffset)
  {
    for (int x=0; x<=b.length-1; x++)
    {
      b[x] = TraceData[traceDataOffset + x];
    }
  }

/**
 * Shuffle the pointers to point to the proper spots in the trace, then load the
 * traces into their arrays.
 */
  private void setTraces()
  {
    int pointers[] = new int[4]; //alphabetical, 0=A, 1=C, 2=G, 3=T
    int datas[] = new int[4];
    char order[] = new char[4];

    datas[0] = DATA9;
    datas[1] = DATA10;
    datas[2] = DATA11;
    datas[3] = DATA12;

    for (int i=0; i<=3; i++)
    {
      order[i]=(char) TraceData[FWO+i];
    }

    for (int i=0; i <=3; i++)
    {
      switch (order[i])
      {
        case 'A': case 'a':
          pointers[0] = datas[i];
          break;
        case 'C': case 'c':
          pointers[1] = datas[i];
          break;
        case 'G': case 'g':
          pointers[2] = datas[i];
          break;
        case 'T': case 't':
          pointers[3] = datas[i];
          break;
        default:
          throw new IllegalArgumentException("Trace contains illegal values.");
      }
    }

    A = new int[TraceLength];
    C = new int[TraceLength];
    G = new int[TraceLength];
    T = new int[TraceLength];

    for (int i=0; i <=3; i++)
    {
      byte[] qq = new byte[TraceLength*2];
      getSubArray(qq, pointers[i]);
      DataInputStream dis = new DataInputStream(new ByteArrayInputStream(qq));
      for (int x=0; x <=TraceLength - 1; x++)
      {
        try
        {
          if (i == 0) A[x] = (int) dis.readShort();
          if (i == 1) C[x] = (int) dis.readShort();
          if (i == 2) G[x] = (int) dis.readShort();
          if (i == 3) T[x] = (int) dis.readShort();
        }catch(IOException e)//This shouldn't happen. If it does something must be seriously wrong.
        {
          throw new IllegalStateException("Unexpected IOException encountered while manipulating internal streams.");
        }
      }
    }
    return;
  }

/**
 * Fetch the sequence from the trace data.
 */
  private void setSeq()
  {
    char tempseq[] = new char[SeqLength];
    for (int x = 0; x <= SeqLength - 1; ++x)
    {
      tempseq[x] = (char) TraceData[PBAS2 + x];
    }
    sequence = new String (tempseq);
  }


/**
 * Fetch the basecalls from the trace data.
 */
  private void setBasecalls()
  {
    Basecalls = new int[SeqLength];
    byte[] qq = new byte[SeqLength*2];
    getSubArray(qq, PLOC);
    DataInputStream dis = new DataInputStream(new ByteArrayInputStream(qq));
    for (int i = 0; i <= SeqLength -1; ++i)
    {
      try
      {
        Basecalls[i]=(int) dis.readShort();
      }catch(IOException e)//This shouldn't happen. If it does something must be seriously wrong.
      {
        throw new IllegalStateException("Unexpected IOException encountered while manipulating internal streams.");
      }
    }
  }

/**
 * Utility method to return an int beginning at <code>pointer</code> in the TraceData array.
 */
  private int getIntAt(int pointer)
  {
    int out = 0;
    byte[] temp = new byte[4];
    getSubArray(temp, pointer);
    try
    {
      DataInputStream dis = new DataInputStream(new ByteArrayInputStream(temp));
      out = dis.readInt();
    }catch(IOException e) //This shouldn't happen. If it does something must be seriously wrong.
    {
      throw new IllegalStateException("Unexpected IOException encountered while manipulating internal streams.");
    }
    return out;
  }

/**
 * Utility method to translate y coordinates from graph space (where up is greater)
 * to image space (where down is greater).
 */
  private int transmute(int ya, int height, double scale)
  {
    return (height - 45 - (int) (ya * scale));
  }

/**
 * Get the maximum height of any of the traces. The data is persisted for performance
 * in the event of multiple calls, but it initialized lazily.
 */
  private int getMaximum()
  {
    if (maximum > 0) return maximum;
    int max = 0;
    for (int x=0; x<=T.length-1; x++)
    {
      if (T[x] > max) max = T[x];
      if (A[x] > max) max = A[x];
      if (C[x] > max) max = C[x];
      if (G[x] > max) max = G[x];
    }
    return max;
  }

  //calculates the necessary scaling to allow the trace to fit vertically
  //in the space specified.
/**
 * Returns the scaling factor necessary to allow all of the traces to fit vertically
 * into the specified space.
 * @param <code>height</code> - the required height in pixels.
 */
  private double calculateScale(int height)
  {
    double newScale = 0.0;
    double max = (double)getMaximum();
    double ht = (double)height;
    newScale = ((ht - 50.0))/max;
    return newScale;
  }

/**
 * Sets up all of the initial pointers to the important records in TraceData.
 */
  private void setIndex()
  {
    int DataCounter, PBASCounter, PLOCCounter, NumRecords;
    byte[] RecNameArray = new byte[4];
    String RecName;

    DataCounter = 0; PBASCounter = 0; PLOCCounter = 0;

    IndexBase = getIntAt(AbsIndexBase + MacJunk);
    NumRecords = getIntAt(AbsIndexBase - 8 + MacJunk);

    for (int record = 0; record <= NumRecords - 1; record++)
    {
      getSubArray(RecNameArray, (IndexBase + (record * 28)));
      RecName = new String (RecNameArray);
      if (RecName.equals("FWO_"))
        FWO = IndexBase + (record * 28) + 20;
      if (RecName.equals("DATA"))
      {
        ++DataCounter;
        if (DataCounter == 9)
          DATA9 = IndexBase + (record * 28) + 20;
        if (DataCounter == 10)
          DATA10 = IndexBase + (record * 28) + 20;
        if (DataCounter == 11)
          DATA11 = IndexBase + (record * 28) + 20;
        if (DataCounter == 12)
          DATA12 = IndexBase + (record * 28) + 20;
      }
      if (RecName.equals("PBAS"))
      {
        ++PBASCounter;
        if (PBASCounter == 2)
          PBAS2 = IndexBase + (record * 28) + 20;
      }
      if (RecName.equals("PLOC"))
      {
        ++PLOCCounter;
        if (PLOCCounter == 2)
          PLOC = IndexBase + (record * 28) + 20;
      }

    } //next record
    TraceLength = getIntAt(DATA12 - 8);
    SeqLength = getIntAt(PBAS2-4);
    PLOC = getIntAt(PLOC) + MacJunk;
    DATA9 = getIntAt(DATA9) + MacJunk;
    DATA10 = getIntAt(DATA10) + MacJunk;
    DATA11 = getIntAt(DATA11) + MacJunk;
    DATA12 = getIntAt(DATA12) + MacJunk;
    PBAS2 = getIntAt(PBAS2) + MacJunk;
  }

/**
 * Test to see if the file is ABI format by checking to see that the first three bytes
 * are "ABI". Also handle the special case where 128 bytes were prepended to the file
 * due to binary FTP from an older macintosh system.
 */
  private boolean isABI()
  {
    char ABI[] = new char[4];

    for (int i=0; i<=2; i++)
    {
      ABI[i]=(char) TraceData[i];
    }
    if (ABI[0] == 'A' && (ABI[1] == 'B' && ABI[2] == 'I'))
    {
      return true;
    }
    else
    {
      for (int i=128; i<=130; i++)
      {
        ABI[i]=(char) TraceData[i];
      }
      if (ABI[0] == 'A' && (ABI[1] == 'B' && ABI[2] == 'I'))
      {
        MacJunk=128;
        return true;
      }
      else
        return false;
    }
  }
}

