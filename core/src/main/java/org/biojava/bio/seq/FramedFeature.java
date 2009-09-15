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

package org.biojava.bio.seq;

import java.io.Serializable;

/**
 * Title:         FramedFeature.<p>
 * Description:  An feature that includes the concept of frame
 * by extending stranded.<p>
 * Copyright:    Copyright (c) 2001<p>
 * @author Mark Schreiber
 * @version 1.0
 */

public interface FramedFeature extends StrandedFeature {

  public static ReadingFrame FRAME_0 = new ReadingFrame("FRAME_0",Frame.FRAME_0);
  public static ReadingFrame FRAME_1 = new ReadingFrame("FRAME_1",Frame.FRAME_1);
  public static ReadingFrame FRAME_2 = new ReadingFrame("FRAME_2",Frame.FRAME_2);

  /**
   * return the reading frame of the feature.
   */
  ReadingFrame getReadingFrame();

  /**
   * @return the Strand that the feature is found on
   */
  Strand getStrand();

  public static class Template extends StrandedFeature.Template{
    public ReadingFrame readingFrame;
  }

  /**
   * A singleton to hold the frame information
   */
  public static class ReadingFrame implements Frame, Serializable{
    private final String text;
    private final int frame;

    private ReadingFrame(String text,int frame){
      this.text = text;
      this.frame = frame;
    }

    public String toString(){
      return text;
    }

    public int getFrame(){
      return frame;
    }
  }
}
